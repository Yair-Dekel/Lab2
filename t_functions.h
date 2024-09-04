#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#define GA_MAXITER		 50			// maximum iterations
#define GA_POPSIZE		 100		    // ga population size
#define GA_ELITRATE		 0.05f		    // elitism rate
#define GA_ELITRATE		 0.05f		    // elitism rate
#define GA_MUTATIONRATE	 0.4f		    // mutation rate
#define GA_MUTATION		 RAND_MAX * GA_MUTATIONRATE
#define GA_TARGET		 std::string("Hello world!")
#define MAX_AGE          3
#define REPRODUCE_THRESH 1


double Non_Unform_Mutation(int generation) {
    double r = -0.01;
    double Pmax = GA_MUTATIONRATE;
    return (2 * Pmax * exp(r * generation)) / (Pmax + Pmax * exp(r * generation));
    
}



template <typename T>
T single_point_crossover(T &parent1, T &parent2)
{
    int cross_point = rand() % parent1.size();
    T child;
    for (int i = 0; i < cross_point; i++){
        child.push_back(parent1[i]);
    }
    for (int i = cross_point; i < parent2.size(); i++){
        child.push_back(parent2[i]);
    }
    return child;
}

template <typename T>
T two_point_crossover(T &parent1, T &parent2)
{
    int cross_point1 = rand() % parent1.size();
    int cross_point2 = rand() % parent1.size();
    T child;
    for (int i = 0; i < cross_point1; i++){
        child.push_back(parent1[i]);
    }
    for (int i = cross_point1; i < cross_point2; i++){
        child.push_back(parent2[i]);
    }
    for (int i = cross_point2; i < parent2.size(); i++){
        child.push_back(parent1[i]);
    }
    return child;
}

template <typename T>
T uniform_crossover(T &parent1, T &parent2)
{
    T child;
    for (int i = 0; i < parent1.size(); i++){
        int j = rand() % 2;
        if (j == 0){
            child.push_back(parent1[i]);
        }
        else{
            child.push_back(parent2[i]);
        }
    }
    return child;
}

template <typename T>
void top_average_and_variance(T &pop, std::vector<double> &average_vec, std::vector<double> &variance_vec)
{
    int best = pop.get_best().get_fitness();
    int sum = pop.calc_fitness_sum();
     
    double average = sum / GA_POPSIZE;
    std::cout << "Top-average: " << best / average << std::endl;
    average_vec.push_back(best / average);

    double variance = pop.calc_fitness_variance(average) / GA_POPSIZE;
    std::cout << "Fitness variance: " << variance << std::endl;
    variance_vec.push_back(variance);

}

template <typename T>
void elitism(T &pop, T &buffer, int esize)
{
    for (int i = 0; i < esize; i++){
        buffer.set(i, pop.get(i));
        if(pop.get_flags().age_F)buffer.increase_age(i);
    }
}

template <typename T>
void mate(T &pop, T &buffer)
{
    int esize = GA_POPSIZE * GA_ELITRATE;
    
    std::vector<int> parents;
    if(pop.get_flags().get_elitism_F()) elitism(pop, buffer, esize);
    else esize = 0;

    if(pop.get_flags().get_SUS_F()){
        pop.parent_selection_SUS(parents, GA_POPSIZE - esize);

        for(int i=esize; i<GA_POPSIZE; i++){
            int parent2 = parents[rand() % parents.size()];
            buffer.set(i, pop.mate(i, parent2)); 
        }
    }
    else{
        for (int i = esize; i < GA_POPSIZE; i++){
            int parent1;
            int parent2;
            if(pop.get_flags().get_RWS_F()){
                pop.parent_selection_RWS(parent1);
                pop.parent_selection_RWS(parent2);
            }
            else if(pop.get_flags().get_tournament_F()){
                int k = rand() % 10 + 2;
                double p = 0.75;
                pop.tournament_selection(k,p,parent1);
                pop.tournament_selection(k,p,parent2);
            }
            else{
                parent1 = rand() % (GA_POPSIZE / 2);
                parent2 = rand() % (GA_POPSIZE / 2);   
            }
            buffer.set(i, pop.mate(parent1, parent2));
        }
    }
}

template <typename T, typename C>
void non_deternistic_crowding(T &pop, T &buffer, const C &individual)
{
    int esize = GA_POPSIZE * GA_ELITRATE;
    elitism(pop, buffer, esize);
    for (int i = esize; i < GA_POPSIZE; i++){
        
        int parent1;
        int parent2;
        pop.parent_selection_RWS(parent1);
        pop.parent_selection_RWS(parent2);

        C parent = pop.get(i);
        C child = pop.mate(parent1, parent2);
        if(boltzmann_replacement(parent, child)){
            buffer.set(i, child);
        }
        else{
            i--;
        }
    }
}

template <typename C>
bool boltzmann_replacement(const C &parent1, const C &offspring) {

    double T = 1.0;
    double expo = exp((-std::abs(static_cast<double>(parent1.get_fitness()) - static_cast<double>(offspring.get_fitness()))) / T);
    double p = expo / (1 + expo);
    double r = static_cast<double>(rand()) / RAND_MAX;
    return r < p;
}
template <typename T>
void swap(T *&population, T *&buffer)
{
    T *temp = population;
    population = buffer;
    buffer = temp;
}

template <typename T>
void increase_age(T &pop)
{
    for (int i = 0; i < GA_POPSIZE; i++){
        pop.increase_age(i);
    }
}


class flags
{    
    public:
    bool RWS_F, SUS_F, PMX_F, CX_F, inversion_F, scramble_F, simple_F, elitism_F, age_F;
    bool linear_scaling_F, tournament_F, mutation_F;
    bool best_fit_F, worst_fit_F, random_F, greedy_F, shuffle_F, fit_F;
    bool single_point_F, two_point_F, uniform_F, rand_crossover_F;

    bool niching_F, crowding_F, clustering_F, simple_mate_F;

    bool learning_F;
    flags();
    flags(const flags &f);
    flags& operator=(const flags &f);
    bool get_elitism_F() const { return elitism_F; }
    bool get_SUS_F() const { return SUS_F; }
    bool get_RWS_F() const { return RWS_F; }
    bool get_tournament_F() const { return tournament_F; }
    
};
flags::flags()
{
    RWS_F = false;
    SUS_F = false;
    PMX_F = false;
    CX_F = false;
    inversion_F = false;
    scramble_F = false;
    simple_F = false;
    elitism_F = false;
    age_F = false;
    linear_scaling_F = false;
    tournament_F = false;
    mutation_F = false;
    best_fit_F = false;
    worst_fit_F = false;
    random_F = false;
    greedy_F = false;
    shuffle_F = true;
    fit_F = true;
    single_point_F = false;
    two_point_F = false;
    uniform_F = false;
    rand_crossover_F = true;
    
    niching_F = true;
    clustering_F = false;
    crowding_F = false;
    simple_mate_F = false;

    learning_F = true;
}

flags::flags(const flags &f)
{
    this->RWS_F = f.RWS_F;
    this->SUS_F = f.SUS_F;
    this->PMX_F = f.PMX_F;
    this->CX_F = f.CX_F;
    this->inversion_F = f.inversion_F;
    this->scramble_F = f.scramble_F;
    this->simple_F = f.simple_F;
    this->elitism_F = f.elitism_F;
    this->age_F = f.age_F;
    this->linear_scaling_F = f.linear_scaling_F;
    this->tournament_F = f.tournament_F;
    this->mutation_F = f.mutation_F;
    this->best_fit_F = f.best_fit_F;
    this->worst_fit_F = f.worst_fit_F;
    this->random_F = f.random_F;
    this->greedy_F = f.greedy_F;
    this->shuffle_F = f.shuffle_F;
    this->fit_F = f.fit_F;
    this->single_point_F = f.single_point_F;
    this->two_point_F = f.two_point_F;
    this->uniform_F = f.uniform_F;
    this->rand_crossover_F = f.rand_crossover_F;

    this->niching_F = f.niching_F;
    this->clustering_F = f.clustering_F;
    this->crowding_F = f.crowding_F;
    this->simple_mate_F = f.simple_mate_F;

    this->learning_F = f.learning_F;
}

flags& flags::operator=(const flags &f)
{
    if (this == &f){
        return *this;
    }
    this->RWS_F = f.RWS_F;
    this->SUS_F = f.SUS_F;
    this->PMX_F = f.PMX_F;
    this->CX_F = f.CX_F;
    this->inversion_F = f.inversion_F;
    this->scramble_F = f.scramble_F;
    this->simple_F = f.simple_F;
    this->elitism_F = f.elitism_F;
    this->age_F = f.age_F;
    this->linear_scaling_F = f.linear_scaling_F;
    this->tournament_F = f.tournament_F;
    this->mutation_F = f.mutation_F;
    this->best_fit_F = f.best_fit_F;
    this->worst_fit_F = f.worst_fit_F;
    this->random_F = f.random_F;
    this->greedy_F = f.greedy_F;
    this->shuffle_F = f.shuffle_F;
    this->fit_F = f.fit_F;
    this->single_point_F = f.single_point_F;
    this->two_point_F = f.two_point_F;
    this->uniform_F = f.uniform_F;
    this->rand_crossover_F = f.rand_crossover_F;

    this->niching_F = f.niching_F;
    this->clustering_F = f.clustering_F;
    this->crowding_F = f.crowding_F;
    this->simple_mate_F = f.simple_mate_F;

    this->learning_F = f.learning_F;
    return *this;
}

template <typename T>
void write_numbers_to_file(const std::string &filename, const std::vector<T> &numbers) {
    std::ofstream outfile(filename);

    if (!outfile) {
        std::cerr << "Error opening file for writing" << std::endl;
        return;
    }

    for (T number : numbers) {
        outfile << number << std::endl;
    }

    outfile.close();
    std::cout << "Numbers written to " << filename << std::endl;
}
