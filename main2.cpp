#pragma warning(disable:4786)		// disable debug warning

#include <iostream>					// for cout etc.
#include <vector>					// for vector class
#include <string>					// for string class
#include <algorithm>				// for sort algorithm
#include <time.h>					// for random seed
#include <math.h>					// for abs()
#include <chrono>
#include <string.h>
#include <cstring>
#include <random>
#include <set>
#include <fstream>
#include <sstream>
#include <limits> 
#include <cmath>  
#include <map>



#include "t_functions.h"

/*#include <QtCharts>
#include <QChartView>
#include <QLineSeries>
#include <QApplication>*/


static int mutation_rate_flag_bin;
static int generation_number_bin;
static int alils_number_bin;
static int global_counter_bin=0;
static int prev_best_fitness_bin = 1000000;
static int curr_best_fitness_bin = 1000000;  
static double SF_GLOB_bin; 


#define BIN_MAX_CAPACITY 100

using namespace std;				// polluting global namespace, but hey...


class bin{
    public:
        int max_capacity;
        int remain_capacity;
        std::vector<int> items;
        bin();
        bin(int max_capacity);
        bin(const bin &b);
        bin& operator=(const bin &b);
        ~bin();
};

bin::bin(){
    max_capacity = BIN_MAX_CAPACITY;
    remain_capacity = max_capacity;
}

bin::bin(int max_capacity){
    this->max_capacity = max_capacity;
    remain_capacity = max_capacity;
}

bin::bin(const bin &b){
    items = b.items;
    max_capacity = b.max_capacity;
    remain_capacity = b.remain_capacity;
}

bin& bin::operator=(const bin &b){
    if (this == &b)
        return *this;
    items = b.items;
    max_capacity = b.max_capacity;
    remain_capacity = b.remain_capacity;
    return *this;
}

bin::~bin(){
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class bins{
    private:
        std::vector<bin> bins_vec;
        int max_capacity;
        int num_items;
        int fitness;
        int age;
        double fitness_share;
    public:
        bins();
        bins(int max_capacity);
        bins(const bins &b);
        bins& operator=(const bins &b);
        ~bins();
        std::vector<bin> get_bins() const;
        void set_bins(std::vector<bin> bins);
        int get_max_capacity() const;
        void set_max_capacity(int max_capacity);
        int get_num_items() const;
        void set_num_items(int num_items);
        int get_fitness() const;
        void set_fitness(int fitness);
        int get_age() const;
        void set_age(int age);
        void init_bins(std::vector<int> items, flags flag);
        void greedy_init(std::vector<int> items);
        void random_init(std::vector<int> items);
        void best_fit_init(std::vector<int> items);
        void worst_fit_init(std::vector<int> items);
        void print_bins() const;
        void calc_fitness();
        bins single_point_crossover(const bins &parent2, std::multiset<int> items_set);
        void push_rest(std::multiset<int> items_set);
        void update_remain_capacity();
        std::vector<int> reduction_to_vector();
        void repair_vector(std::vector<int> items);

        void set_fitness_share(double fitness_share){
            this->fitness_share = fitness_share;
        }
        double get_fitness_share(){
            return fitness_share;
        }


};

bins::bins(){
    max_capacity = BIN_MAX_CAPACITY;
    num_items = 0;
    fitness = ~0;
    age = 0;
}

bins::bins(int max_capacity){
    this->max_capacity = max_capacity;
    num_items = 0;
    fitness = ~0;
    age = 0;
}

bins::bins(const bins &b){
    bins_vec = b.bins_vec;
    max_capacity = b.max_capacity;
    num_items = b.num_items;
    fitness = b.fitness;
    age = b.age;
}

bins& bins::operator=(const bins &b){
    if (this == &b)
        return *this;
    bins_vec = b.bins_vec;
    max_capacity = b.max_capacity;
    num_items = b.num_items;
    fitness = b.fitness;
    age = b.age;
    return *this;
}

bins::~bins(){
}

std::vector<bin> bins::get_bins() const{
    return bins_vec;
}

void bins::set_bins(std::vector<bin> bins){
    bins_vec = bins;
}

int bins::get_max_capacity() const{
    return max_capacity;
}

void bins::set_max_capacity(int max_capacity){
    this->max_capacity = max_capacity;
}

int bins::get_num_items() const{
    return num_items;
}

void bins::set_num_items(int num_items){
    this->num_items = num_items;
}

int bins::get_fitness() const{
    return fitness;
}

void bins::set_fitness(int fitness){
    this->fitness = fitness;
}

int bins::get_age() const{
    return age;
}

void bins::set_age(int age){
    this->age = age;
}

void bins::init_bins(std::vector<int> items, flags flag){
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_int_distribution<> dist(1, 100);

 /*   std::cout << "flag.shuffle_F: " << flag.shuffle_F << std::endl;
    std::cout << "flag.random_F: " << flag.random_F << std::endl;
    std::cout << "flag.greedy_F: " << flag.greedy_F << std::endl;
    std::cout << "flag.best_fit_F: " << flag.best_fit_F << std::endl;
    std::cout << "flag.worst_fit_F: " << flag.worst_fit_F << std::endl;
    std::cout << "flag.fit_F: " << flag.fit_F << std::endl;
*/

    // Shuffle the vector
    if(flag.shuffle_F)
        std::shuffle(items.begin(), items.end(), g);

    int i = dist(g);

    if(flag.random_F)
        i = 1;
    else if(flag.greedy_F)
        i = 25;
    else if(flag.best_fit_F)
        i = 50;
    else if(flag.worst_fit_F)
        i = 75;
    
    if(flag.fit_F)
        i = dist(g);

    //std::cout << "i: " << i << std::endl;

    if (i > 0 && i < 25){
        //std::cout << "in random!!!!" << std::endl;
        this->random_init(items);}
    else if(i >= 25 && i < 50){
        //std::cout << "in greedy!!!!" << std::endl;
        this->greedy_init(items);}
    else if(i >= 50 && i < 75){
        //std::cout << "in best_fit!!!!" << std::endl;
        this->best_fit_init(items);}
    else if(i >= 75 && i < 100){
        //std::cout << "in worst_fit!!!!" << std::endl;
        this->worst_fit_init(items);}
    else if(i == 100){
        
        //std::cout << "in sort!!!!";
        int j = dist(g) % 2;
        if(j == 0)
            std::sort(items.begin(), items.end());
        else
            std::sort(items.begin(), items.end(), std::greater<int>());
        j= dist(g) % 4;
        if(j == 0)
            this->random_init(items);
        else if(j == 1)
            this->greedy_init(items);
        else if(j == 2)
            this->best_fit_init(items);
        else if(j == 3)
            this->worst_fit_init(items);
    }

    //this->print_bins();

    this->calc_fitness();
}

void bins::greedy_init(std::vector<int> items){
    bins_vec.clear();
    bin b(max_capacity);
    bins_vec.push_back(b);
    for (int i = 0; i < items.size(); i++){
        for(int j = 0; j < bins_vec.size(); j++){
            if (bins_vec[j].remain_capacity >= items[i]){
                bins_vec[j].items.push_back(items[i]);
                bins_vec[j].remain_capacity -= items[i];
                break;
            }
            else if (j == bins_vec.size()-1){
                bin b(max_capacity);
                bins_vec.push_back(b);
                bins_vec.back().items.push_back(items[i]);
                bins_vec.back().remain_capacity -= items[i];
                break;
            }
        }
    }
    num_items = items.size();
    this->calc_fitness();
    //std::cout << "in greedy - Fitness: " << fitness << std::endl;
}

void bins::random_init(std::vector<int> items){

    bins_vec.clear();
    bin b(max_capacity);
    bins_vec.push_back(b);
    for (int i = 0; i < items.size(); i++){
        int j = rand() % bins_vec.size();
        int count = bins_vec.size();
        while (bins_vec[j].remain_capacity < items[i] && count--){
            j = rand() % bins_vec.size();
        }
        if (bins_vec[j].remain_capacity >= items[i]){
            bins_vec[j].items.push_back(items[i]);
            bins_vec[j].remain_capacity -= items[i];
        }
        else{
            bin b(max_capacity);
            bins_vec.push_back(b);
            bins_vec.back().items.push_back(items[i]);
            bins_vec.back().remain_capacity -= items[i];
        }
    }
    num_items = items.size();
    this->calc_fitness();
    //std::cout << "in random - Fitness: " << fitness << std::endl;
}

void bins::best_fit_init(std::vector<int> items){
    bins_vec.clear();
    bin b(max_capacity);
    bins_vec.push_back(b);
    for (int i = 0; i < items.size(); i++){
        int best_fit_indx = 0;
        // Find the bin with the less remaining capacity
        for(int j = 0; j < bins_vec.size(); j++){
            if(bins_vec[j].remain_capacity < bins_vec[best_fit_indx].remain_capacity && bins_vec[j].remain_capacity >= items[i]){
                best_fit_indx = j;
            }
        }
        if (bins_vec[best_fit_indx].remain_capacity >= items[i]){ 
            bins_vec[best_fit_indx].items.push_back(items[i]);
            bins_vec[best_fit_indx].remain_capacity -= items[i];
        }
        else {
            bin b(max_capacity);
            bins_vec.push_back(b);
            bins_vec.back().items.push_back(items[i]);
            bins_vec.back().remain_capacity -= items[i];
        }
    }
    num_items = items.size();
    this->calc_fitness();
    //std::cout << "in best fit - Fitness: " << fitness << std::endl;

}

void bins::worst_fit_init(std::vector<int> items){
    bins_vec.clear();
    bin b(max_capacity);
    bins_vec.push_back(b);
    for (int i = 0; i < items.size(); i++){
        int worst_fit_indx = 0;
        for(int j = 0; j < bins_vec.size(); j++){
            if(bins_vec[j].remain_capacity > bins_vec[worst_fit_indx].remain_capacity){
                worst_fit_indx = j;
            }
        }
        if (bins_vec[worst_fit_indx].remain_capacity >= items[i]){ 
            bins_vec[worst_fit_indx].items.push_back(items[i]);
            bins_vec[worst_fit_indx].remain_capacity -= items[i];
        }
        else {
            bin b(max_capacity);
            bins_vec.push_back(b);
            bins_vec.back().items.push_back(items[i]);
            bins_vec.back().remain_capacity -= items[i];
        }
    }
    num_items = items.size();
    this->calc_fitness();
    //std::cout << "in worst fit - Fitness: " << fitness << std::endl;
}               

void bins::print_bins() const{
    for (int i = 0; i < bins_vec.size(); i++){
        std::cout << "Bin " << i+1 << ": ";
        for (int j = 0; j < bins_vec[i].items.size(); j++){
            std::cout << bins_vec[i].items[j] << " ";
        }
        std::cout << "||";
    }
    std::cout << std::endl << "Fitness: " << fitness << std::endl;
}

void bins::calc_fitness(){
    fitness = 0;
    fitness += age*10;
    for (int i = 0; i < bins_vec.size(); i++){
        int sum = 0;
        for (int j = 0; j < bins_vec[i].items.size(); j++){
            sum += bins_vec[i].items[j];
        }
        if (sum < max_capacity) fitness++;
            
        fitness += abs(sum - max_capacity);
    }
}

bins bins::single_point_crossover(const bins &parent2, std::multiset<int> items_set){
    int cross_point = rand() % bins_vec.size();
    int cross_point2 = rand() % parent2.bins_vec.size();
    auto start_erase = bins_vec.begin() + cross_point;
    
    this->bins_vec.erase(start_erase, this->bins_vec.end());


    for (int i = cross_point; i < parent2.bins_vec.size(); i++){
        this->bins_vec.push_back(parent2.bins_vec[i]);
    }



    for (int i = 0; i < this->bins_vec.size(); i++){
        for (int j = 0; j < this->bins_vec[i].items.size(); j++){
            
            size_t is_erase = items_set.erase(this->bins_vec[i].items[j]);
            
            if (is_erase == 0){
                this->bins_vec[i].items.erase(this->bins_vec[i].items.begin() + j);
                j--;

            }
        }
        
        if(this->bins_vec[i].items.size() == 0){
            this->bins_vec.erase(this->bins_vec.begin() + i);
            i--;
        }
    }

    this->update_remain_capacity();

    this->push_rest(items_set);

    this->calc_fitness();
    return *this;
}

void bins::push_rest(std::multiset<int> items_set){
    for(auto item : items_set){
        for(int j = 0; j < bins_vec.size(); j++){
            if (bins_vec[j].remain_capacity >= item){
                bins_vec[j].items.push_back(item);
                bins_vec[j].remain_capacity -= item;
                break;
            }
            else if (j == bins_vec.size()-1){
                bin b(max_capacity);
                bins_vec.push_back(b);
                bins_vec.back().items.push_back(item);
                bins_vec.back().remain_capacity -= item;
                break;
            }
        }
    }
}

void bins::update_remain_capacity(){
    for (int i = 0; i < bins_vec.size(); i++){
        bins_vec[i].remain_capacity = max_capacity;
        for (int j = 0; j < bins_vec[i].items.size(); j++){
            bins_vec[i].remain_capacity -= bins_vec[i].items[j];
        }
    }
}

std::vector<int> bins::reduction_to_vector(){
    std::vector<int> items;
    for (int i = 0; i < bins_vec.size(); i++){
        for (int j = 0; j < bins_vec[i].items.size(); j++){
            items.push_back(bins_vec[i].items[j]);
        }
    }
    return items;
}

void bins::repair_vector(std::vector<int> items){
    bins_vec.clear();
    bin b(max_capacity);
    bins_vec.push_back(b);
    for (int i = 0; i < items.size(); i++){
        if (bins_vec.back().remain_capacity >= items[i]){
            bins_vec.back().items.push_back(items[i]);
            bins_vec.back().remain_capacity -= items[i];
        }
        else{
            bin b(max_capacity);
            bins_vec.push_back(b);
            bins_vec.back().items.push_back(items[i]);
            bins_vec.back().remain_capacity -= items[i];
        }
    }
    num_items = items.size();
    this->calc_fitness();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class bin_vec
{
    private:
        std::vector<bins> pop;
        std::vector<int> linear_scaled_fitnesses;
        std::multiset<int> items_set;
        int capacity;
    public:
        flags flag;
        bin_vec();
        bin_vec(std::vector<int> items,std::multiset<int> item_set, int size, int capacity, flags flag);
        bin_vec(std::multiset<int> item_set, int capacity, flags flag);    
        bin_vec(const bin_vec &b);
        bin_vec& operator=(const bin_vec &b);
        ~bin_vec();
        int get_capacity() const;
        void set_capacity(int capacity);
        void sort_by_fitness();
        int calc_fitness_sum();
        double calc_fitness_variance(double mean);
        std::vector<bins> get_vec() const;
        bins get_best() const;

        double calc_relative_fitness(const bins& individual);
        double set_mutation_probability_individual(bins& individual);
        double calc_SF();

        bins mate(int parent1, int parent2);
        void set(int i, bins b);
        bins get(int i) const;
        std::vector<int> repair_vector(std::vector<int> items);
        void increase_age(int i);
        
        void parent_selection_SUS(vector<int> &parents, int parents_num);
        void parent_selection_RWS(int &i1);
        void tournament_selection(int K, double P, int &selected_index);
        void liniar_scaling();

        int calc_number_of_different_alleles();

        void set_flags(flags flag){
            this->flag = flag;
        }
        flags get_flags(){
            return flag;
        }

        int distance_between_two_bins(bins b1, bins b2);
        std::vector<std::vector<int>> distance_matrix();
        std::vector<std::vector<double>> shared_distance_matrix(std::vector<std::vector<int>> &distance_matrix, double sigma_share);
        std::vector<double> f_share_fitnesses(std::vector<std::vector<double>> &shared_distance_matrix);

        int parent_selection_RWS(std::vector<double> fitnesses);

        void Niche(bin_vec &buffer, double sigma_share);
        int count_niches(double sigma_share);

        double k_means_clustering(int k);
        double k_means_clustering(int k, const std::vector<std::vector<int>>& dist_matrix);
        std::vector<std::vector<int>> k_means_clustering2(int k, const std::vector<std::vector<int>>& dist_matrix);
        int find_best_k(int max_k);
        int find_best_k2(int max_k, bin_vec &buffer);

        void non_deterministic_crowding(bin_vec &buffer, double sigma_share);

};

bin_vec::bin_vec()
{
    pop.resize(GA_POPSIZE);
}

bin_vec::bin_vec(std::vector<int> items,std::multiset<int> item_set, int size, int capacity, flags flag)
{
    pop.resize(size);
    this->items_set = item_set;
    this->capacity = capacity;
    this->flag = flag;
    for (int i = 0; i < size; i++){
        bins b(capacity);
        b.init_bins(items, flag);
        pop[i] = b;
        pop[i].set_fitness_share(pop[i].get_fitness());
    }
}

bin_vec::bin_vec(std::multiset<int> item_set, int capacity, flags flag)
{
    pop.resize(GA_POPSIZE);
    this->flag = flag;
    this->items_set = item_set;
    this->capacity = capacity;
}

bin_vec::bin_vec(const bin_vec &b)
{
    pop = b.pop;
    items_set = b.items_set;
    capacity = b.capacity;
}

bin_vec& bin_vec::operator=(const bin_vec &b)
{
    if (this == &b)
        return *this;
    pop = b.pop;
    items_set = b.items_set;
    capacity = b.capacity;
    return *this;
}

bin_vec::~bin_vec()
{
}

int bin_vec::get_capacity() const
{
    return capacity;
}

void bin_vec::set_capacity(int capacity)
{
    this->capacity = capacity;
}

void bin_vec::sort_by_fitness()
{
    std::sort(pop.begin(), pop.end(), [](const bins &a, const bins &b) -> bool{
        return a.get_fitness() < b.get_fitness();
    });
}

int bin_vec::calc_fitness_sum()
{
    int total_fitness = std::accumulate(pop.begin(), pop.end(), 0, [](int sum, const bins& b) {
    return sum + b.get_fitness();
    });
    return total_fitness;
}

double bin_vec::calc_fitness_variance(double mean)
{
    double fitness_variance = std::accumulate(pop.begin(), pop.end(), 0, [mean](double sum, const bins& b) {
    return sum + (b.get_fitness() - mean) * (b.get_fitness() - mean);
    });
    return fitness_variance;
}

std::vector<bins> bin_vec::get_vec() const
{
    return pop;
}

bins bin_vec::get_best() const
{
    return pop[0];
}

/******************************************************** */


/// the relative fitness calculation
double bin_vec::calc_relative_fitness(const bins& individual) {
    int fitness_sum = calc_fitness_sum();
    int count = pop.size();
    
    double mean = static_cast<double>(fitness_sum) / count;
    double Rf = static_cast<double>(individual.get_fitness()) / mean;
    
    /// normalize
    
    double SF = 0;
    
    SF = SF_GLOB_bin;

    double norm_Rf = Rf / SF;
    
    return norm_Rf;

}

double bin_vec::calc_SF(){
    int fitness_sum = calc_fitness_sum();
    int count = pop.size();
    double mean = static_cast<double>(fitness_sum) / count;
    double SF = 0;
    for(bins itm: pop){
        SF = SF + (itm.get_fitness()/ mean);
    } 
    return SF;
}


/// based on the relative fitness we are setting the mutation probability value 
double bin_vec::set_mutation_probability_individual(bins& individual) {
    
    double max_mutation_probability = GA_MUTATIONRATE;
    double relative_fitness = calc_relative_fitness(individual);
    
    double individual_mutation_probability = max_mutation_probability * (1.0 - relative_fitness);
    
    
    return individual_mutation_probability;

}



/**************************************************** */
bins bin_vec::mate(int parent1, int parent2)
{
    
    std::multiset<int> temp_set = items_set;

    bins p1 = pop[parent1];
    bins p2 = pop[parent2];

    vector<int> items1 = p1.reduction_to_vector();
    vector<int> items2 = p2.reduction_to_vector();

    vector<int> items;

    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_int_distribution<> dist(1, 3);
    int i = dist(g);

    if(flag.single_point_F||(flag.rand_crossover_F&&i == 1))
        items = single_point_crossover(items1, items2);
    else if(flag.two_point_F||(flag.rand_crossover_F&&i == 2))
        items = two_point_crossover(items1, items2);
    else
        items = uniform_crossover(items1, items2);        
     
    
    items = repair_vector(items);

    bins child;
    flag.shuffle_F = false;
    flag.fit_F = true;

    child.init_bins(items, flag);

    double mutation_rate;
    
    if(mutation_rate_flag_bin == 1){
        mutation_rate = GA_MUTATIONRATE;
    }
    else if(mutation_rate_flag_bin == 2){
        mutation_rate = Non_Unform_Mutation(generation_number_bin);
    }
    else if(mutation_rate_flag_bin == 3){
        
        int alleles = alils_number_bin;
        
        int threshold_smaller = 500;
        int threshold_biger = 700;
        
        mutation_rate = GA_MUTATIONRATE;

        if(alleles < threshold_smaller){
            mutation_rate = 0.8;
        }
        else if(alleles > threshold_biger){
            mutation_rate = 0.2;

        }
    }
    else if(mutation_rate_flag_bin == 4){
        mutation_rate = GA_MUTATIONRATE;

        if(global_counter_bin >= 5){

            mutation_rate = 0.9f;
            

            if((global_counter_bin > 25)){
                
                mutation_rate = GA_MUTATIONRATE;
                global_counter_bin = 0;
            }
        }

    }
    else if(mutation_rate_flag_bin == 5){
        double prob_mutation = set_mutation_probability_individual(child);

        mutation_rate = prob_mutation;
    }

    if(flag.mutation_F && (mutation_rate*RAND_MAX < rand())){
        std::random_device rd;
        std::mt19937 g(rd());
        std::uniform_int_distribution<> dist(1, 100);
        int i = dist(g);
        if(i > 0 && i < 10){
            std::vector<int> items = child.reduction_to_vector();
            int j = rand() % items.size();
            int k = rand() % items.size();
            std::swap(items[j], items[k]);
            child.repair_vector(items);
        }
    }

    return child;
}

void bin_vec::set(int i, bins b)
{
    pop[i] = b;
}

bins bin_vec::get(int i) const
{
    return pop[i];
}

std::vector<int> bin_vec::repair_vector(std::vector<int> items)
{
    std::vector<int> new_items;
    std::multiset<int> items_set_temp = this->items_set;
    for (int i = 0; i < items.size(); i++){
        if (items_set_temp.count(items[i]) > 0){
            new_items.push_back(items[i]);
            items_set_temp.erase(items_set_temp.find(items[i]));
        }
    }
    for(auto item : items_set_temp){
        new_items.push_back(item);
    }
    return new_items;
}

void bin_vec::increase_age(int i)
{
    pop[i].set_age(pop[i].get_age() + 1);
}

void bin_vec::parent_selection_SUS(vector<int> &parents, int parents_num) {
	int tsize = 9;
	int max = pop[GA_POPSIZE-1].get_fitness();
	max++; 				// Avoiding fitness = 0
	int total_fitness = 0;
	for (const auto &element : pop) { // Calculate the total fitness
		total_fitness += max-element.get_fitness();
	}
	// Select points
	int r1 = rand() % total_fitness;
	vector<int> points;

	for(int i = 0; i < parents_num; i++){
		int point = (r1 + i * total_fitness / parents_num) % total_fitness;
		points.push_back(point);
	}
    // Sort the points
	sort(points.begin(), points.end());

	// Select parents by accumulated fitness
	int current_member = 0;
	int comulative_sum = 0;
	for (int point: points) {
		while(comulative_sum + max - pop[current_member].get_fitness() < point){
			comulative_sum += max - pop[current_member].get_fitness();
			current_member++;
		}
		parents.push_back(current_member);
	}
}

void bin_vec::parent_selection_RWS(int &i1) {
	    
    int tsize = 9;

	int max =  linear_scaled_fitnesses[GA_POPSIZE-1];
	max++; // Avoiding fitness = 0
	int total_fitness = 0;
	for (const auto &element : linear_scaled_fitnesses) {
		total_fitness += max-element;
	}
	int r1 = rand() % total_fitness;
	int sum = 0;
	for (int i = 0; i < GA_POPSIZE; i++) {
		sum += max-linear_scaled_fitnesses[i];
		if (sum >= r1) {
			i1 = i;
			break;
		}
	}
}

void bin_vec::tournament_selection(int K, double P, int &selected_index) {
    int tournament_winner_index = 0;
    unsigned int best_fitness = std::numeric_limits<unsigned int>::max(); ///Initialize to maximum possible value
    ///Perform K tournaments
    for (int i = 0; i < K; ++i) {        
        ///Randomly select an individual from the population
        int index = rand() % pop.size();
        bins competitor = pop[index];

        ///Determine if the competitor wins the current round of the tournament based on probability P
        if (    (((double)rand() / RAND_MAX) < P) 
            &&  (competitor.get_fitness() < best_fitness)) {
            
            best_fitness = competitor.get_fitness();
            tournament_winner_index = index;
        }
    }
    ///Return the index of the tournament winner
    selected_index = tournament_winner_index;
}

void bin_vec::liniar_scaling(){
    unsigned int max_value = 0;
    unsigned int min_value = numeric_limits<unsigned int>::max(); // Initialize to maximum possible value
    unsigned int sum = 0;

    linear_scaled_fitnesses.resize(pop.size());

    // Find max, min, and sum of fitness values
    int i=1;
    
    sum = calc_fitness_sum();
    min_value = pop[0].get_fitness();
    max_value = pop[pop.size()-1].get_fitness();

    // Calculate average fitness
    double avrg = sum / pop.size();

    // Calculate scaling parameters
    double a = static_cast<double>(1 - i) / (max_value - min_value);
    double b = 1 - a * avrg;

    // Scale fitness values and create a new scaled population
    for(int i = 0; i < pop.size(); i++){
        linear_scaled_fitnesses[i] = a * pop[i].get_fitness() + b;
    }
}

int bin_vec::calc_number_of_different_alleles(){
    
    std::set<std::set<std::set<int>>> alleles;
    for (const auto& individual : pop) {
        std::set<std::set<int>> bins;
        for (const auto& bin : individual.get_bins()) {
            std::set<int> bin_items(bin.items.begin(), bin.items.end());

            bins.insert(bin_items);
        }
        alleles.insert(bins);
    }

    return alleles.size();
}

int bin_vec::distance_between_two_bins(bins b1, bins b2){
    int distance = 0;
    std::vector<int> items1 = b1.reduction_to_vector();
    std::vector<int> items2 = b2.reduction_to_vector();
    
    for(int i = 0; i < items1.size(); i++){
        if(items1[i] != items2[i]){
            int temp = items2[i];
            // Find the index of the item in the second vector
            int index = -1;
            for(int j = i; j < items2.size(); j++){
                if(items2[j] == items1[i]){
                    index = j;
                    break;
                }
            }
            // Swap the items
            items2[index] = temp;
            
            distance++;
        }
    }

    return distance;
}

std::vector<std::vector<int>> bin_vec::distance_matrix(){
    std::vector<std::vector<int>> distance_matrix(pop.size(), std::vector<int>(pop.size(), 0));
    for(int i = 0; i < pop.size(); i++){
        for(int j = i+1; j < pop.size(); j++){
            distance_matrix[i][j] = distance_between_two_bins(pop[i], pop[j]);
            distance_matrix[j][i] = distance_matrix[i][j];
        }
    }
    return distance_matrix;
}

std::vector<std::vector<double>> bin_vec::shared_distance_matrix(std::vector<std::vector<int>> &distance_matrix, double sigma_share){
    std::vector<std::vector<double>> shared_distance_matrix(pop.size(), std::vector<double>(pop.size(), 0));
    for(int i = 0; i < pop.size(); i++){
        for(int j = i; j < pop.size(); j++){
            if (distance_matrix[i][j] <= sigma_share){                
                shared_distance_matrix[i][j] = 1 - distance_matrix[i][j] / sigma_share;
                shared_distance_matrix[j][i] = shared_distance_matrix[i][j];
            }
        }
    }
    return shared_distance_matrix;
}

std::vector<double> bin_vec::f_share_fitnesses(std::vector<std::vector<double>> &shared_distance_matrix){
    std::vector<double> f_share_fitnesses;

    for(int i = 0; i < pop.size(); i++){
        double f_share = 0;
        for(int j = 0; j < pop.size(); j++){
            f_share += shared_distance_matrix[i][j];
        }
        double raw_fit = pop[i].get_fitness();
        f_share_fitnesses.push_back(raw_fit / f_share);
    }
    return f_share_fitnesses;
}

int bin_vec::parent_selection_RWS(std::vector<double> fitnesses){
    double max =  *std::max_element(fitnesses.begin(), fitnesses.end());
    max++; // Avoiding fitness = 0
    double total_fitness = 0;
    for (const auto &element : fitnesses) {
        total_fitness += max-element;
    }
    // random number between 0 and total_fitness
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, total_fitness);
    double r1 = dis(gen);

    double sum = 0;
    int index = 0;
    for (int i = 0; i < GA_POPSIZE; i++) {
        sum += max-fitnesses[i];
        if (sum >= r1) {
            return i;
            index = i;
        }
    }
    return index;
}

/*void bin_vec::Niche(bin_vec &buffer, double sigma_share){
    std::vector<std::vector<int>> distance_matrix = this->distance_matrix();
    std::vector<std::vector<double>> shared_distance_matrix = this->shared_distance_matrix(distance_matrix, sigma_share);
    std::vector<double> f_share_fitnesses = this->f_share_fitnesses(shared_distance_matrix);
    std::vector<int> indexes;
    
    int counter = 0;
    for (int i = 0; i < pop.size(); i++){
        indexes.push_back(i);
    }
    
    for (int i = 0; i < f_share_fitnesses.size(); i++){
        std::cout << f_share_fitnesses[i] << " ";
    }
  

    for(int i = 0; i < buffer.get_vec().size(); i++){
        int index = parent_selection_RWS(f_share_fitnesses);
        //std::cout << "index: " << index << std::endl;

        std::shuffle(indexes.begin(), indexes.end(), std::mt19937(std::random_device()()));
        for (int j = 0; j < indexes.size(); j++){
            if (indexes[j] != index){
                if (distance_matrix[index][indexes[j]] < sigma_share){
                    //std::cout << "indexes[j]: " << indexes[j] << std::endl;
                    bins child = mate(index, indexes[j]);
                    buffer.set(i, child);
                    shared_distance_matrix[index][indexes[j]] += 1 - (distance_matrix[index][indexes[j]] / sigma_share);
                    
                    shared_distance_matrix[indexes[j]][index] = shared_distance_matrix[index][indexes[j]];
                    
                    f_share_fitnesses[index] = pop[index].get_fitness() / shared_distance_matrix[index][indexes[j]];
                    break;
                }
                else if (j == pop.size()-1){
                    i--;
                    counter++;
                    if (counter > 1000){
                        std::cout << "Error: Niche algorithm is stuck in a loop" << std::endl;
                        sigma_share+=10;
                        counter = 0;
                        if (sigma_share < 0){
                            std::cout << "Error: sigma_share is negative" << std::endl;
                            exit(1);
                        }
                    }
                
                }
                else{
                    //std::cout<<"Error: 111111111111111111111111111"<<endl;
                }
            }

        }

    }
}*/

void bin_vec::Niche(bin_vec &buffer, double sigma_share) {

    int e_size = pop.size()*GA_ELITRATE;
    elitism(*this, buffer, e_size);
    // Precompute the distance and shared distance matrices
    std::vector<std::vector<int>> distance_matrix = this->distance_matrix();
    std::vector<std::vector<double>> shared_distance_matrix = this->shared_distance_matrix(distance_matrix, sigma_share);
    std::vector<double> f_share_fitnesses = this->f_share_fitnesses(shared_distance_matrix);

    std::vector<int> indexes(pop.size());
    std::iota(indexes.begin(), indexes.end(), 0); // Fill indexes with values 0 to pop.size()-1

    int counter = 0;
    for (int i = e_size; i < buffer.get_vec().size(); i++) {
        int index = parent_selection_RWS(f_share_fitnesses);

        std::shuffle(indexes.begin(), indexes.end(), std::mt19937(std::random_device()()));

        bool mate_found = false;
        for (int j = 0; j < indexes.size(); j++) {
            if (indexes[j] != index) {
                if (distance_matrix[index][indexes[j]] < sigma_share) {
                    // Mate selected, produce offspring
                    bins child = mate(index, indexes[j]);
                    buffer.set(i, child);

                    // Update shared distance matrix and fitness values
                    double shared_value = 1 - (distance_matrix[index][indexes[j]] / sigma_share);
                    shared_distance_matrix[index][indexes[j]] += shared_value;
                    shared_distance_matrix[indexes[j]][index] = shared_distance_matrix[index][indexes[j]];

                    // Recompute shared fitness for the selected index
                    f_share_fitnesses[index] = pop[index].get_fitness() / shared_distance_matrix[index][indexes[j]];

                    mate_found = true;
                    break;
                }
            }
        }

        if (!mate_found) {
            // If no mate found within sigma_share, re-try or adjust sigma_share
            i--;
            counter++;
            if (counter > 1000) {
                std::cout << "Error: Niche algorithm is stuck in a loop" << std::endl;
                sigma_share += 10; // Adjust sigma_share cautiously
                counter = 0;
            }
        }
    }
    std::cout << "Sigma share: " << sigma_share << std::endl;
    std::cout << "Niches: " << count_niches(sigma_share) << std::endl;
}

int bin_vec::count_niches(double sigma_share) {
    std::vector<std::vector<int>> distance_matrix = this->distance_matrix();
    std::vector<bool> visited(pop.size(), false);
    int niche_count = 0;

    // Iterate through each individual in the population
    for (int i = 0; i < pop.size(); ++i) {
        if (!visited[i]) {
            niche_count++;
            // Start a new niche
            std::set<int> current_niche;
            current_niche.insert(i);
            visited[i] = true;

            // Find all individuals within sigma_share distance and mark them as visited
            for (int j = 0; j < pop.size(); ++j) {
                if (!visited[j] && distance_matrix[i][j] < sigma_share) {
                    current_niche.insert(j);
                    visited[j] = true;
                }
            }
            std::cout << "Niche " << niche_count << ": ";
            for (int member : current_niche) {
                std::cout << member << " ";
            }
            std::cout << std::endl;
        }
    }

    return niche_count;
}


// The k-means clustering method
double bin_vec::k_means_clustering(int k) {
    // Step 1: Initialize centroids randomly from the population
    std::vector<bins> centroids;
    std::sample(pop.begin(), pop.end(), std::back_inserter(centroids), k, std::mt19937{std::random_device{}()});

    std::vector<int> assignments(pop.size(), -1);  // Cluster assignments for each bin
    bool changed;
    double total_distance = 0.0;

    do {
        changed = false;
        total_distance = 0.0;

        // Step 2: Assign each bin to the nearest centroid
        for (size_t i = 0; i < pop.size(); ++i) {
            int nearest_centroid = 0;
            int min_distance = distance_between_two_bins(pop[i], centroids[0]);

            for (int j = 1; j < k; ++j) {
                int distance = distance_between_two_bins(pop[i], centroids[j]);
                if (distance < min_distance) {
                    min_distance = distance;
                    nearest_centroid = j;
                }
            }

            // Accumulate the distance for total distance calculation
            total_distance += min_distance;

            // If the assignment changes, mark that a change occurred
            if (assignments[i] != nearest_centroid) {
                assignments[i] = nearest_centroid;
                changed = true;
            }
        }

        // Step 3: Recalculate centroids
        if (changed) {
            std::vector<std::vector<int>> clusters(k);
            for (size_t i = 0; i < pop.size(); ++i) {
                clusters[assignments[i]].push_back(i);
            }

            for (int j = 0; j < k; ++j) {
                if (!clusters[j].empty()) {
                    // Randomly select a new centroid from the cluster members
                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::uniform_int_distribution<> dist(0, clusters[j].size() - 1);
                    centroids[j] = pop[clusters[j][dist(gen)]];
                }
            }
        }
    } while (changed);

    // Output the minimum total distance after clustering
    std::cout << "Minimum total distance: " << total_distance << std::endl;

    return total_distance;

    // At this point, assignments represent the cluster each bin belongs to
    // You can further process these clusters as needed
}


double bin_vec::k_means_clustering(int k, const std::vector<std::vector<int>>& dist_matrix) {
    std::vector<bins> centroids;
    std::sample(pop.begin(), pop.end(), std::back_inserter(centroids), k, std::mt19937{std::random_device{}()});

    std::vector<int> assignments(pop.size(), -1);
    bool changed;
    double total_distance = 0.0;

    do {
        changed = false;
        total_distance = 0.0;

        // Step 2: Assign each bin to the nearest centroid
        for (size_t i = 0; i < pop.size(); ++i) {
            int nearest_centroid = 0;
            int min_distance = dist_matrix[i][0];

            for (int j = 1; j < k; ++j) {
                int distance = dist_matrix[i][j];
                if (distance < min_distance) {
                    min_distance = distance;
                    nearest_centroid = j;
                }
            }

            total_distance += min_distance;

            if (assignments[i] != nearest_centroid) {
                assignments[i] = nearest_centroid;
                changed = true;
            }
        }

        // Step 3: Recalculate centroids
        if (changed) {
            std::vector<std::vector<int>> clusters(k);
            for (size_t i = 0; i < pop.size(); ++i) {
                clusters[assignments[i]].push_back(i);
            }

            for (int j = 0; j < k; ++j) {
                if (!clusters[j].empty()) {
                    // Choose a random bin as the new centroid
                    std::uniform_int_distribution<int> distribution(0, clusters[j].size() - 1);
                    std::random_device rd;
                    std::mt19937 gen(rd());
                    int random_index = distribution(gen);
                    centroids[j] = pop[clusters[j][random_index]];
                }
            }
        }
    } while (changed);

    return total_distance;
}

std::vector<std::vector<int>> bin_vec::k_means_clustering2(int k, const std::vector<std::vector<int>>& dist_matrix) {
    std::vector<bins> centroids;
    std::sample(pop.begin(), pop.end(), std::back_inserter(centroids), k, std::mt19937{std::random_device{}()});

    std::vector<int> assignments(pop.size(), -1);
    bool changed;
    double total_distance = 0.0;
    std::vector<std::vector<int>> final_clusters(k);


    do {
        changed = false;
        total_distance = 0.0;

        // Step 2: Assign each bin to the nearest centroid
        for (size_t i = 0; i < pop.size(); ++i) {
            int nearest_centroid = 0;
            int min_distance = dist_matrix[i][0];

            for (int j = 1; j < k; ++j) {
                int distance = dist_matrix[i][j];
                if (distance < min_distance) {
                    min_distance = distance;
                    nearest_centroid = j;
                }
            }

            total_distance += min_distance;

            if (assignments[i] != nearest_centroid) {
                assignments[i] = nearest_centroid;
                changed = true;
            }
        }

        // Step 3: Recalculate centroids
        if (changed) {
            std::vector<std::vector<int>> clusters(k);
            for (size_t i = 0; i < final_clusters.size(); ++i) {
                final_clusters[i].clear();
            }
            for (size_t i = 0; i < pop.size(); ++i) {
                clusters[assignments[i]].push_back(i);
                final_clusters[assignments[i]].push_back(i);
            }

            for (int j = 0; j < k; ++j) {
                if (!clusters[j].empty()) {
                    // Choose a random bin as the new centroid
                    std::uniform_int_distribution<int> distribution(0, clusters[j].size() - 1);
                    std::random_device rd;
                    std::mt19937 gen(rd());
                    int random_index = distribution(gen);
                    centroids[j] = pop[clusters[j][random_index]];
                }
            }
        }
    } while (changed);

    return final_clusters;
}

int bin_vec::find_best_k(int max_k) {
    std::vector<double> distances;

    for (int k = 2; k <= max_k; ++k) {
        double total_distance = k_means_clustering(k);
        distances.push_back(total_distance);
        std::cout << "k = " << k << ", Total Distance = " << total_distance << std::endl;
    }

    // Find the "elbow" point
    double best_k = 2;
    double min_difference = std::numeric_limits<double>::max();

    for (int i = 2; i < distances.size() - 1; ++i) {
        double diff1 = distances[i] - distances[i - 1];
        double diff2 = distances[i + 1] - distances[i];
        if (diff1 > diff2 && diff1 < min_difference) {
            min_difference = diff1;
            best_k = i + 1;
        }
    }

    std::cout << "Best k based on elbow method: " << best_k << std::endl;
    return best_k;
}


int bin_vec::find_best_k2(int max_k, bin_vec &buffer) {
    std::map<int, double> k_to_distance;
    auto dist_matrix = distance_matrix();  // Precompute distance matrix

    for (int k = 2; k <= max_k; ++k) {
        double distance = k_means_clustering(k, dist_matrix);
        k_to_distance[k] = distance;
        //std::cout << "k = " << k << ", Total Distance = " << distance << std::endl;
    }

    // Determine the "elbow" point
    int best_k = 2;
    double max_diff = 0.0;

    for (int k = 3; k <= max_k; ++k) {
        double diff1 = k_to_distance[k - 1] - k_to_distance[k];
        double diff2 = k_to_distance[k] - k_to_distance[k - 2];
        if (diff1 > diff2 && diff1 > max_diff) {
            max_diff = diff1;
            best_k = k;
        }
    }

    std::cout << "Best k based on elbow method: " << best_k << std::endl;

    std::vector<std::vector<int>> clusters = k_means_clustering2(best_k, dist_matrix);
    
    int esize = pop.size()*GA_ELITRATE;
    elitism(*this, buffer, esize);

    for (int i = esize; i < pop.size(); i++){
        int random_individual = rand() % pop.size();
        int cluster_index = 0;
        for (int j = 0; j < clusters.size(); j++){
            if (std::find(clusters[j].begin(), clusters[j].end(), random_individual) != clusters[j].end()){
                cluster_index = j;
                break;
            }
        }
        bins child = mate(clusters[cluster_index][rand() % clusters[cluster_index].size()], random_individual);
        buffer.set(i, child);
    }

    return best_k;
}

void bin_vec::non_deterministic_crowding(bin_vec &buffer, double T){

    int e_size = pop.size()*GA_ELITRATE;
    elitism(*this, buffer, e_size);
    
    for(int i = e_size; i < buffer.get_vec().size(); i++){
        int index = 0;
        liniar_scaling();
        parent_selection_RWS(index);

        bins child1 = mate(index, i);

        bins child2 = mate(i, index);
        
        double parent1_fitness = pop[index].get_fitness();
        double parent2_fitness = pop[i].get_fitness();
        double child1_fitness = child1.get_fitness();
        double child2_fitness = child2.get_fitness();
        if (child1_fitness > parent1_fitness) {
            double P_replace = exp((child1_fitness - parent1_fitness) / T);
            if ((rand() / double(RAND_MAX)) < P_replace) {
                buffer.set(i, child1);
            }
            else {
                buffer.set(i, pop[i]);
            }
        } else if(child2_fitness > parent2_fitness) {
            double P_replace = exp((child2_fitness - parent2_fitness) / T);
            if ((rand() / double(RAND_MAX)) < P_replace) {
                buffer.set(i, child2);
            }
            else {
                buffer.set(i, pop[index]);
            }
            buffer.set(i, pop[index]);
        }
        else{
            buffer.set(i, pop[i]);
        }       
        

    }
}


///********************************************************************************************************************************************************************

void read_the_file(std::ifstream &infile, std::string &line, int &num_problems,
 std::string &problem_identifier, int &bin_capacity, int &num_items, int &num_bins,
 std::vector<std::vector<int>> &items_of_items, std::vector<int> &capacities, std::multiset<int> &items_set)
{
    // Read the first line to get the number of problems
    if (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (!(iss >> num_problems)) {
            std::cerr << "Error reading number of problems" << std::endl;
            exit(1);
        }
    }

    for(int i = 0; i < num_problems; i++){
        // Read the second line to get the problem identifier
        if (std::getline(infile, line)) {
            problem_identifier = line;
        }

        // Read the third line to get the bin capacity, number of items, and number of bins
        if (std::getline(infile, line)) {
            std::istringstream iss(line);
            if (!(iss >> bin_capacity >> num_items >> num_bins)) {
                std::cerr << "Error reading bin capacity, number of items, and number of bins" << std::endl;
                exit(1);
            }
        }
        capacities.push_back(bin_capacity);
        std::vector<int> items;
        // Read the remaining lines to get the items
        for(int j = 0; j < num_items; j++){
            std::getline(infile, line);
            std::istringstream iss(line);
            int item;
            if (iss >> item) {
                items.push_back(item);
            }
        }
        items_of_items.push_back(items);
    }
}

/////**********************************************************************************/////
/***/
///this function recieves an array of data and saves this data to a file.
void write_numbers_to_file(const std::string &filename, const std::vector<int> &numbers) {
    std::ofstream outfile(filename);

    if (!outfile) {
        std::cerr << "Error opening file for writing" << std::endl;
        return;
    }

    for (int number : numbers) {
        outfile << number << std::endl;
    }

    outfile.close();
    std::cout << "Numbers written to " << filename << std::endl;
}


/////**********************************************************************************/////



int main (int argc, char* argv[])
{
    flags flag;
    int i =1;
    int index = 5;
    int niche_size = 100, num_clusters = 30;

    while (i < argc){
        if (strcmp(argv[i], "-RWS") == 0){
            flag.RWS_F = true;
        }
        else if (strcmp(argv[i], "-SUS") == 0){
            flag.SUS_F = true;
        }
        else if (strcmp(argv[i], "-PMX") == 0){
            flag.PMX_F = true;
        }
        else if (strcmp(argv[i], "-CX") == 0){
            flag.CX_F = true;
        }
        else if (strcmp(argv[i], "-inversion") == 0){
            flag.inversion_F = true;
        }
        else if (strcmp(argv[i], "-scramble") == 0){
            flag.scramble_F = true;
        }
        else if (strcmp(argv[i], "-simple") == 0){
            flag.simple_F = true;
        }
        else if (strcmp(argv[i], "-elitism") == 0){
            flag.elitism_F = true;
        }
        else if (strcmp(argv[i], "-age") == 0){
            flag.age_F = true;
        }
        else if (strcmp(argv[i], "-linear_scaling") == 0){
            flag.linear_scaling_F = true;
        }
        else if (strcmp(argv[i], "-tournament") == 0){
            flag.tournament_F = true;
        }
        else if (strcmp(argv[i], "-mutation") == 0){
            flag.mutation_F = true;
        }
        else if (strcmp(argv[i], "-index") == 0){
            index = atoi(argv[i+1]);
            i++;
        }
        else if (strcmp(argv[i], "-best_fit") == 0){
            flag.best_fit_F = true;
        }
        else if (strcmp(argv[i], "-worst_fit") == 0){
            flag.worst_fit_F = true;
        }
        else if (strcmp(argv[i], "-random") == 0){
            flag.random_F = true;
        }
        else if (strcmp(argv[i], "-greedy") == 0){
            flag.greedy_F = true;
        }
        else if (strcmp(argv[i], "-nshuffle") == 0){
            flag.shuffle_F = false;
        }
        else if (strcmp(argv[i], "-nfit") == 0){
            flag.fit_F = false;
        }
        else if (strcmp(argv[i], "-single_point") == 0){
            flag.single_point_F = true;
        }
        else if (strcmp(argv[i], "-two_point") == 0){
            flag.two_point_F = true;
        }
        else if (strcmp(argv[i], "-uniform") == 0){
            flag.uniform_F = true;
        }
        else if (strcmp(argv[i], "-nrand_crossover") == 0){
            flag.rand_crossover_F = false;
        }
        else if (strcmp(argv[i], "-niche") == 0){
            flag.niching_F = true;
        }
        else if (strcmp(argv[i], "-crowding") == 0){
            flag.crowding_F = true;
            flag.niching_F = false;
        }
        else if (strcmp(argv[i], "-clustering") == 0){
            flag.clustering_F = true;
            flag.niching_F = false;
        }
        else if (strcmp(argv[i], "-nniche") == 0){
            flag.niching_F = false;
            flag.simple_mate_F = true;
        }
        else if (strcmp(argv[i], "-niche_size") == 0){
            niche_size = atoi(argv[i+1]);
            i++;
        }
        else if (strcmp(argv[i], "-num_clusters") == 0){
            num_clusters = atoi(argv[i+1]);
            i++;
        }

    

        else{
            std::cout << "Invalid flag" << std::endl;
            return 1;
        }
        i++;
    }


    // Random seed
    // Seed with a real random value, if available
    random_device r;
    // Choose a random number in the range [0, 10]
    uniform_int_distribution<int> uniform_dist(0, 10);
    default_random_engine e1(r());
    for(int i = 0; i < uniform_dist(e1); i++)
    {
        rand();
    }
    
    
    std::ifstream infile("binpack.txt"); // Open the file
    if (!infile.is_open()) { // Check if the file is opened successfully
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }

    std::string line;
    int num_problems;
    std::string problem_identifier;
    int bin_capacity;
    int num_items;
    int num_bins;
    std::vector<std::vector<int>> items_of_items;
    std::vector<int> capacities;
    std::multiset<int> items_set;


    read_the_file(infile, line, num_problems, problem_identifier, bin_capacity, num_items, num_bins, items_of_items, capacities, items_set);

   
    infile.close(); // Close the file

    // Output the results
    std::cout << "Number of problems: " << num_problems << std::endl;
    std::cout << "Problem identifier: " << problem_identifier << std::endl;
    std::cout << "Bin capacity: " << bin_capacity << std::endl;
    std::cout << "Number of items: " << num_items << std::endl;
    std::cout << "Number of bins: " << num_bins << std::endl;
    std::cout << "Items: ";

    std::cout << "index: " << index << endl;

    std::vector<double> average;
    std::vector<double> variance;

    std::vector<int> items = items_of_items[index];

    for (int item : items) {
        std::cout << item << " ";
    }
    std::cout << std::endl;

    items_set.insert(items.begin(), items.end());

    
    
    bin_vec pop(items, items_set, GA_POPSIZE, capacities[index], flag);
    bin_vec buffer(items_set, capacities[index], flag);
    bin_vec *population = &pop;
    bin_vec *buffer_ptr = &buffer;
    population->flag = flag;
    buffer_ptr->flag = flag;

    ///****************************************************************** */
    int mutation_rate_flag_;
    std::cout << "Enter mutation rate flag number: " << endl;
    std::cout << "1. basic" << endl;
    std::cout << "2. Non-Unform" << endl;
    std::cout << "3. Adaptive" << endl;
    std::cout << "4. Triggered Hyper" << endl;
    std::cout << "5. Self-Adaptive" << endl << endl;    
    std::cin >> mutation_rate_flag_;

    std::cout << "flag is " << mutation_rate_flag_ << endl;
    
    mutation_rate_flag_bin = mutation_rate_flag_;

    ///****************************************************************** */
    
    std::vector<int> numbers;

    std::vector<int> fitness_vec;
    for (int i=0; i<GA_MAXITER; i++){

        population->sort_by_fitness();

        /****** */
        generation_number_bin = i + 1;
        alils_number_bin = population->calc_number_of_different_alleles();
        SF_GLOB_bin = population->calc_SF();
        /****** */

        std::cout << "Generation: " << i+1 << " Fitness: " << population->get_best().get_fitness() << std::endl;
        
        int lell = population->calc_number_of_different_alleles();
        std::cout << "Number of different alleles: " << lell << std::endl;
        numbers.push_back(lell);

        population->get_best().print_bins();

        top_average_and_variance(*population, average, variance);
        
        /******** */
        curr_best_fitness_bin = population->get_best().get_fitness();  
        if(prev_best_fitness_bin == curr_best_fitness_bin){
            global_counter_bin = global_counter_bin + 1;
        }else{
            global_counter_bin = 0;
        }
        /******* */

        fitness_vec.push_back(population->get_best().get_fitness());
        if (population->get_best().get_fitness()==0){
            break;
        }

        if (flag.simple_mate_F) mate(*population, *buffer_ptr);

        if (flag.niching_F) population->Niche(*buffer_ptr, niche_size);

        if (flag.clustering_F) int best_k = population->find_best_k2(num_clusters, *buffer_ptr);
        
        if (flag.crowding_F) population->non_deterministic_crowding(*buffer_ptr, 0.1);

        /****** */
        prev_best_fitness_bin = curr_best_fitness_bin;   
        /****** */       
        swap(population, buffer_ptr);      
    }

    write_numbers_to_file("allels_BINS.txt", numbers);
    write_numbers_to_file("fitness_BINS.txt", fitness_vec);
    write_numbers_to_file("average_BINS.txt", average);
    write_numbers_to_file("variance_BINS.txt", variance);


    return 0;

}