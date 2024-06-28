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
#include <set>
#include <random>
#include <fstream>
#include <numeric>
#include <functional>

#include "t_functions.h"



#define BALD_STRING  "01110100110101001010"
#define MAX_LEARNING 1000

using namespace std;

class baldwin
{
private:
    double fitness;
    std::string str; 

public:
    baldwin();
    double get_fitness() { return fitness; }
    void calc_fitness(int learning);
    void init_str();
    std::string get_str() { return str; }
    void set_str(std::string s) { str = s; }
    void learn();
    double correct_position_percentage();
    double incorrect_position_percentage();
};

baldwin::baldwin()
{
    str.resize(strlen(BALD_STRING));
}

void baldwin::calc_fitness(int learning)
{
    fitness = 1 + double(19 * learning) / MAX_LEARNING;
}

void baldwin::init_str(){
     // Random number generation setup
    random_device rd; // Obtain a random number from hardware
    mt19937 gen(rd()); // Seed the generator
    uniform_real_distribution<> dis(0.0, 1.0); // Define the range

    // Initialize the str with the same length as BOLD_SRING
    str.resize(strlen(BALD_STRING));

    std::set<int> indexes;
    for (int i = 0; i < str.size(); i++){
        indexes.insert(i);
    }
    int unknowns = str.size() / 2;
    int right = str.size() / 4;
    int wrong = str.size() - unknowns - right;

    for (int i = 0; i < right; i++){
        int index = rand() % str.size();
        if (indexes.find(index) != indexes.end()){
            str[index] = BALD_STRING[index];
            indexes.erase(index);
        }
        else{
            i--;
        }
    }
    for (int i = 0; i < wrong; i++){
        int index = rand() % str.size();
        if (indexes.find(index) != indexes.end()){
            if (BALD_STRING[index] == '1'){
                str[index] = '0';
            }
            else{
                str[index] = '1';
            }
            indexes.erase(index);
        }
        else{
            i--;
        }
    }
    for (int i = 0; i < unknowns; i++){
        int index = rand() % str.size();
        if (indexes.find(index) != indexes.end()){
            str[index] = '?';
            indexes.erase(index);
        }
        else{
            i--;
        }
    }

    learn();
}

void baldwin::learn()
{
    // Random number generation setup
    random_device rd; // Obtain a random number from hardware
    // Choose a random number in the range [0, 1]
    uniform_int_distribution<int> uniform_dist(0, 1);
    default_random_engine e1(rd());

    int learning = MAX_LEARNING;

    for (int i = 0; i < MAX_LEARNING; i++)
    {
        std::string new_str = str;
        for (int j = 0; j < new_str.size(); j++)
        {
            if (new_str[j] == '?')
            {
                uniform_dist(e1) == 1 ? new_str[j] = '1' :new_str[j] = '0';
            }
        }
        if (new_str == BALD_STRING)
        {
            learning = i;
            break;
        }
    }

    calc_fitness(MAX_LEARNING - learning);
}

double baldwin::correct_position_percentage()
{
    double correct = 0;
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == BALD_STRING[i])
        {
            correct++;
        }
    }
    return correct / str.size();
}

double baldwin::incorrect_position_percentage()
{
    double incorrect = 0;
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] != BALD_STRING[i] && str[i] != '?')
        {
            incorrect++;
        }
    }
    return incorrect / str.size();
}

/*void baldwin::correct_incorrect_position_percentage(&correct, &incorrect)
{
    double correct = 0;
    double incorrect = 0;
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == BALD_STRING[i])
        {
            correct++;
        }
        else
        {
            incorrect++;
        }
    }
    correct /= str.size();
    incorrect /= str.size();
}*/

class baldwin_vec
{
private:
    vector<baldwin> pop;
public:
    flags flag;
    baldwin_vec();
    void init_pop();
    void sort_by_fitness();
    baldwin &get_best();
    flags get_flags() { return flag; }
    void parent_selection_SUS(vector<int> &parents, int n);
    void parent_selection_RWS(int &parent);
    void tournament_selection(int k, double p, int &parent);
    void set(int i, baldwin b) { pop[i] = b; }
    baldwin mate(int p1, int p2);
    baldwin& get(int i) { return pop[i]; }
    void increase_age(int i);
    baldwin single_point_crossover(int p1, int p2);
    baldwin two_point_crossover(int p1, int p2);
    baldwin uniform_crossover(int p1, int p2);
    double mean_correct_position_percentage();
    double mean_incorrect_position_percentage();
};

baldwin_vec::baldwin_vec()
{
    pop.resize(GA_POPSIZE);

    flag.age_F = false;
    flag.best_fit_F = false;
    flag.CX_F = false;
    flag.elitism_F = false;
    flag.fit_F = false;
    flag.greedy_F = false;
    flag.inversion_F = false;
    flag.linear_scaling_F = false;
    flag.mutation_F = false;
    flag.PMX_F = false;
    flag.random_F = false;
    flag.RWS_F = true;
    flag.scramble_F = false;
    flag.shuffle_F = false;
    flag.simple_F = false;
    flag.SUS_F = false;
    flag.tournament_F = false;
    flag.two_point_F = false;
    flag.uniform_F = false;
    flag.worst_fit_F = false;
    flag.single_point_F = false;
    flag.rand_crossover_F = false;
}

void baldwin_vec::init_pop()
{
    for (int i = 0; i < pop.size(); ++i)
    {
        pop[i].init_str();
        std::cout << "Individual " << i << " : " << pop[i].get_str() << std::endl;
    }
    
}

void baldwin_vec::sort_by_fitness()
{
    sort(pop.begin(), pop.end(), [](baldwin &a, baldwin &b) { return a.get_fitness() > b.get_fitness(); });
}

baldwin &baldwin_vec::get_best()
{
    return pop[0];
}

void baldwin_vec::parent_selection_SUS(vector<int> &parents, int n)
{
    vector<double> fitness(pop.size());
    vector<double> cum_fitness(pop.size());
    vector<double> r(n);

    // Calculate the fitness of each individual
    for (size_t i = 0; i < pop.size(); ++i)
    {
        fitness[i] = pop[i].get_fitness();
    }

    // Calculate the cumulative fitness
    partial_sum(fitness.begin(), fitness.end(), cum_fitness.begin());

    // Calculate the average fitness
    double avg_fitness = cum_fitness.back() / pop.size();

    // Calculate the selection pressure
    double selection_pressure = 2.0;

    // Calculate the distance between the pointers
    double pointer_distance = avg_fitness / selection_pressure;

    // Generate n random numbers between 0 and 1
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    for (size_t i = 0; i < n; ++i)
    {
        r[i] = dis(gen);
    }

    // Select the parents
    parents.clear();
    for (size_t i = 0; i < n; ++i)
    {
        double pointer = r[i] * pointer_distance;
        for (size_t j = 0; j < pop.size(); ++j)
        {
            if (cum_fitness[j] >= pointer)
            {
                parents.push_back(j);
                break;
            }
        }
    }
}

void baldwin_vec::parent_selection_RWS(int &parent)
{

    vector<double> fitness(pop.size());
    vector<double> cum_fitness(pop.size());

    // Calculate the fitness of each individual
    for (size_t i = 0; i < pop.size(); ++i)
    {
        fitness[i] = pop[i].get_fitness();
    }

    // Calculate the cumulative fitness
    partial_sum(fitness.begin(), fitness.end(), cum_fitness.begin());

    // Generate a random number between 0 and the total fitness
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, cum_fitness.back());
    double r = dis(gen);

    // Select the parent
    for (size_t i = 0; i < pop.size(); ++i)
    {
        if (cum_fitness[i] >= r)
        {
            parent = i;
            break;
        }
    }
}

void baldwin_vec::tournament_selection(int k, double p, int &parent)
{
    return;
}

void baldwin_vec::increase_age(int i)
{
    return;
}

baldwin baldwin_vec::mate(int p1, int p2)
{
    random_device rd;
    default_random_engine e1(rd());
    uniform_int_distribution<> dist(1,3);
    int cross = dist(e1);

    if (flag.single_point_F) cross = 1;
    else if (flag.two_point_F) cross = 2;
    else if (flag.uniform_F) cross = 3;

    baldwin child;

    if (cross == 1) child = single_point_crossover(p1, p2);
    else if (cross == 2) child = two_point_crossover(p1, p2);
    else if (cross == 3) child = uniform_crossover(p1, p2);  


    child.learn();

    //std::cout << child.get_str() << std::endl;


    return child;

}

baldwin baldwin_vec::single_point_crossover(int p1, int p2)
{
    int cross_point = rand() % pop[p1].get_str().size();
    baldwin child;

    std::string str1 = pop[p1].get_str().substr(0, cross_point);
    std::string str2 = pop[p2].get_str().substr(cross_point);

    child.set_str(str1 + str2);

    return child;
}

baldwin baldwin_vec::two_point_crossover(int p1, int p2)
{
    int cross_point1 = rand() % pop[p1].get_str().size();
    int cross_point2 = rand() % pop[p1].get_str().size();

    if (cross_point1 > cross_point2) swap(cross_point1, cross_point2);
    baldwin child;

    std::string str1 = pop[p1].get_str().substr(0, cross_point1);
    std::string str2 = pop[p2].get_str().substr(cross_point1, cross_point2 - cross_point1);
    std::string str3 = pop[p1].get_str().substr(cross_point2);

    child.set_str(str1 + str2 + str3);

    return child;
}

baldwin baldwin_vec::uniform_crossover(int p1, int p2)
{
    baldwin child;

    std::string str1 = pop[p1].get_str();
    std::string str2 = pop[p2].get_str();
    std::string str3;
    str3.resize(str1.size());

    for (int i = 0; i < str1.size(); i++){
        int j = rand() % 2;
        if (j == 0) str3[i] = str1[i];
        else        str3[i] = str2[i];
    }
    child.set_str(str3);

    return child;
}

double baldwin_vec::mean_correct_position_percentage()
{
    double sum = 0;
    for (int i = 0; i < pop.size(); i++){
        sum += pop[i].correct_position_percentage();
    }
    return sum / pop.size();
}

double baldwin_vec::mean_incorrect_position_percentage()
{
    double sum = 0;
    for (int i = 0; i < pop.size(); i++){
        sum += pop[i].incorrect_position_percentage();
    }
    return sum / pop.size();
}





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


int main()
{
        
    baldwin_vec pop;
    baldwin_vec buffer;
    baldwin_vec *population = &pop;
    baldwin_vec *buffer_ptr = &buffer;
    
    pop.init_pop();

    std::vector<double> correct;
    std::vector<double> incorrect;
    
    for (int i=0; i<GA_MAXITER; i++){

        population->sort_by_fitness();
        
        std::cout << "Generation: " << i+1 << " Fitness: " << population->get_best().get_fitness() << std::endl;
        std::cout << "Best individual: " << population->get_best().get_str() << std::endl;

        incorrect.push_back(population->mean_incorrect_position_percentage());
        write_numbers_to_file("incorrect.txt", incorrect);
        mate(*population, *buffer_ptr);
       
        swap(population, buffer_ptr);      
    }


    return 0;

    
}