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

using namespace std;

class baldwin
{
private:
    double fitness;
    std::string str; 

public:
    baldwin();
    double get_fitness() { return fitness; }
    double calc_fitness();
    void init_str();
};

baldwin::baldwin()
{
    str.resize(strlen(BALD_STRING));
}

double baldwin::calc_fitness()
{
    fitness = 0;
    for (size_t i = 0; i < str.size(); ++i)
    {
        if (str[i] == BALD_STRING[i])
        {
            fitness++;
        }
    }
    return fitness;
}

void baldwin::init_str(){
     // Random number generation setup
    random_device rd; // Obtain a random number from hardware
    mt19937 gen(rd()); // Seed the generator
    uniform_real_distribution<> dis(0.0, 1.0); // Define the range

    // Initialize the str with the same length as BOLD_SRING
    str.resize(strlen(BALD_STRING));

    for (size_t i = 0; i < str.size(); ++i)
    {
        double prob = dis(gen); // Generate a random number between 0 and 1

        if (prob < 0.50) 
        {
            str[i] = '?';
        }
        else if (prob < 0.75) 
        {
            str[i] = '1';
        }
        else 
        {
            str[i] = '0';
        }
    }
}

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


int main()
{
        
    baldwin_vec pop;
    baldwin_vec buffer;
    baldwin_vec *population = &pop;
    baldwin_vec *buffer_ptr = &buffer;
    
    pop.init_pop();
    
    for (int i=0; i<GA_MAXITER; i++){

        population->sort_by_fitness();
        
        std::cout << "Generation: " << i+1 << " Fitness: " << population->get_best().get_fitness() << std::endl;
        
        if (population->get_best().get_fitness()==0){
            break;
        }

        mate(*population, *buffer_ptr);
       
        swap(population, buffer_ptr);      
    }


    return 0;

    
}