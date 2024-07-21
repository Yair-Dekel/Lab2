#pragma warning(disable:4786)		// disable debug warning

#include <iostream>					// for cout etc.
#include <vector>					// for vector class
#include <string>					// for string class
#include <algorithm>				// for sort algorithm
#include <time.h>					// for random seed
#include <math.h>					// for abs()
#include <cmath>
#include <chrono>
#include <string.h>
#include <cstring>
#include <set>
#include <random>
#include <fstream>
#include <numeric>
#include <functional>

#include "t_functions.h"

/*#include <QtCharts>
#include <QChartView>
#include <QLineSeries>
#include <QApplication>*/


// #define GA_MAXITER		 200			// maximum iterations
// #define GA_ELITRATE		 0.05f		    // elitism rate
// #define GA_MUTATIONRATE	 0.4f		    // mutation rate
// #define GA_MUTATION		 RAND_MAX * GA_MUTATIONRATE
// #define GA_TARGET		 std::string("Hello world!")
// #define MAX_AGE          3
// #define REPRODUCE_THRESH 1
// #define BIN_MAX_CAPACITY 10

using namespace std;				// polluting global namespace, but hey...



/////**********************************************************************************/////

static int mutation_rate_flag;
static vector<double> preset_data;
static double value_to_insrt = 0.0 ;
static int value_to_insrtcounter = 0;
static int curr_fitness;
static int global_counter = 0;
static int alils_number = 0;

static int prev_best_fitness=1000;
static int curr_best_fitness=1000;

static double SF_GLOB; 


class sudoku
{
private:
    std::vector<int> board[9];
    bool steel_numbers[9][9]={false};
    unsigned int fitness;
    int age;
    
public:
    sudoku();
    sudoku(const sudoku &sud);
    sudoku& operator=(const sudoku &sud);
    ~sudoku();
    void set_board(std::vector<std::vector<int>> board);
    std::vector<std::vector<int>> get_board() const;
    void set_board(int i, int j, int value);
    int get_board_cell(int i, int j) const;
    void set_steel_numbers(int i, int j, bool value);
    bool get_steel_numbers(int i, int j) const;
    unsigned int get_fitness() const;
    void set_fitness(int k);
    int get_age() const;
    void set_age(int age);
    void calc_fitness();
    void print_board();
    void mutate();
    void inversion_mutation();
    void scramble_mutation();
    void simple_mutation();
    void swap_cell(int &a, int &b);
    void kill_and_replace();
    void init_board();
};

sudoku::sudoku()
{
    for(int i=0;i<9;i++){
        board[i].resize(9);   
    }
    fitness = 0;
    age = 0;
}

sudoku::sudoku(const sudoku &sud)
{

    for(int i=0;i<9;i++){
        this->board[i].resize(9);
        for(int j=0;j<9;j++){
            this->board[i][j] = sud.get_board_cell(i,j);
            this->steel_numbers[i][j] = sud.get_steel_numbers(i,j);
        }
    }
    this->fitness = sud.get_fitness();
    this->age = sud.age;
}

sudoku& sudoku::operator=(const sudoku &sud)
{
    if (this == &sud){
        return *this;
    }
    for(int i=0;i<9;i++){
        this->board[i].resize(9);
        for(int j=0;j<9;j++){
            this->board[i][j] = sud.get_board_cell(i,j);
            this->steel_numbers[i][j] = sud.get_steel_numbers(i,j);
        }
    }
    this->fitness = sud.get_fitness();
    this->age = sud.age;
    return *this;
}

sudoku::~sudoku()
{
}

void sudoku::set_board(std::vector<std::vector<int>> board)
{
    for(int i=0;i<9;i++){
        for(int j=0;j<9;j++){
            this->board[i][j] = board[i][j];
            if (this->board[i][j] != 0){
				this->steel_numbers[i][j] = true;
			}
        }
    }
}

std::vector<std::vector<int>> sudoku::get_board() const
{
    std::vector<std::vector<int>> res(9, std::vector<int>(9, 0)); // Initialize res with correct size and initial value
    for(int i=0; i<9; i++){
        for(int j=0; j<9; j++){
            res[i][j] = board[i][j];
        }
    }
    return res;
}

void sudoku::set_board(int i, int j, int value)
{
    board[i][j] = value;
}

int sudoku::get_board_cell(int i, int j) const
{
    return board[i][j];
}

void sudoku::set_steel_numbers(int i, int j, bool value)
{
    steel_numbers[i][j] = value;
}

bool sudoku::get_steel_numbers(int i, int j) const
{
    return steel_numbers[i][j];
}

unsigned int sudoku::get_fitness() const
{
    return fitness;
}

void sudoku::set_fitness(int k){
    this->fitness = k;
}

int sudoku::get_age() const
{
    return age;
}

void sudoku::set_age(int age)
{
    this->age = age;
}

void sudoku::calc_fitness()
{
    fitness = 0;
    for (int i = 0; i < 9; i++){
        std::set<int> col = {1,2,3,4,5,6,7,8,9};
        for (int j = 0; j < 9; j++){
            col.erase(board[j][i]);
        }
        fitness += col.size();
    }
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            std::set<int> square = {1,2,3,4,5,6,7,8,9};
            for (int k = 0; k < 3; k++){
                for (int l = 0; l < 3; l++){
                    square.erase(board[i*3+k][j*3+l]);
                }
            }
            fitness += square.size();
        }
    }
}

void sudoku::print_board()
{
    for (int i = 0; i < 9; i++){
        for (int j = 0; j < 9; j++){
            std::cout << board[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void sudoku::mutate()
{
    int i = rand() % 3;
    if (i == 0)
        this->inversion_mutation();
    else if (i == 1)
        this->scramble_mutation();
    else
        this->simple_mutation();
}

void sudoku::inversion_mutation()
{
    int tsize = 9;
	int line = rand() % tsize; //choose a random line
	int count = 0;             //count the number of cells that are not steel numbers
	std::vector<int> indexes;  //store the indexes of cells that are not steel numbers
	for(int i=0; i<tsize; i++){
		if(!this->get_steel_numbers(line, i)){
			count++;
			indexes.push_back(i);
		}
	}

	if(count < 2){ //if there are less than 2 cells that are not steel numbers, return
		return;
	}
	int invers_num = (rand() % count) + 1; //choose a random number of cells to invert
	int start = rand() % (count - invers_num + 1); //choose a random starting index
	for(int i=0; i<invers_num/2; i++){

		this->swap_cell(this->board[line][indexes[start+i]], this->board[line][indexes[start+invers_num-i-1]]);
    }
}

void sudoku::scramble_mutation()
{
	int tsize = 9;
	int line = rand() % tsize;
	int count = 0;
	std::set<int> indexes;
	std::set<int> values;
	for(int i=0; i<tsize; i++){
		if(!this->steel_numbers[line][i]){
			count++;
			indexes.insert(i);
			values.insert(this->board[line][i]);
		}
	}
	if(count < 2){
		return;
	}
	while (indexes.size() > 0){
		int rand_num = rand() % indexes.size();
		auto it = indexes.begin();
		for (int i = 0; i < rand_num; i++){
			it++;
		}
		int rand_num2 = rand() % values.size();
		auto it2 = values.begin();
		for (int i = 0; i < rand_num2; i++){
			it2++;
		}
		this->board[line][*it] = *it2;
		indexes.erase(it);
		values.erase(it2);
	}
}    

void sudoku::swap_cell(int &a, int &b)
{
    int temp = a;
    a = b;
    b = temp;
}

void sudoku::simple_mutation()
{
    int tsize = 9;
    int line = rand() % tsize;
    int j = rand() % tsize;
    int i = rand() % tsize;

    while(this->get_steel_numbers(line,i) || this->get_steel_numbers(line,j)){
        j = rand() % tsize;
        i = rand() % tsize;
        line = rand() % tsize;
    }
    this->swap_cell(this->board[line][i], this->board[line][j]);

}

void sudoku::kill_and_replace()
{
    int i = rand() % 3;
    if (i == 0)
        this->init_board();
    else 
        this->mutate();
    this->set_age(0);
    this->calc_fitness();
}

void sudoku::init_board()
{
    int k=0;   
	for (int i = 0; i < 9; i++){
		std::set<int> row = {1,2,3,4,5,6,7,8,9};
		for (int j = 0; j < 9; j++){
			if (this->get_steel_numbers(i, j)){
				row.erase(this->board[i][j]);
			}
		}
		for (int j = 0; j < 9; j++){
			if (!this->get_steel_numbers(i, j)){
				int rand_num = rand() % row.size();
				auto it = row.begin();
				for (int l = 0; l < rand_num; l++){
					it++;
				}
				this->set_board(i, j, *it);
				row.erase(it);
			}
		}
	}
}


/////**********************************************************************************/////

//population

class sudoku_vec
{
private:
    std::vector<sudoku> vec;
    std::vector<int> linear_scaled_fitnesses;
    flags flag;
    int generation_number;

public:
    sudoku_vec(sudoku s, int size);
    sudoku_vec(int size);
    ~sudoku_vec();
    unsigned int get_all_fitnesses();
    void sort_by_fitness();
    sudoku get_best();
    sudoku mate(int parent1, int parent2);
    void set(int i, sudoku s);
    sudoku get(int i);

/************************************************************** */
    void set_generation_number(int GN);
    double calc_relative_fitness(const sudoku& individual);
    double set_mutation_probability_individual(sudoku& individual);
    double calc_SF();
/************************************************************** */

    sudoku crossover_PMX(int parent1, int parent2);
    sudoku crossover_CX(int parent1, int parent2);
    void increase_age(int i);
    
    void parent_selection_RWS(int &i1);
    void parent_selection_SUS(vector<int> &parents, int parents_num);  
    void tournament_selection(int K, double P, int &selected_index);  
    void liniar_scaling();

    double calc_average_gene_distance();
    int calc__number_of_different_alleles();
    int calc_fitness_sum();
    void set_flags(flags flag);
    const flags& get_flags() const;
    double calc_fitness_variance(double mean);

};

void sudoku_vec::set_generation_number(int GN){
    generation_number = GN;
}


sudoku_vec::sudoku_vec(sudoku sud, int size)
{
    vec.resize(size);
    linear_scaled_fitnesses.resize(size);
    std::cout << "Constructor size: " << size << std::endl;
	
    for (int i=0; i<size; i++) {
		sudoku citizen;
		
		int k=0;
        
        std::vector<std::vector<int>> sud_board = sud.get_board();
		for (int i = 0; i < 9; i++){
			std::set<int> row = {1,2,3,4,5,6,7,8,9};
			for (int j = 0; j < 9; j++){
				if (sud_board[i][j] != 0){
					row.erase(sud_board[i][j]);
					citizen.set_board(i, j, sud_board[i][j]);
				}
			}

			for (int j = 0; j < 9; j++){
				if (citizen.get_board_cell(k, j) == 0){
					int rand_num = rand() % row.size();
					auto it = row.begin();
					for (int l = 0; l < rand_num; l++){
						it++;
					}
					citizen.set_board(k, j, *it);
					row.erase(it);
				}
			}
			k++;
		}
		for(int i=0; i<9; i++){
			for(int j=0; j<9; j++){
				citizen.set_steel_numbers(i,j ,sud.get_steel_numbers(i,j));
			}
		}
        citizen.calc_fitness();
		vec[i] = citizen;
	}
    
    liniar_scaling();
	
}

sudoku_vec::sudoku_vec(int size)
{
    vec.resize(size);
    linear_scaled_fitnesses.resize(size);

}

sudoku_vec::~sudoku_vec()
{
}

unsigned int sudoku_vec::get_all_fitnesses()
{
    unsigned int fitness = 0;
    for (int i = 0; i < vec.size(); i++){
        fitness += vec[i].get_fitness();
    }
    return fitness;
}

void sudoku_vec::sort_by_fitness()
{
    std::sort(vec.begin(), vec.end(), [](sudoku a, sudoku b) { return a.get_fitness() < b.get_fitness(); });
}

sudoku sudoku_vec::get_best()
{
    return vec[0];
}



//--------------------------------------------------------------------------------------------------------------


/// the relative fitness calculation
double sudoku_vec::calc_relative_fitness(const sudoku& individual) {
    int fitness_sum = calc_fitness_sum();
    int count = vec.size();
    
    double mean = static_cast<double>(fitness_sum) / count;
    double Rf = static_cast<double>(individual.get_fitness()) / mean;
    
    /// normalize
    
    double SF = 0;
    
    SF = SF_GLOB;

    double norm_Rf = Rf / SF;
    
    return norm_Rf;

}

double sudoku_vec::calc_SF(){
    int fitness_sum = calc_fitness_sum();
    int count = vec.size();
    double mean = static_cast<double>(fitness_sum) / count;
    double SF = 0;
    for(sudoku itm: vec){
        SF = SF + (itm.get_fitness()/ mean);
    } 
    return SF;
}


/// based on the relative fitness we are setting the mutation probability value 
double sudoku_vec::set_mutation_probability_individual(sudoku& individual) {
    
    double max_mutation_probability = GA_MUTATIONRATE;
    double relative_fitness = calc_relative_fitness(individual);
    
    double individual_mutation_probability = max_mutation_probability * (1.0 - relative_fitness);
    
    
    return individual_mutation_probability;

}


//--------------------------------------------------------------------------------------------------------------

sudoku sudoku_vec::mate(int parent1, int parent2)
{

    int i1 = rand() % 2;

    sudoku child;
    if (i1 == 0)
        child = crossover_CX(parent1,parent2);
    else
        child = crossover_PMX(parent1,parent2);
    
    
    if(mutation_rate_flag == 1){

        if (rand() < GA_MUTATION && flag.mutation_F){
            child.mutate();
        }

    }
    else if(mutation_rate_flag == 2){
        
        double mutation_rate =  Non_Unform_Mutation(generation_number);
        

        int mutation_threshold = mutation_rate * RAND_MAX;

        if (rand() < mutation_threshold && flag.mutation_F){ 
            child.mutate();
        }
    
    }
    else if(mutation_rate_flag == 3){
        
        int alleles = alils_number;
        

        int threshold_smaller = 500;
        int threshold_biger = 700;
        double mutation_rate = GA_MUTATIONRATE;
        if(alleles < threshold_smaller){
            mutation_rate = 0.8;
        }
        else if(alleles > threshold_biger){
            mutation_rate = 0.2;

        }

        int mutation_threshold = mutation_rate * RAND_MAX;
        
        if (rand() < mutation_threshold && flag.mutation_F){
            child.mutate();
        } 

    }
    else if(mutation_rate_flag == 4){
          
        double mutation_rate = GA_MUTATIONRATE;

        if(global_counter >= 5){

            mutation_rate = 0.9f;
            

            if((global_counter > 25)){
                
                mutation_rate = GA_MUTATIONRATE;
                global_counter = 0;
            }
        }

        int mutation_threshold = mutation_rate * RAND_MAX;
        
        if (rand() < mutation_threshold && flag.mutation_F){
           
            child.mutate();
        } 
        
    }
    else if(mutation_rate_flag == 5){
        double prob_mutation = set_mutation_probability_individual(child);
        double mutation_threshold = prob_mutation * RAND_MAX;

        
        if (rand() < mutation_threshold && flag.mutation_F) 
            child.mutate();  

    }

    return child;
}




void sudoku_vec::set(int i, sudoku s)
{
    vec[i] = s;
}

sudoku sudoku_vec::get(int i)
{
    return vec[i];
}

sudoku sudoku_vec::crossover_PMX(int i1, int i2)
{
    int tsize = 9;
    
    
    while (vec[i1].get_age() < REPRODUCE_THRESH && vec[i1].get_age()!=0){
        i1 = rand() % GA_POPSIZE/2;
    }
    while (vec[i2].get_age() < REPRODUCE_THRESH && vec[i2].get_age()!=0){
        i2 = rand() % GA_POPSIZE/2;
    }
    sudoku parent1 = vec[i1];
    sudoku parent2 = vec[i2];
    
    
    sudoku child;

    child = parent1;
	
	for(int i=0; i<tsize; i++){
		int spos = rand() % tsize;
        int tries = 20;
		while (parent1.get_steel_numbers(i,spos)&& tries--)
		{
			spos = rand() % tsize;
		}
		int first = parent1.get_board_cell(i,spos); // Save the value to be replaced
		int second = parent2.get_board_cell(i,spos); // Save the value to be replaced
		for(int j=0; j<tsize; j++){
			if(parent1.get_board_cell(i,j) == second && !parent1.get_steel_numbers(i,j)){ 
                child.set_board(i,spos,second); 
                child.set_board(i,j,first); 
                break;
            }
        }
	
    }
    child.calc_fitness();
    return child;
}

sudoku sudoku_vec::crossover_CX(int parent1_ind, int parent2_ind)
{   

    int tsize = 9;
    int i1 = parent1_ind;
    int i2 = parent2_ind;

    sudoku parent1 = vec[i1];
    sudoku parent2 = vec[i2];
    sudoku offspring = parent1;

	for (int i=0; i<tsize; i++){
		int random = rand() % tsize;
        int tries = 20;
		while (parent1.get_steel_numbers(i,random) && tries--)//choose a random number that is not steel number
		{
			random = rand() % tsize;
		}
		
		int start = parent1.get_board_cell(i,random);
		while(parent2.get_board_cell(i,random) != start){
			offspring.set_board(i, random, parent2.get_board_cell(i,random));
			for(int j=0; j<tsize; j++){
				if(parent1.get_board_cell(i,j) == offspring.get_board_cell(i,random)){
					random = j;
					break;
				}
			}
		}
		offspring.set_board(i,random, parent2.get_board_cell(i,random));
	}
	
    offspring.calc_fitness();
    return offspring;
}

void sudoku_vec::increase_age(int i){
    
    vec[i].set_age(vec[i].get_age() + 1); //age++
    if(vec[i].get_age() > MAX_AGE){
        vec[i].kill_and_replace();
    }
}

void sudoku_vec::liniar_scaling(){
    unsigned int max_value = 0;
    unsigned int min_value = numeric_limits<unsigned int>::max(); // Initialize to maximum possible value
    unsigned int sum = 0;

    // Find max, min, and sum of fitness values
    int i=1;
    for(const auto &ind : vec){
        if (ind.get_fitness() > max_value) {
            max_value = ind.get_fitness();
        }
        if (ind.get_fitness() < min_value) {
            min_value = ind.get_fitness();
        }
        sum += ind.get_fitness();
        i++;
    }

    // Calculate average fitness
    unsigned int avrg = sum / vec.size();

    // Calculate scaling parameters
    float a = static_cast<float>(1 - i) / (max_value - min_value);
    float b = 1 - a * avrg;

    // Scale fitness values and create a new scaled population
    for(int i = 0; i < vec.size(); i++){
        linear_scaled_fitnesses[i] = a * vec[i].get_fitness() + b;
    }
    
}

void sudoku_vec::parent_selection_RWS(int &i1) {
	    
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

void sudoku_vec::parent_selection_SUS(vector<int> &parents, int parents_num) {
	int tsize = 9;
	int max = vec[GA_POPSIZE-1].get_fitness();
	max++; 				// Avoiding fitness = 0
	int total_fitness = 0;
	for (const auto &element : vec) { // Calculate the total fitness
		total_fitness += max-element.get_fitness();
	}
	// Select points
	int r1 = rand() % total_fitness;
	vector<int> points;

	for(int i = 0; i < parents_num; i++){
		int point = (r1 + i * total_fitness / parents_num) % total_fitness;
		points.push_back(point);
	}
	sort(points.begin(), points.end());

	// Select parents by accumulated fitness
	int current_member = 0;
	int comulative_sum = 0;
	for (int point: points) {
		while(comulative_sum + max - vec[current_member].get_fitness() < point){
			comulative_sum += max - vec[current_member].get_fitness();
			current_member++;
		}
		parents.push_back(current_member);
	}
}

void sudoku_vec::tournament_selection(int K, double P, int &selected_index) {
    int tournament_winner_index = 0;
    unsigned int best_fitness = std::numeric_limits<unsigned int>::max(); ///Initialize to maximum possible value
    ///Perform K tournaments
    for (int i = 0; i < K; ++i) {        
        ///Randomly select an individual from the population
        int index = rand() % vec.size();
        const sudoku& competitor = vec[index];

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

double sudoku_vec::calc_average_gene_distance(){
    std::vector<int> pairwise_distances;
    int pair_count = 0;    
    
    for (int i = 0; i < vec.size(); i++){
        for (int j = i + 1; j < vec.size(); j++){
            pair_count++;
            for (int k = 0; k < 9; k++){
                for(int l=0; l<9; l++){
                    pairwise_distances.push_back(std::abs(vec[i].get_board_cell(k,l) - vec[j].get_board_cell(k,l)));
                }
            }
        }
    }
    
    return std::accumulate(pairwise_distances.begin(), pairwise_distances.end(), 0.0) / pair_count;
}

int sudoku_vec::calc__number_of_different_alleles() {
    
    std::set< std::vector<std::vector<int>> > alleles;

    for (const auto& individual : vec) {
        alleles.insert(individual.get_board());
    }

    return alleles.size();

}

int sudoku_vec::calc_fitness_sum()
{
    int total_fitness = std::accumulate(vec.begin(), vec.end(), 0, [](int sum, const sudoku& b) {
    return sum + b.get_fitness();
    });
    return total_fitness;
}

void sudoku_vec::set_flags(flags flag)
{
    this->flag = flag;
}

const flags& sudoku_vec::get_flags() const
{
    return flag;
}

double sudoku_vec::calc_fitness_variance(double mean)
{
    double fitness_variance = std::accumulate(vec.begin(), vec.end(), 0, [mean](double sum, const sudoku& b) {
    return sum + (b.get_fitness() - mean) * (b.get_fitness() - mean);
    });
    return fitness_variance;
}

/////**********************************************************************************/////

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





int main(int argc, char* argv[])
{
    flags flag;
    int i = 1;
    int index = 3;
    int pop_size = GA_POPSIZE;
    
    
    while (i < argc){
        if (strcmp(argv[i], "-RWS") == 0){
            flag.RWS_F = true;
            flag.SUS_F = false;
            flag.tournament_F = false;
            flag.linear_scaling_F = true;
        }
        else if (strcmp(argv[i], "-SUS") == 0){
            flag.SUS_F = true;
            flag.RWS_F = false;
            flag.tournament_F = false;
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
        else if (strcmp(argv[i], "-pop_size") == 0){
            pop_size = atoi(argv[i+1]);
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
    
    
    
    std::vector<std::vector<int>> board1 = {
    {0,0,0,2,6,0,7,0,1},
	{6,8,0,0,7,0,0,9,0},
	{1,9,0,0,0,4,5,0,0},
	{8,2,0,1,0,0,0,4,0},
	{0,0,4,6,0,2,9,0,0},
	{0,5,0,0,0,3,0,2,8},
	{0,0,9,3,0,0,0,7,4},
	{0,4,0,0,5,0,0,3,6},
	{7,0,3,0,1,8,0,0,0}
    };

    std::vector<std::vector<int>> board2 = {
    {1,0,0,4,8,9,0,0,6},
	{7,3,0,0,7,0,0,4,0},
	{0,0,0,0,0,1,2,9,5},
	{0,0,7,1,2,0,6,0,0},
	{5,0,0,7,0,3,0,0,8},
	{0,0,6,0,9,5,7,0,0},
	{9,1,4,6,0,0,0,0,0},
	{0,2,0,0,0,0,0,3,7},
	{8,0,0,5,1,2,0,0,4}
    };

    std::vector<std::vector<int>> board3 = {
	{0,2,0,6,0,8,0,0,0},
	{5,8,0,0,0,9,7,0,0},
	{0,0,0,0,4,0,0,0,0},
	{3,7,0,0,0,0,5,0,0},
	{6,0,0,0,0,0,0,0,4},
	{0,0,8,0,0,0,0,1,3},
	{0,0,0,0,2,0,0,0,0},
	{0,0,9,8,0,0,0,3,6},
	{0,0,0,3,0,6,0,9,0}
    };

    std::vector<std::vector<int>> board4 = {
	{0,0,0,6,0,0,4,0,0},
	{7,0,0,0,0,3,6,0,0},
	{0,0,0,0,9,1,0,8,0},
	{0,0,0,0,0,0,0,0,0},
	{0,5,0,1,8,0,0,0,3},
	{0,0,0,3,0,6,0,4,5},
	{0,4,0,2,0,0,0,6,0},
	{9,0,3,0,0,0,0,0,0},
	{0,2,0,0,0,0,1,0,0}
    };

    std::vector<std::vector<int>> board5 = {
	{2,0,0,3,0,0,0,0,0},
	{8,0,4,0,6,2,0,0,3},
	{0,1,3,8,0,0,2,0,0},
	{0,0,0,0,2,0,3,9,0},
	{5,0,7,0,0,0,6,2,1},
	{0,3,2,0,0,6,0,0,0},
	{0,2,0,0,0,9,1,4,0},
	{6,0,1,2,5,0,8,0,9},
	{0,0,0,0,0,1,0,0,2}
    };

    std::vector<std::vector<int>> board6 = {
    {0,2,0,0,0,0,0,0,0},
	{0,0,0,6,0,0,0,0,3},
    {0,7,4,0,8,0,0,0,0},
    {0,0,0,0,0,3,0,0,2},
    {0,8,0,0,4,0,0,1,0},
    {6,0,0,5,0,0,0,0,0},
    {0,0,0,0,1,0,7,8,0},
    {5,0,0,0,0,9,0,0,0},
    {0,0,0,0,0,0,0,4,0}
    };
    
    std::vector<std::vector<int>> board;
    if (index == 1){
        board = board1;
    }
    else if (index == 2){
        board = board2;
    }
    else if (index == 3){
        board = board3;
    }
    else if (index == 4){
        board = board4;
    }
    else if (index == 5){
        board = board5;
    }
    else if (index == 6){
        board = board6;
    }
    else{
        std::cout << "Invalid index" << std::endl;
        return 1;
    }
    
    std::vector<double> average_gene_distance;
    std::vector<double> variance;
    std::vector<double> distance;
    std::vector<int> different_alleles;
    std::vector<int> best_fitness_vec;

    sudoku s;
    s.set_board(board);
    sudoku_vec sv(s, pop_size);
    sudoku_vec buffer(pop_size);
    sudoku_vec *population = &sv;
    sudoku_vec *buffer_ptr = &buffer;
    population->set_flags(flag);
    buffer_ptr->set_flags(flag);

    //******************************************************************************* */
    int mutation_rate_flag_;
    std::cout << "Enter mutation rate flag number: " << endl;
    std::cout << "1. basic" << endl;
    std::cout << "2. Non-Unform" << endl;
    std::cout << "3. Adaptive" << endl;
    std::cout << "4. Triggered Hyper" << endl;
    std::cout << "5. Self-Adaptive" << endl << endl;    
    std::cin >> mutation_rate_flag_;

    std::cout << "flag is " << mutation_rate_flag_ << endl; 
    //******************************************************************************* */
    
    for (int i=0; i<  GA_MAXITER; i++){
        population->sort_by_fitness();
        if(flag.linear_scaling_F)population->liniar_scaling();
        
        std::cout << "Generation: " << i+1 << " Fitness: " << population->get_best().get_fitness() << std::endl;
        population->get_best().print_board();
        population->set_generation_number(i+1);

        top_average_and_variance(*population, average_gene_distance, variance);
        
        //different_alleles.push_back(population->calc__number_of_different_alleles());
        //distance.push_back(population->calc_average_gene_distance());
        //best_fitness_vec.push_back(population->get_best().get_fitness());
        //cout << "Average gene distance: " << population->calc_average_gene_distance() << endl;
        //cout << "Number of different alleles: " << population->calc__number_of_different_alleles() << endl;
        
        curr_best_fitness = population->get_best().get_fitness();

            
        if(prev_best_fitness == curr_best_fitness){
            global_counter = global_counter + 1;
        }
        
        else{
            global_counter = 0;
        }


        if (population->get_best().get_fitness() == 0){
            break;
        }
        mutation_rate_flag = mutation_rate_flag_;
        
        alils_number = population->calc__number_of_different_alleles();
        SF_GLOB = population->calc_SF();

        //mate(*population, *buffer_ptr);
        non_deternistic_crowding(*population, *buffer_ptr, population->get_best());

        preset_data.push_back(population->get_best().get_fitness());
        
        prev_best_fitness = curr_best_fitness;   

        swap(population, buffer_ptr);
        
                
    }

    //write_numbers_to_file("average_gene_distance_sudoku.txt", average_gene_distance);
    //write_numbers_to_file("variance_sudoku.txt", variance);
    //write_numbers_to_file("different_alleles_sudoku.txt", different_alleles);
    //write_numbers_to_file("distance_sudoku.txt", distance);
    //write_numbers_to_file("best_fitness_sudoku.txt", best_fitness_vec);
    write_numbers_to_file("Non_Unform_mutation.txt", preset_data);


    return 0;
}