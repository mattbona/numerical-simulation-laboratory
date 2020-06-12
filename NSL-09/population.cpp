/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "population.h"

using namespace std;

// Constructor
Population :: Population(){};
// Destructor
Population :: ~Population(){};
// Methods
std::vector<Chromosome> Population :: GetPopulation(){
        return population;
};

void Population :: InitializePopulation(int my_population_size, int my_number_of_cities,
                                        std::vector<city> *my_p_world, Random *my_p_rnd){
        population_size = my_population_size;
        // Initialize population with chromosomes that have
        // uniform distributed path permutation
        population.resize(population_size);
        for(int i=0; i< population_size; i++)
                population[i].InitializeChromosome(my_number_of_cities,
                                                   my_p_world,my_p_rnd);
};

void Population :: SortPopulation(){
        vector<double> distances(population_size);
        vector<int> sorted_indexes(population_size);
        vector<Chromosome> new_population(population_size);
        // Compute L1 distances
        for(int i=0; i<population_size; i++){
                distances[i]= population[i].GetPathL1Distance();
        }
        // Sort indexes according to the L1 distance of the single path
        sorted_indexes = get_sorted_indexes_vector(distances);
        // Order the population
        for(int i=0; i<population_size; i++){
                new_population[i] = population[sorted_indexes[i]];
        }
        population = new_population;
};

// Functions
vector<int> get_sorted_indexes_vector(const vector<double> &v){
        // Initialize original index locations
        vector<int> idx(v.size());
        iota(idx.begin(), idx.end(), 0);
        // Sort indexes based on comparing values in v
        sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

        return idx;
};
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
