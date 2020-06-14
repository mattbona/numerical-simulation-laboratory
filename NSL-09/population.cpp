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
// Copy constructor
Population :: Population(const Population& my_population){
        this->population_size = my_population.population_size;
        this->p_rnd = my_population.p_rnd;
        this->population = my_population.population;
        this->permutation_probability = my_population.permutation_probability;
        this->shift_probability = my_population.shift_probability;

};
// Methods
std::vector<Chromosome> Population :: GetPopulation(){
        return population;
};
void Population :: InitializePopulation(int my_population_size, int my_number_of_cities,
                                        std::vector<city> *my_p_world, Random *my_p_rnd){
        population_size = my_population_size;
        p_rnd = my_p_rnd;
        // Initialize population with chromosomes that have
        // uniform distributed path permutation
        population.resize(population_size);
        for(int i=0; i< population_size; i++)
                population[i].InitializeChromosome(my_number_of_cities,
                                                   my_p_world,p_rnd);
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
void Population :: SetMutationProbabilities(double my_permutation_probability,
                                            double my_shift_probability){
        permutation_probability = my_permutation_probability;
        shift_probability = my_shift_probability;

};
void Population :: MutateChromosomes(int k, int l){
        vector<Chromosome> mutated_chromosome(2);
        mutated_chromosome[0] = population[k];
        mutated_chromosome[1] = population[l];

        // Chromosome 1
        if(p_rnd->Rannyu() < permutation_probability)
                mutated_chromosome[0].PermutePath();
        if(p_rnd->Rannyu() < shift_probability)
                mutated_chromosome[0].ShiftPath();

        // Chromosome 2
        if(p_rnd->Rannyu() < permutation_probability)
                mutated_chromosome[1].PermutePath();
        if(p_rnd->Rannyu() < shift_probability)
                mutated_chromosome[1].ShiftPath();

        // Check if mutated paths fullfil the bonds
        mutated_chromosome[0].CheckPath();
        mutated_chromosome[1].CheckPath();

        population[k] = mutated_chromosome[0];
        population[l] = mutated_chromosome[1];
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
