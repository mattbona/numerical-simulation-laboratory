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

std::vector<Chromosome> Population :: GetPopulation(){
        return population;
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
