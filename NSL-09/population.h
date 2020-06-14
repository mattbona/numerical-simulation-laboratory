/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __POPULATION__
#define __POPULATION__

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort
#include "random.h"
#include "city.h"
#include "chromosome.h"

class Population{
private:
        int population_size;
        Random* p_rnd;
        std::vector<Chromosome> population;
        double permutation_probability;
        double block_permutation_probability;
        double shift_probability;
        double partial_shift_probability;
        double inversion_probability;
        double crossover_probability;
public:
        // Constructor
        Population();
        // Destructor
        ~Population();
        // Copy constructor
        Population(const Population&);
        // Methods
        std::vector<Chromosome> GetPopulation();
        void SetPopulation(std::vector<Chromosome>);
        void InitializePopulation(int, int, std::vector<city>*, Random*);
        void SortPopulation();
        void SetMutationProbabilities(double, double, double, double, double);
        std::vector<Chromosome> GetMutatedChromosomes(Chromosome,Chromosome);
        void SetCrossoverProbability(double);
        std::vector<Chromosome> GetCrossoveredChromosomes(int, int);
};
//Functions
std::vector<int> get_sorted_indexes_vector(const std::vector<double> &);

#endif
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
