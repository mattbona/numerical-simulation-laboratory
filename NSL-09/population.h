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
#include "random.h"
#include "city.h"
#include "chromosome.h"

class Population{

private:
        int population_size;
        std::vector<Chromosome> population;
public:
        // Constructor
        Population();
        // Destructor
        ~Population();
        // Methods
        void InitializePopulation(int, int, std::vector<city>*, Random*);
        std::vector<Chromosome> GetPopulation();

};

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
