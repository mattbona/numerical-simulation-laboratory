/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __CHROMOSOME__
#define __CHROMOSOME__

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "city.h"

class Chromosome{

private:
        int number_of_cities;
        std::vector<city> *p_world;
        std::vector<int> path;
        Random *p_rnd;
public:
        // Constructor
        Chromosome();
        // Destructor
        ~Chromosome();
        // Methods
        void InitializeChromosome(int, std::vector<city> *, Random *);
        void PermutePath();
        void CheckPath();
        std::vector<int> GetPath();
        double GetPathL1Distance();
        double GetPathL2Distance();

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
