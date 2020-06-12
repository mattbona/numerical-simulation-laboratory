/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "LSN_Exercise_091.h"
#include "chromosome.h"
#include "population.h"

using namespace std;

int main()
{
        // Inizialization
        Input();


        Population path_population;
        path_population.InitializePopulation(population_size, number_of_cities,
                                             &world, &rnd);

        for(int j=0; j<population_size; j++){
                for(int i=0; i<number_of_cities; i++)
                        cout << path_population.GetPopulation()[j].GetPath()[i] << endl;
                cout << path_population.GetPopulation()[j].GetPathL2Distance() << endl;
                cout << endl;
        }

        return 0;
}

//##############################################################################
void Input(void){
        int seed[4];
        int p1, p2;
        ifstream Primes("Primes");
        if (Primes.is_open()){
            Primes >> p1 >> p2 ;
        } else cerr << "PROBLEM: Unable to open Primes" << endl;
        Primes.close();

        ifstream input("seed.in");
        if (input.is_open()){
                    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                    rnd.SetRandom(seed,p1,p2);
                    input.close();
        } else cerr << "PROBLEM: Unable to open seed.in" << endl;

        cout << "Solving the Traveling Salesman Problem (TSP)" << endl;
        cout << "via a Genetic Algorithm.                    " << endl << endl;

        ifstream ReadInput, ReadConfig;
        //Read informations from input.dat
        ReadInput.open("input.dat");
        ReadInput >> number_of_generations;
        ReadInput >> population_size;

        // Read city positions from world_config.dat
        ReadConfig.open("world_config.dat");
        ReadConfig >> number_of_cities;
        world.resize(number_of_cities);
        for(int i=0; i<number_of_cities; i++){
                ReadConfig >> world[i].x;
                ReadConfig >> world[i].y;
        }

        cout << "Simulation parameters: " << endl;
        cout << "Number of generations = " << number_of_generations << endl;
        cout << "Population size = " << population_size << endl;
        cout << "Number of cities = " << number_of_cities << endl<<endl;

        ReadInput.close();
        ReadConfig.close();
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
