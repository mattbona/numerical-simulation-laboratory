#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"

using namespace std;

int main(int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;

	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;


	int M = 1000000;			// Total number of throws
	int N = 100;                 	// Number of blocks
	int L = static_cast<int>(M/N);  // Number of throws in each block, please use for M a multiple of N
	double *random_vec = new double[M]();		// Define random vector

	for(int i=0; i < M; i++){	// Load the vector with random number distributed uniformly
		random_vec[i] = rnd.Rannyu();
	}

	// Estimate of I sampling the uniform distr in [0,1]
	double *average1 = new double[N]();		// Define average vector
	double *average_sqr1 = new double[N]();		// Define average squared vector
	
	for(int i=0; i < N; i++){       // Compute the average of my observable and the aveˆ2 to calculate the variance
		double sum = 0;
		for(int j=0; j < L; j++){
			int k = j + i*L;
			sum +=  M_PI/2*cos(M_PI/2*random_vec[k]);
		}
		average1[i] = sum/L;
		average_sqr1[i] = pow(average1[i],2);
	}
	
	prog_average_std_dev_block_method("data/EX02_1(1).dat", average1, average_sqr1, N);	

         // Estimate of I using importance sampling with w(x)=pi/2*sin(pi/2*x)
         double *average2 = new double[N]();             // Define average vector
         double *average_sqr2 = new double[N]();         // Define average squared vector
 
         for(int i=0; i < N; i++){       // Compute the average of my observable and the aveˆ2 to calculate the variance
                 double sum = 0;
                 for(int j=0; j < L; j++){
                         int k = j + i*L;
                         sum += ( (M_PI)/3 )*( cos(M_PI/2*random_vec[k])/(1 - pow(random_vec[k],2.)) );
                 }
                 average2[i] = sum/L;
                 average_sqr2[i] = pow(average2[i],2);
         }
 
         prog_average_std_dev_block_method("data/EX02_1(2).dat", average2, average_sqr2, N);


	return 0;
}
