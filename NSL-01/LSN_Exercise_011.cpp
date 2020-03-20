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


	int M = 100000;			// Total number of throws
	int N = 100;                 	// Number of blocks
	int L = static_cast<int>(M/N);  // Number of throws in each block, please use for M a multiple of N
	double *random_vec = new double[M]();		// Define random vector
	double *average = new double[N]();		// Define average vector
	double *average_sqr = new double[N]();		// Define average squared vector
	double *prog_average = new double[N]();		// Define progressive average vector
	double *prog_average_sqr = new double[N]();	// Define progressive average squared vector
	double *prog_error = new double[N]();		// Define progressive error vector

	ofstream out_file;
	out_file.open("ex1.dat");

	for(int i=0; i < M; i++){	// Load the vector with random number distributed uniformely
		random_vec[i] = rnd.Rannyu();
	}

	for(int i=0; i < N; i++){       // Compute the average of my observable and the aveˆ2 to calculate the variance
		double sum = 0;
		for(int j=0; j < L; j++){
			int k = i*j + L;
			sum = sum + random_vec[k];
		}
		average[i] = sum/L;
		average_sqr[i] = pow(average[i],2);
	}
	
	for(int i=0; i < N; i++){       // Compute and save the progressive average and error over the blocks 

		for(int j=0; j < i+1; j++){
			prog_average[i] += average[j];
			prog_average_sqr[i] += average_sqr[j];
		}

		prog_average[i] = prog_average[i]/(i+1);
		prog_average_sqr[i] = prog_average_sqr[i]/(i+1);
		prog_error[i] = std_dev(prog_average[i],prog_average_sqr[i],i);
		out_file << prog_average[i] << "\t" << prog_error[i] << endl;
	}

	out_file.close();
	return 0;
}
