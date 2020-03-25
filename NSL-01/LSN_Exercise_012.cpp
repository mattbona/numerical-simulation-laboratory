#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"
#include "central_limit.h"

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

	// Create and load a vector with M uniformly distr random number 
	int M = 1000000;			
	double *random_vec = new double[M]();	// Define random vector

	for(int i=0; i < M; i++){		// Load the vector
		random_vec[i] = rnd.Rannyu();
	}

	int cases = 2;
	central_limit uniform;
	uniform.set_uniform_cases(cases);

	cout << uniform.get_uniform_prob() << endl;
	
	return 0;
}
