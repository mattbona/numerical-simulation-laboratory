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

	// Create and load a vector with M uniformly distr random number
        int M = 1E4;			// Total number of throws
	int N = 2;

	double exponential_decay_rate = 1;
        double lorentzian_mean = 0;
        double lorentzian_width = 1;

        double *average_standard_dice = new double[M]();
        double *average_exponential_dice = new double[M]();
        double *average_lorentzian_dice = new double[M]();

        ofstream out_file;
        out_file.open("results/EX012_N2.dat");
	for(int i=0; i < M; i++){
		double sum_standard_dice = 0;
                double sum_exponential_dice = 0;
                double sum_lorentzian_dice = 0;
		for(int j=0; j < N; j++){
			sum_standard_dice +=  rnd.Rannyu();
                        sum_exponential_dice += rnd.Exponential(exponential_decay_rate);
                        sum_lorentzian_dice += rnd.Cauchy(lorentzian_width,lorentzian_mean);
                }
		average_standard_dice[i] = sum_standard_dice/N;
                average_exponential_dice[i] = sum_exponential_dice/N;
                average_lorentzian_dice[i] = sum_lorentzian_dice/N;
                out_file << average_standard_dice[i] << " " << average_exponential_dice[i] << " " << average_lorentzian_dice[i] << endl;
	}

        out_file.close();
	return 0;
}
