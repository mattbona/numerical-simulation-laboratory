#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"
#include "asset.h"

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
	
	// Estimate European call-option and put-option price by sampling directly the final asset price S(T) for a GBM(mu,sigma)
	double mu=0, sigma=1, asset_t0=100, T=1;
	double gaussian_var = rnd.Gauss(mu,sigma);	
	GBM asset(mu,sigma);
	asset.SetAssetValue(asset_t0);
	asset.SetGaussianVar(gaussian_var);
	cout <<"S(0)= " << asset.GetAssetValue() << endl;
	asset.UpdateAssetValue(T);
	cout <<"After T = " << T << " S(T) = " << asset.GetAssetValue() << endl;
	return 0;
}
