#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"
#include "asset.h"
#include "option.h"

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
    int block_number = 1E2;
    int iteration_per_block = 1E4;
    int tot_random_number = block_number*iteration_per_block;

    double risk_free_rate = 0.1;
    double drift = risk_free_rate, volatility = 0.25;
    double asset_initial_price = 100;
    double expire_date = 1;
    double initial_date = 0;
    double strike_price = 100;

    GBM asset(drift,volatility);
    double *random_gauss_vec = new double[tot_random_number]();		// Define random vector
	for(int i=0; i <  tot_random_number; i++){	// Load the vector with random number distributed uniformly
		random_gauss_vec[i] = rnd.Gauss(0,1);
	}

    European call_european;
    call_european.SetStrikePrice(strike_price);
    European put_european;
    put_european.SetStrikePrice(strike_price);

/*
double sum=0;
for(int i=10;i<1E5; i++){
    asset.SetAssetPrice(asset_initial_price);
    asset.SetGBMGaussianVar(random_gauss_vec[i]);
    asset.UpdateAssetPrice(initial_date, expire_date);
    put_european.SetAssetPrice(asset.GetAssetPrice());
    put_european.UpdatePutOptionProfit();
    sum+=exp(-1*risk_free_rate*expire_date)*put_european.GetPutOptionProfit();
    cout <<"Asset price:\t"<<asset.GetAssetPrice() <<"\tOption price:\t"<<exp(-1*risk_free_rate*expire_date)*put_european.GetPutOptionProfit() << endl;
};
    cout<<"Average option price:\t" << sum/1E5 << endl;
*/

    double *average1 = new double[block_number]();		// Define average vector
    double *average_sqr1 = new double[block_number]();		// Define average squared vector
    double *average2 = new double[block_number]();		// Define average vector
    double *average_sqr2 = new double[block_number]();		// Define average squared vector
    for(int i=0; i < block_number; i++){       // Compute the average of my observable and the aveˆ2 to calculate the variance
        double sum1 = 0;
        double sum2 = 0;
        for(int j=0; j < iteration_per_block; j++){
            int k = j + i*iteration_per_block;
            asset.SetGBMGaussianVar(random_gauss_vec[k]);

            asset.SetAssetPrice(asset_initial_price);
            asset.UpdateAssetPrice(initial_date, expire_date);

            call_european.SetAssetPrice(asset.GetAssetPrice());
            call_european.UpdateCallOptionProfit();
            sum1 += exp(-1*risk_free_rate*expire_date)*call_european.GetCallOptionProfit();
            put_european.SetAssetPrice(asset.GetAssetPrice());
            put_european.UpdatePutOptionProfit();
            sum2 += exp(-1*risk_free_rate*expire_date)*put_european.GetPutOptionProfit();
        }
        average1[i] = sum1/iteration_per_block;
        average_sqr1[i] = pow(average1[i],2);

        average2[i] = sum2/iteration_per_block;
        average_sqr2[i] = pow(average2[i],2);
    }

    prog_average_std_dev_block_method("data/EX03_1(1).dat", average1, average_sqr1, block_number);
    prog_average_std_dev_block_method("data/EX03_1(2).dat", average2, average_sqr2, block_number);

	return 0;
}
