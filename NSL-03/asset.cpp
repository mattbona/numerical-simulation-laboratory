#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "asset.h"

using namespace std;

Asset :: Asset(){asset_price = 0;}

Asset :: ~Asset(){}

// Class Asset Methods
double Asset :: GetAssetPrice(){
	return asset_price;
};

void Asset :: SetAssetPrice(double my_asset_price){
	asset_price = my_asset_price;
};

GBM :: GBM(double my_mu, double my_sigma){
	mu = my_mu;
	sigma = my_sigma;
};

GBM :: ~GBM(){}

// Class GBM Methods
void GBM :: SetGaussianVar(double my_gaussian_var){
	gaussian_var = my_gaussian_var;
};

double GBM :: GetMu(){
	return mu;
};

double GBM :: GetSigma(){
	return sigma;
};

void GBM :: UpdateAssetPrice(double t){
	double old_asset_price = this->GetAssetPrice();
	double new_asset_price = old_asset_price*exp( (mu - 0.5*pow(sigma,2.) )*t + sigma*gaussian_var*sqrt(t) );
	this->SetAssetPrice(new_asset_price);
};
