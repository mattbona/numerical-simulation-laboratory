#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "asset.h"
 
using namespace std;

Asset :: Asset(){asset_value = 0;}
 
Asset :: ~Asset(){}

// Class Asset Methods
double Asset :: GetAssetValue(){
	return asset_value;
};

void Asset :: SetAssetValue(double my_asset_value){
	asset_value = my_asset_value;
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

void GBM :: UpdateAssetValue(double t){
	double old_asset_value = this->GetAssetValue();
	double new_asset_value = old_asset_value*exp( (mu - 0.5*pow(sigma,2.) )*t + sigma*gaussian_var  );
	this->SetAssetValue(new_asset_value);
};
