#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "central_limit.h"

using namespace std;

central_limit :: central_limit(){}

central_limit :: ~central_limit(){}

void central_limit :: set_uniform_cases(int x){
	cases = x;
}

double central_limit :: get_uniform_prob(){
	return 1./cases;
}
/*
double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}
*/
