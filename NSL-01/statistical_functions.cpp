/*################################### 
File containing the implementation of 
different statistical function.

LIST of available functions:

- standard deviation;
- block method for uncertainty;

###################################*/

#include <iostream>
#include <fstream>
#include <cmath>
#include "statistical_functions.h"

using namespace std;

double std_dev(double AV, double AV2, int n){
    return sqrt( ( AV2 - pow(AV,2) )/n );
};
