#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random_walk.h"
#include "random.h"

using namespace std;

random_walk :: random_walk(){}

random_walk :: ~random_walk(){}

void random_walk :: set_step_lenght(int my_step_lenght){
    step_lenght = my_step_lenght;
};

void random_walk :: set_steps_number(int my_steps_number){
    steps_number = my_steps_number;
};

void random_walk :: set_prob_backw(double my_prob_backw){
    prob_backw = my_prob_backw;
};

double random_walk :: euclidean_distance(double* pointA, double* pointB){
    return sqrt( pow(pointA[0]-pointB[0],2.) + pow(pointA[1]-pointB[1],2.) + pow(pointA[2]-pointB[2],2.) );
};

void random_walk :: square_lattice(double* random_vec, double* origin, double* walker_head, int repetition, double* coin_toss, double *sum_RW_distance_vec, double *sum_sqr_RW_distance_vec){

    for(int j=0; j < steps_number; j++){                  // Cycle over the steps of one trajectory
        int k = j + repetition*steps_number;
        double RW_distance = 0;

        if(random_vec[k]> 0 && random_vec[k]<= 1/3.){               // Select the x-direction
            if(coin_toss[k] <= prob_backw){ walker_head[0] += -step_lenght; }
            else { walker_head[0] += step_lenght; }
        } else if(random_vec[k]> 1/3. && random_vec[k]<= 2/3.){     // Select the y-direction
            if(coin_toss[k] <= prob_backw){ walker_head[1] += -step_lenght; }
            else { walker_head[1] += step_lenght; }
        } else if(random_vec[k]> 2/3. && random_vec[k]<= 1){        // Select the z-direction
            if(coin_toss[k] <= prob_backw){ walker_head[2] += -step_lenght; }
            else { walker_head[2] += step_lenght; }
        }

        RW_distance = euclidean_distance(origin, walker_head);
        sum_RW_distance_vec[j] += RW_distance;
        sum_sqr_RW_distance_vec[j] += pow(RW_distance,2.);
    }
};
