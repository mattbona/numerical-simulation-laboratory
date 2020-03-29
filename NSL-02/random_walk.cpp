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
    /*
    Method to set the step lenght of the random walk.
    */
    step_lenght = my_step_lenght;
};

void random_walk :: set_steps_number(int my_steps_number){
    /*
    Method to set the number of steps of the random walk.
    */
    steps_number = my_steps_number;
};

void random_walk :: set_prob_backw(double my_prob_backw){
    /*
    Method to set the probability of the walker to do a step backward.
    */
    prob_backw = my_prob_backw;
};

double random_walk :: euclidean_distance(double* pointA, double* pointB){
    /*
    Method that return the euclidean distance provided the x,y,z coordinate
    of two given point in space.
    */
    return sqrt( pow(pointA[0]-pointB[0],2.) + pow(pointA[1]-pointB[1],2.) + pow(pointA[2]-pointB[2],2.) );
};

void random_walk :: square_lattice(double* random_vec, double* origin, double* walker_head, int repetition, double* coin_toss, double *sum_RW_distance_vec, double *sum_sqr_RW_distance_vec){
    /*
    Method that will update two given vectors with the sum over all repetition of the walk of the
    distance traveled and the sum over all repetition of the walk of the square distance traveled on a square lattice.

    This method needs:
    - double *random_vec: a pointer to a vector of double containing $steps_number uniformly distributed random number
                          that will be used to choose in which direction take the step;
    - double *origin: a pointer to a vector of double containing the three coordinate of the origin;
    - double *walker_head: a pointer to a vector of double containing th three coordinate of the head of the walker;
    - int repetition: a int with the number of time that the walk will be repeated;
    - double *coin_toss: a pointer to a vector of double containing $steps_number uniformly distributed random number
                         that will be used to decide if take the step backward;
    - double *sum_RW_distance_vec: a pointer to the vector of double containing
                                   the $steps_number sum of distances over all trajectories;
    - double *sum_sqr_RW_distance_vec: a pointer to the vector of double containing
                                       the $steps_number sum of square distances over all trajectories;
    */
    for(int j=0; j < steps_number; j++){                  // Cycle over the steps of one trajectory
        int k = j + repetition*steps_number;
        double RW_distance = 0;

        if(random_vec[k]> 0 && random_vec[k]<= 1/3.){               // Select the x-direction
            if(coin_toss[k] <= prob_backw){ walker_head[0] += -step_lenght; }   // Decide if take the step backward
            else { walker_head[0] += step_lenght; }
        } else if(random_vec[k]> 1/3. && random_vec[k]<= 2/3.){     // Select the y-direction
            if(coin_toss[k] <= prob_backw){ walker_head[1] += -step_lenght; }   // Decide if take the step backward
            else { walker_head[1] += step_lenght; }
        } else if(random_vec[k]> 2/3. && random_vec[k]<= 1){        // Select the z-direction
            if(coin_toss[k] <= prob_backw){ walker_head[2] += -step_lenght; }   // Decide if take the step backward
            else { walker_head[2] += step_lenght; }
        }

        RW_distance = euclidean_distance(origin, walker_head);
        sum_RW_distance_vec[j] += RW_distance;
        sum_sqr_RW_distance_vec[j] += pow(RW_distance,2.);
    }
};

void random_walk :: continuum(double* random_theta_vec, double* random_phi_vec, double* origin, double* walker_head, int repetition, double* coin_toss, double *sum_RW_distance_vec, double *sum_sqr_RW_distance_vec){
    /*
    Method that will update two given vectors with the sum over all repetition of the walk of the
    distance traveled and the sum over all repetition of the walk of the square distance traveled
    in the continuum.

    This method needs:
    - double* random_theta_vec: a pointer to a vector of double containing $steps_number uniformly distributed random number
                                over [0, pi] that will be used to choose in which theta-direction take the step;
    - double* random_phi_vec: a pointer to a vector of double containing $steps_number uniformly distributed random number
                                over [0, 2pi] that will be used to choose in which phi-direction take the step;
    - double *origin: a pointer to a vector of double containing the three coordinate of the origin;
    - double *walker_head: a pointer to a vector of double containing th three coordinate of the head of the walker;
    - int repetition: a int with the number of time that the walk will be repeated;
    - double *coin_toss: a pointer to a vector of double containing $steps_number uniformly distributed random number
                         that will be used to decide if take the step backward;
    - double *sum_RW_distance_vec: a pointer to the vector of double containing
                                   the $steps_number sum of distances over all trajectories;
    - double *sum_sqr_RW_distance_vec: a pointer to the vector of double containing
                                       the $steps_number sum of square distances over all trajectories;
    */

    for(int j=0; j < steps_number; j++){          // Cycle over the steps of one trajectory
        int k = j + repetition*steps_number;
        double RW_distance = 0;

        if(coin_toss[k] <= prob_backw){           // Decide if take the step backward
            walker_head[0] += -step_lenght*sin(random_theta_vec[k])*cos(random_phi_vec[k]); // Update x-position
            walker_head[1] += -step_lenght*sin(random_theta_vec[k])*sin(random_phi_vec[k]); // Update y-position
            walker_head[2] += -step_lenght*cos(random_theta_vec[k]);                        // Update z-position
        } else {
            walker_head[0] += step_lenght*sin(random_theta_vec[k])*cos(random_phi_vec[k]);
            walker_head[1] += step_lenght*sin(random_theta_vec[k])*sin(random_phi_vec[k]);
            walker_head[2] += step_lenght*cos(random_theta_vec[k]);
        }

        RW_distance = euclidean_distance(origin, walker_head);
        sum_RW_distance_vec[j] += RW_distance;
        sum_sqr_RW_distance_vec[j] += pow(RW_distance,2.);
    }
};
