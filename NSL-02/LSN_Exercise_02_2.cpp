#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"
#include "random_walk.h"

using namespace std;

int main(int argc, char *argv[]){

    // Initialize the random class rnd
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

    int N = 1E3;                    // Number of steps;
    int L = 1E4;                    // Number of repetition for the trajectory;
    int M = L*N;	        	    // Number of random values needed;
    double *coin = new double[M](); // Define and load a vector with random number distributed uniformly
	for(int i=0; i < M; i++){       // used to decide if take the step forward or backward;
        coin[i] = rnd.Rannyu();
    }
    double prob_backw = 0.5;        // Define the probability of making a step backward;
    int lattice_constant = 1;       // Define lattice constant;
    double *walker_head = new double[3]; // Define the vector containing the xyz coordinate of the head of the walker;
    double* origin = new double[3];      // Define  load and the vector containing
    origin[0] = 0;                       // the xyz coordinate of the origin;
    origin[1] = 0;
    origin[2] = 0;
    random_walk r_walker;
    r_walker.set_step_lenght(lattice_constant);
    r_walker.set_steps_number(N);
    r_walker.set_prob_backw(prob_backw);

    // Estimate of the 3D distance on a cubic lattice done by a RW
    double *random_vec = new double[M]();   // Define and load a vector filled with random
	for(int i=0; i < M; i++){	            // uniformly distributed number used to select
		random_vec[i] = rnd.Rannyu();       // the xyz direction over which take the step;
	}
    double *sum_RW_distance_vec1 = new double[N]();		    // Define a vector containing the sum of all the distances over all trajectories
    double *sum_sqr_RW_distance_vec1 = new double[N]();		// Define a vector containing the sum of all the square distances over all trajectories

	for(int i=0; i < L; i++){      // Cycle over the trajectories
        walker_head[0] = origin[0];        // The walker will always begin in the origin
        walker_head[1] = origin[1];
        walker_head[2] = origin[2];
        r_walker.square_lattice(random_vec, origin, walker_head, i, coin, sum_RW_distance_vec1, sum_sqr_RW_distance_vec1);
    }

    ofstream out_file1;                     // Open a file on which will be printed the average and the
    out_file1.open("data/EX02_2(1).dat");   // standard dev. of the distances traveled by the walker on
    for(int j=0; j < N; j++){               // a square lattice over all the trajectories;
       out_file1 << sum_RW_distance_vec1[j]/L << " " << std_dev(sum_RW_distance_vec1[j]/L,sum_sqr_RW_distance_vec1[j]/L,L) << endl;
    }
    out_file1.close();
/*
    // Estimate of the 3D distance in the continuum done by a RW
    double *random_theta_vec = new double[M]();		// Define random vector to select theta direction
    double theta_min = 0;
    double theta_max = M_PI;
	for(int i=0; i < M; i++){	// Load the vector used to select the direction with random number distributed uniformly over [0,pi]
		random_theta_vec[i] = rnd.Rannyu()*(theta_max-theta_min) + theta_min;
	}
    double *random_phi_vec = new double[M]();		// Define random vector to select phi direction
    double phi_min = 0;
    double phi_max = 2*M_PI;
	for(int i=0; i < M; i++){	// Load the vector used to select the direction with random number distributed uniformly over [0,2pi]
		random_phi_vec[i] = rnd.Rannyu()*(phi_max-phi_min) + phi_min;
	}
    double *sum_RW_distance_vec2 = new double[N]();		    // Define average vector
    double *sum_sqr_RW_distance_vec2 = new double[N]();		// Define average squared vector

    for(int i=0; i < L; i++){                      // Cycle over the trajectories
        walker_head[0] = 0;
        walker_head[1] = 0;
        walker_head[2] = 0;
        for(int j=0; j < N; j++){                  // Cycle over the steps of one trajectory
            int k = j + i*N;
            double RW_distance = 0;

            if(coin[k] <= prob_backw){
                walker_head[0] += -lattice_constant*sin(random_theta_vec[k])*cos(random_phi_vec[k]);
                walker_head[1] += -lattice_constant*sin(random_theta_vec[k])*sin(random_phi_vec[k]);
                walker_head[2] += -lattice_constant*cos(random_theta_vec[k]);
            } else {
                walker_head[0] += lattice_constant*sin(random_theta_vec[k])*cos(random_phi_vec[k]);
                walker_head[1] += lattice_constant*sin(random_theta_vec[k])*sin(random_phi_vec[k]);
                walker_head[2] += lattice_constant*cos(random_theta_vec[k]);
            }

            RW_distance = euclidean_distance(origin, walker_head);
		    sum_RW_distance_vec2[j] += RW_distance;
		    sum_sqr_RW_distance_vec2[j] += pow(RW_distance,2.);
        }
    }

    ofstream out_file2;
    out_file2.open("data/EX02_2(2).dat");
    for(int j=0; j < N; j++){
       out_file2 << sum_RW_distance_vec2[j]/L << " " << std_dev(sum_RW_distance_vec2[j]/L,sum_sqr_RW_distance_vec2[j]/L,L) << endl;
   }
    out_file2.close();
*/
    return 0;
}
