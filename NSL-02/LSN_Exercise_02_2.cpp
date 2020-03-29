#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"

using namespace std;

double euclidean_distance(double* pointA, double* pointB);

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

	int M = 1E6;	        		// Total number of throws
	int N = 1E2;                 	// Number of steps
	int L = static_cast<int>(M/N);  // Number of repetition for the trajectory, please use for M a multiple of N
    double *coin = new double[M]();

	for(int i=0; i < M; i++){   // Load the vector used to select forward or backward step with random number distributed uniformly
        coin[i] = rnd.Rannyu();
    }
    double *origin = new double[3];
    origin[0] = 0;
    origin[1] = 0;
    origin[2] = 0;
    double *walker_head = new double[3];
    int lattice_constant = 1;                      // Define lattice constant
    double prob_backw = 0.5;                       // Define the probability of making a step backward

/*
    // Estimate of the 3D distance on a cubic lattice done by a RW
    double *random_vec = new double[M]();		// Define random vector to select xyz direction
	for(int i=0; i < M; i++){	// Load the vector used to select the direction with random number distributed uniformly
		random_vec[i] = rnd.Rannyu();
	}
    double *sum_RW_distance_vec1 = new double[N]();		    // Define average vector
    double *sum_sqr_RW_distance_vec1 = new double[N]();		// Define average squared vector

	for(int i=0; i < L; i++){                      // Cycle over the trajectories
        walker_head[0] = 0;
        walker_head[1] = 0;
        walker_head[2] = 0;

		for(int j=0; j < N; j++){                  // Cycle over the steps of one trajectory
			int k = j + i*N;
		    double RW_distance = 0;

            if(random_vec[k]> 0 && random_vec[k]<= 1/3.){               // Select the x-direction
                if(coin[k] <= prob_backw){ walker_head[0] += -lattice_constant; }
                else { walker_head[0] += lattice_constant; }
            } else if(random_vec[k]> 1/3. && random_vec[k]<= 2/3.){     // Select the y-direction
                if(coin[k] <= prob_backw){ walker_head[1] += -lattice_constant; }
                else { walker_head[1] += lattice_constant; }
            } else if(random_vec[k]> 2/3. && random_vec[k]<= 1){        // Select the z-direction
                if(coin[k] <= prob_backw){ walker_head[2] += -lattice_constant; }
                else { walker_head[2] += lattice_constant; }
            }

            RW_distance = euclidean_distance(origin, walker_head);
		    sum_RW_distance_vec1[j] += RW_distance;
		    sum_sqr_RW_distance_vec1[j] += pow(RW_distance,2.);
        }

    }

    ofstream out_file1;
    out_file1.open("data/EX02_2(1).dat");
    for(int j=0; j < N; j++){
       out_file1 << sum_RW_distance_vec1[j]/L << " " << std_dev(sum_RW_distance_vec1[j]/L,sum_sqr_RW_distance_vec1[j]/L,L) << endl;
   }
    out_file1.close();
*/
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

    return 0;
}

double euclidean_distance(double* pointA, double* pointB){

    return sqrt( pow(pointA[0]-pointB[0],2.) + pow(pointA[1]-pointB[1],2.) + pow(pointA[2]-pointB[2],2.) );

};
