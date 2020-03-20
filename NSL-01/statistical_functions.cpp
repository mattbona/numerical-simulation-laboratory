/*################################### 
File containing the implementation of 
different statistical function.

LIST of available functions:

- standard deviation: std_dev();
- block method for uncertainty: prog_average_std_dev_block_method();

###################################*/

#include "statistical_functions.h"

double std_dev(double average, double sqrd_average, int sample_number){
	/*
	Returns the standard deviation normalized 
	over the square root of the number of samples. 
	*/

    return sqrt( ( sqrd_average - pow(average,2) )/sample_number );
};

void prog_average_std_dev_block_method(const string& output_file, double* average_array, double* sqrd_average_array, int number_blocks){
	/*
	Prints a file in data/ containing the progressive 
	average and std-dev of the measures over the 
	Monte Carlo blocks.
	*/

        ofstream out_file;
        out_file.open(output_file);
        double *prog_average = new double[number_blocks]();         // Define progressive average vector
        double *prog_average_sqr = new double[number_blocks]();     // Define progressive average squared vector
        double *prog_error = new double[number_blocks]();           // Define progressive error vector

	for(int i=0; i < number_blocks; i++){ 

		for(int j=0; j < i+1; j++){
        		prog_average[i] += average_array[j];
        		prog_average_sqr[i] += sqrd_average_array[j];
		}

		prog_average[i] = prog_average[i]/(i+1);
		prog_average_sqr[i] = prog_average_sqr[i]/(i+1);
		prog_error[i] = std_dev(prog_average[i],prog_average_sqr[i],i);

		out_file << prog_average[i] << " " << prog_error[i] << endl;
	}
	
	out_file.close();
};	
