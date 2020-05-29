/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "LSN_Exercise_081.h"

using namespace std;

int main()
{
  // Inizialization
  Input();
  // Equilibration
  cout << "Equilibration..." << endl << endl;
  for(int ieq=0; ieq<neqstep; ieq++)
  {
      Move();
  }

  // Simulation
  for(int iblk=1; iblk <= nblk; ++iblk)
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Single 1D quantum particle         " << endl;
  cout << "Variational Monte Carlo simulation." << endl << endl;
  cout << "External confinement potential: v(x) = x^4 - 5/2 * x^2" << endl;
  cout << "Approximate ground state function:" << endl;
  cout << "\texp( -(x-mu)**2 / 2*sigma**2 ) + exp( -(x+mu)**2 / 2*sigma**2 )" << endl << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> mu;
  ReadInput >> sigma;
  ReadInput >> delta;
  ReadInput >> neqstep;
  ReadInput >> nblk;
  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of equilibration steps = " << neqstep << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  ih = 0; // Hamiltonian
  n_props = 1; //Number of observables

}

void Move(void)
{

}

void Measure()
{

}

void Reset(int iblk) //Reset block averages
{

   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{

    ofstream H_ave;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    H_ave.open("results/h_ave.dat",ios::app);

    stima_ham = 0; //Potential energy
    glob_av[ih] += stima_ham;
    glob_av2[ih] += stima_ham*stima_ham;
    err_ham=Error(glob_av[ih],glob_av2[ih],iblk);

//Pressure
    H_ave << " "<< iblk << " " << glob_av[ih]/(double)iblk << " "<< err_ham << endl;

    cout << "----------------------------" << endl << endl;

    H_ave.close();

}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
