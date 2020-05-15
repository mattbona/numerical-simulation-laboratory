/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include "statistical_functions.h"

using namespace std;

int main(){
  Input();             //Inizialization

  cout << "Start the equilibration..." << endl << endl;
  for(int i=0; i<eqsteps; ++i)
      Move();

//  int nconf = 1;
  int steps_per_block = int(nstep/nblocks);
  cout << "Start the simulation..." << endl;
  for(int iblock=0; iblock <= nblocks ; ++iblock){
      if(iblock%10 == 0) cout << "Number of block: " << iblock << endl;
      for(int istep=1; istep <=steps_per_block; ++istep){
         Move();           //Move particles with Verlet algorithm
         if(istep%10 == 0){
            Measure(iblock);     //Properties measurement
//            ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
//            nconf += 1;
         }
      }
      av_Epot[iblock] /= steps_per_block/10;
      av_EKin[iblock] /= steps_per_block/10;
      av_Etot[iblock] /= steps_per_block/10;
      av_Temp[iblock] /= steps_per_block/10;
      av2_Epot[iblock] += av_Epot[iblock]*av_Epot[iblock];
      av2_EKin[iblock] += av_EKin[iblock]*av_EKin[iblock];
      av2_Etot[iblock] += av_Etot[iblock]*av_Etot[iblock];
      av2_Temp[iblock] += av_Temp[iblock]*av_Temp[iblock];

  }

  prog_average_std_dev_block_method("results/e_pot.dat", av_Epot, av2_Epot, nblocks);
  prog_average_std_dev_block_method("results/e_kin.dat", av_EKin, av2_EKin, nblocks);
  prog_average_std_dev_block_method("results/e_tot.dat", av_Etot, av2_Etot, nblocks);
  prog_average_std_dev_block_method("results/temp.dat", av_Temp, av2_Temp, nblocks);

  ConfFinal();         //Write final configuration to restart

  return 0;
}

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator

  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> eqsteps;
  ReadInput >> restart;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Number of equilibration steps = " << eqsteps << endl << endl;
  if(restart) cout << "Velocities'll be rescaled to match the desired temperature using the configuration file: old.0." << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  //Read old configuration
  bool old_conf_file_is_open = 0;
  if(restart){
      ifstream ReadOldConf("old.0");
      old_conf_file_is_open = ReadOldConf.is_open();
      if(old_conf_file_is_open) {
          cout << "Read initial configuration from file old.0 " << endl;
          for (int i=0; i<npart; ++i){
              ReadOldConf >> xold[i] >> yold[i] >> zold[i];
              xold[i] = xold[i] * box;
              yold[i] = yold[i] * box;
              zold[i] = zold[i] * box;
          }
          ReadOldConf.close();
      } else cout << "No old.0 file found. Run the simulation with restart=0 at least one time to get it!" << endl;
  }

    //Prepare initial velocities
    if(restart && old_conf_file_is_open){

        Move();

        double t = 0.0;
        for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
        stima_temp = (2.0 / 3.0) * t/(double)npart;

        cout << "Preparing velocities according to the old configuration (mesured T = "<< stima_temp << ") rescaling them in agreement with desired temperature T = " << temp << endl << endl;

        double fs=sqrt(temp/stima_temp);

        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;

            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
        return;

    }else {
       cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
       double sumv[3] = {0.0, 0.0, 0.0};
       for (int i=0; i<npart; ++i){
         vx[i] = rand()/double(RAND_MAX) - 0.5;
         vy[i] = rand()/double(RAND_MAX) - 0.5;
         vz[i] = rand()/double(RAND_MAX) - 0.5;

         sumv[0] += vx[i];
         sumv[1] += vy[i];
         sumv[2] += vz[i];
       }
       for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
       double sumv2 = 0.0, fs;
       for (int i=0; i<npart; ++i){
         vx[i] = vx[i] - sumv[0];
         vy[i] = vy[i] - sumv[1];
         vz[i] = vz[i] - sumv[2];

         sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
       }
       sumv2 /= (double)npart;

       fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
       for (int i=0; i<npart; ++i){
         vx[i] *= fs;
         vy[i] *= fs;
         vz[i] *= fs;

         xold[i] = Pbc(x[i] - vx[i] * delta);
         yold[i] = Pbc(y[i] - vy[i] * delta);
         zold[i] = Pbc(z[i] - vz[i] * delta);
       }
    }
   return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure(int iblock){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    av_Epot[iblock] += stima_pot;
    av_EKin[iblock] += stima_kin;
    av_Etot[iblock] += stima_etot;
    av_Temp[iblock] += stima_temp;

    return;
}


void ConfFinal(void){ //Write final configuration
    ofstream WriteConf;
    ofstream WriteOldConf;

    cout << "Print final configuration to file config.final and old.final" << endl << endl;
    WriteConf.open("config.final");
    WriteOldConf.open("old.final");

    for (int i=0; i<npart; ++i){
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    }
    for (int i=0; i<npart; ++i){
        WriteOldConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    }

    WriteOldConf.close();
    WriteConf.close();
    return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
