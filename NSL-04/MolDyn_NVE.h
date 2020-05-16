/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp,stima_pres;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
int nblocks = 100;
double delta;

//Prepare array for blocking method
double *av_Epot  = new double[nblocks]();
double *av_EKin  = new double[nblocks]();
double *av_Etot  = new double[nblocks]();
double *av_Temp  = new double[nblocks]();
double *av_Pres  = new double[nblocks]();

double *av2_Epot  = new double[nblocks]();
double *av2_EKin  = new double[nblocks]();
double *av2_Etot  = new double[nblocks]();
double *av2_Temp  = new double[nblocks]();
double *av2_Pres  = new double[nblocks]();

//equilibration
bool restart;
int eqsteps;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(int);
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
