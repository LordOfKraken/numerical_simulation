/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __NVT__
#define __NVT__

//parameters, observables
const int m_props=1000;
int n_props, ih;
double bin_size,nbins,x_min,x_max,hist_norm;
double hist[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_h,err_h;

//configuration
double x,mu,sigma;

// simulation
int nstep, nblk;
double step,delta;

//functions
void Input(int);
void Reset(int);
void Averages(int);
void Move(double);
void ConfFinal(void);
void Measure(void);
double Error(double,double,int);
double metropolis(double, double);
double rate(int,double);
double find_step(double , double , double , double);
Random rng_load();


#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
