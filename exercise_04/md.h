/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <vector>
#include <tuple>

using namespace std;
//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie,ip;
double stima_epot, stima_ekin, stima_etot, stima_temp, stima_pres;
int igofr;

double blk_av[m_props];
double glob_av[m_props],glob_av2[m_props];
double blk_norm, bin_size;
int nbins;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;
int old;


// simulation
int nstep, nblock, iprint, seed, eqstep;
double delta;

//blocking_method
//vector<double> v_epot, v_ekin, v_etot, v_temp, v_pres;

//functions
void Input(string,string);
void Move(void);
void ConfFinal(string);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);
void Averages(int);
void Reset(int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
