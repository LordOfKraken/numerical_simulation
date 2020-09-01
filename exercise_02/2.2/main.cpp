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
#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include "rng/random.h"

using namespace std;

Random rng_load();
double Error(double, double, int);

vector <double> s_lat, s2_lat;
vector <double> s_con, s2_con;

int main (int argc, char *argv[])
{
  Random rnd = rng_load();

  double l_lat=1;
  double l_con=1;
  int n_walk = 10000;
  int n_step = 100;

  int dir;
  double theta, phi;
  double x_lat, y_lat, z_lat;
  double x_con, y_con, z_con;
  double r2_lat, r2_con;

  // contains the sum and sum2 of the distance at each step.
  s_lat.resize(n_step);
  s2_lat.resize(n_step);
  s_con.resize(n_step);
  s2_con.resize(n_step);

  for(int j=0; j< n_walk; j++)
  {
    // every walk starts from 0
    x_lat=0, y_lat=0, z_lat=0;
    x_con=0, y_con=0, z_con=0;

    r2_lat=0, r2_con=0;
    for(int i=0; i < n_step; i++)
    {
      //lattice walk section
      dir=int(rnd.Rannyu(0,6));
      switch(dir)
      {
        case 0:
        x_lat++;
        break;
        case 1:
        y_lat++;
        break;
        case 2:
        z_lat++;
        break;
        case 3:
        x_lat--;
        break;
        case 4:
        y_lat--;
        break;
        case 5:
        z_lat--;
        break;
      }
      r2_lat = pow(x_lat*l_lat,2) + pow(y_lat*l_lat,2) + pow(z_lat*l_lat,2);
      s_lat[i] += sqrt(r2_lat);
      if(j==0)
        cout << i << " " << sqrt(r2_lat) << " " << s_lat[i] << endl;
      s2_lat[i] += r2_lat;
      if(j==0)
        cout << i << " " << r2_lat << " " << s2_lat[i] << endl;

      //continuum walk section
      phi = rnd.Rannyu(0,2.*M_PI);
      theta = acos(1 - 2.*rnd.Rannyu());

      x_con += l_con * cos(phi) * sin(theta);
      y_con += l_con * sin(phi) * sin(theta);
      z_con += l_con * cos(theta);

      r2_con = pow(x_con,2) + pow(y_con,2) + pow(z_con,2);
      s_con[i] += sqrt(r2_con);
      s2_con[i] += r2_con;
    }
  }

/********************************
* Saving to file the last simulation
********************************/

  ofstream l_out("data/lattice.out");
  if (l_out.is_open())
  {
    for (int i = 0; i < n_step; i++)
      l_out << i + 1 << " " << sqrt(s2_lat[i]/n_walk) << " " << Error(s_lat[i], s2_lat[i] , n_walk) << endl;
  }
  else
    cerr <<"Unable to open lattice output file: data saving failed" <<endl;
  l_out.close();

  ofstream c_out("data/continuum.out");
  if (c_out.is_open())
  {
    for (int i = 0; i < n_step; i++)
      c_out << i + 1 << " " << sqrt(s2_con[i]/n_walk) << " " << Error(s_con[i] , s2_con[i], n_walk) << endl;
  }
  else
    cerr <<"Unable to open continuum output file: data saving failed" <<endl;
  c_out.close();
  return 0;
}

/***********************************
*   Function implementation
***********************************/

double Error(double sum, double sum2, int iblk)
{
  return sqrt((sum2/iblk - pow(sum/iblk,2))/iblk);
}

/*********************************/

Random rng_load() {

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("rng/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("rng/seed.in");
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

   rnd.SaveSeed();
   return rnd;
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
