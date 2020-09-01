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
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exp(double mean) {
  double n = Rannyu();
  return -1.*log(1. - n)/mean;
}

double Random :: Lorentz(double mean, double gamma) {
  double n=Rannyu();
  return mean + gamma * tan(M_PI*(n-0.5));
}

double Random :: Importance()
{
  // Sampling from the inverse CDF of a linear PDF 2*(1-x)
  return 1 + sqrt(1-Rannyu());
}

point Random :: Walk(double step)
{
  double theta, phi;
  theta =2*M_PI*Rannyu();
  phi = acos( 1. - 2. * Rannyu());
  point p(step * sin(theta) * cos(phi), step * sin(theta) * sin(phi), step * cos(phi));
  return p;
}

point Random :: Walk(double step, point p)
{
  double theta, phi;
  theta =2*M_PI*Rannyu();
  phi = acos( 1. - 2. * Rannyu());
  point a(step * sin(theta) * cos(phi), step * sin(theta) * sin(phi), step * cos(phi));
  return p+a;
}

point Random :: Walk_XYZ(double step)
{
  point p(Rannyu(-step,step),Rannyu(-step,step),Rannyu(-step,step));
  return p;
}

point Random :: Walk_XYZ(double step, point p)
{
  point p_new(Rannyu(-step,step),Rannyu(-step,step),Rannyu(-step,step));
  return p+p_new;
}

point Random :: WalkGauss(double mean, double sigma) {
   double x=sqrt(-2.*log(1.-Rannyu()))*cos(2.*M_PI*Rannyu());
   double y=sqrt(-2.*log(1.-Rannyu()))*cos(2.*M_PI*Rannyu());
   double z=sqrt(-2.*log(1.-Rannyu()))*cos(2.*M_PI*Rannyu());
   x = mean + x * sigma;
   y = mean + y * sigma;
   z = mean + z * sigma;
   point p(x,y,z);
   return p;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
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
