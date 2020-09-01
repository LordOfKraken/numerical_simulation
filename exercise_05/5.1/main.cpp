/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "rng/random.h"

using namespace std;

double error(double, double, int);
Random rng_load();
double r(point);
double metropolis(double, double);

point uniform_step(point, double );
point normal_step(point, double, double);

double psi_1s(point);
double psi_2p(point);

// step_type = 0 --> Uniform distribution, move in a cube
// step_type = 1 --> Normal distribution
// both steps are centered on the previous position

point walk_1s(point ,double ,int step_type = 0);
point walk_2p(point ,double ,int step_type = 0);

// run many times, and return the acceptance rate
double rate_2p(int,double,int);
double rate_1s(int,double,int);

// run the rate function many times with different steps,
// using bisection to find the optimal one.
// the "step" it's the range in which each cartesian cooridinate in uniformly sampled
double find_step_1s(double , double , double , double ,int);
double find_step_2p(double , double , double , double ,int);
// for the normal distribution, the "step" is the standard deviation



//const double a0 = 5.29177E-11;
const double a0 = 1;
Random rnd = rng_load();

int main(int argc, char *argv[])
{
    // Acceptance ratio
    cout << "Finding the optimal step to achieve 50% accept rate"<<endl;
    cout << "Uniform sampling"<< endl;
    double step_1s = find_step_1s(0.5,5,0.5,0.02,0);
    double step_2p = find_step_2p(0.5,5,0.5,0.02,0);
    //double step_1s = 1.20;
    //double step_2p = 3.03;

    cout << "Step 1s = " << step_1s << endl;
    cout << "Step 2p = " << step_2p << endl<<endl;

    cout << "Normal sampling"<<endl;

    double step_1s_g = find_step_1s(0,2,0.5,0.02,1);
    double step_2p_g = find_step_2p(0,2,0.5,0.02,1);
    //double step_1s_g = 0.75;
    //double step_2p_g = 1.875;

    cout << "mu = 0 for both states" << endl;
    cout << "Sigma 1s = " << step_1s_g << endl;
    cout << "Sigma 2p = " << step_2p_g << endl;

    // Equilibration

    cout << "Evolve the two states starting from two different points, and save the distance"<<endl;
    cout << "Analyze the data to find the equilibration steps required"<<endl<<endl;

    // Starting points
    point p_close = {1,1,1};
    point p_far = {100,100,100};

    // First point - Uniform sampling
    point p_1s_close = walk_1s(p_close,step_1s);
    point p_1s_far = walk_1s(p_far,step_1s);
    point p_2p_close = walk_2p(p_close,step_2p);
    point p_2p_far = walk_2p(p_far,step_2p);

    // First point - Normal sampling
    point p_1s_close_g = walk_1s(p_close,step_1s_g,1);
    point p_1s_far_g = walk_1s(p_far,step_1s_g,1);
    point p_2p_close_g = walk_2p(p_close,step_2p_g,1);
    point p_2p_far_g = walk_2p(p_far,step_2p_g,1);

    //point p_1s_close_g = p_close;

    // output files
    ofstream out_close("data/r_close.out");
    ofstream out_far("data/r_far.out");

    if(out_far.is_open() == false)
    {
      cerr << "ERROR: Unable to open output file, maybe data directory is missing?"<<endl;
      exit(1);
    }

    // equilibration steps

    for(int i = 0; i < 1000; i++)
    {
      p_1s_close=walk_1s(p_1s_close,step_1s);
      p_2p_close=walk_2p(p_2p_close,step_2p);
      p_1s_close_g=walk_1s(p_1s_close_g,step_1s_g,1);
      p_2p_close_g=walk_2p(p_2p_close_g,step_2p_g,1);

      out_close << r(p_1s_close) << " " << r(p_1s_close_g) << " " << r(p_2p_close)<< " " << r(p_2p_close_g) <<endl;
      //out_close << r(p_1s_close_g) << endl;
      p_1s_far=walk_1s(p_1s_far,step_1s);
      p_2p_far=walk_2p(p_2p_far,step_2p);
      p_1s_far_g=walk_1s(p_1s_far_g,step_1s_g,1);
      p_2p_far_g=walk_2p(p_2p_far_g,step_2p_g,1);

      out_far << r(p_1s_far) << " " << r(p_1s_far_g) << " " << r(p_2p_far)<< " " << r(p_2p_far_g) <<endl;
    }

    out_far.close();
    out_close.close();

    // the system is equilibrated now

    // Block size
    // study autocorrelation of the data after equilibration.

    ofstream out_close_eq("data/r_close_eq.out");
    ofstream out_far_eq("data/r_far_eq.out");

    for(int i = 0; i < 2000; i++)
    {
      p_1s_close=walk_1s(p_1s_close,step_1s);
      p_2p_close=walk_2p(p_2p_close,step_2p);
      p_1s_close_g=walk_1s(p_1s_close_g,step_1s_g,1);
      p_2p_close_g=walk_2p(p_2p_close_g,step_2p_g,1);

      out_close_eq << r(p_1s_close) << " " << r(p_1s_close_g) << " " << r(p_2p_close)<< " " << r(p_2p_close_g) <<endl;
      //out_close_eq << r(p_1s_close_g) << endl;
      p_1s_far=walk_1s(p_1s_far,step_1s);
      p_2p_far=walk_2p(p_2p_far,step_2p);
      p_1s_far_g=walk_1s(p_1s_far_g,step_1s_g,1);
      p_2p_far_g=walk_2p(p_2p_far_g,step_2p_g,1);

      out_far_eq << r(p_1s_far) << " " << r(p_1s_far_g) << " " << r(p_2p_far)<< " " << r(p_2p_far_g) <<endl;
    }

    out_close_eq.close();
    out_far_eq.close();

    // Data blocking

    int block_size = 50;
    int n_block =10000;
    int iprint = 1000;

    cout << "Starting data blocking with " << n_block << " blocks off size "<< block_size << endl;

    double r_av_1s, r_av_1s_g, r_av_2p, r_av_2p_g;

    double sum_1s = 0;
    double sum_1s_g = 0;
    double sum_2p = 0;
    double sum_2p_g = 0;

    double sum2_1s = 0;
    double sum2_1s_g = 0;
    double sum2_2p = 0;
    double sum2_2p_g = 0;

    ofstream out_r_1s("data/r_av_1s.out");
    ofstream out_r_1s_g("data/r_av_1s_g.out");
    ofstream out_r_2p("data/r_av_2p.out");
    ofstream out_r_2p_g("data/r_av_2p_g.out");

    ofstream out_p_1s("data/points_1s.out");
    ofstream out_p_1s_g("data/points_1s_g.out");
    ofstream out_p_2p("data/points_2p.out");
    ofstream out_p_2p_g("data/points_2p_g.out");

    point p_1s = p_1s_close;
    point p_2p = p_2p_close;
    point p_1s_g = p_1s_close_g;
    point p_2p_g = p_2p_close_g;

    for ( int j = 0; j < n_block; ++j)
    {
      r_av_1s = 0;
      r_av_2p = 0;
      r_av_1s_g = 0;
      r_av_2p_g = 0;

      if ((j+1)%iprint == 0)
      {
        cout << " Current block number = " << j+1<< endl;
      }

      for ( int i = 0; i < block_size; i++)
      {
        p_1s = walk_1s(p_1s, step_1s);
        p_2p = walk_2p(p_2p, step_2p);
        p_1s_g = walk_1s(p_1s_g, step_1s_g, 1);
        p_2p_g = walk_2p(p_2p_g, step_2p_g, 1);

        r_av_1s += r(p_1s);
        r_av_2p += r(p_2p);
        r_av_1s_g += r(p_1s_g);
        r_av_2p_g += r(p_2p_g);
      }

      out_p_1s << p_1s.x << " " << p_1s.y << " " << p_1s.z << endl;
      out_p_2p << p_2p.x << " " << p_2p.y << " " << p_2p.z << endl;
      out_p_1s_g << p_1s_g.x << " " << p_1s_g.y << " " << p_1s_g.z << endl;
      out_p_2p_g << p_2p_g.x << " " << p_2p_g.y << " " << p_2p_g.z << endl;

      r_av_1s /= block_size;
      sum_1s += r_av_1s;
      sum2_1s += pow(r_av_1s, 2);

      r_av_1s_g /= block_size;
      sum_1s_g += r_av_1s_g;
      sum2_1s_g += pow(r_av_1s_g, 2);

      r_av_2p /= block_size;
      sum_2p += r_av_2p;
      sum2_2p += pow(r_av_2p, 2);

      r_av_2p_g /= block_size;
      sum_2p_g += r_av_2p_g;
      sum2_2p_g += pow(r_av_2p_g, 2);

      out_r_1s << sum_1s * 1./(j + 1) << " " << error(sum_1s, sum2_1s, j+1) << endl;
      out_r_1s_g << sum_1s_g * 1./(j + 1) << " " << error(sum_1s_g, sum2_1s_g, j+1) << endl;
      out_r_2p << sum_2p * 1./(j + 1) << " " << error(sum_2p, sum2_2p, j+1) << endl;
      out_r_2p_g << sum_2p_g * 1./(j + 1) << " " << error(sum_2p_g, sum2_2p_g, j+1) << endl;
    }
    return 0;
}


/***********************************
*   Function implementation
***********************************/

double r(point p)
{
  return pow(p.x*p.x + p.y*p.y + p.z*p.z,0.5);
}

point walk_1s(point p_1s,double step,int step_type)
{
  point p_new;
  if(step_type == 1)
    p_new = normal_step(p_1s,0,step);
  else
    p_new = uniform_step(p_1s,step);

  //cout << "punto ( " << p_1s.x << " , " << p_1s.y << " , " << p_1s.z << " )" << endl;
  //cout << " " << psi_1s(p_new*a0) << " " << psi_1s(p_1s*a0) << endl;
  if (rnd.Rannyu() <= metropolis(psi_1s(p_new*a0), psi_1s(p_1s*a0)))
  {
    p_1s = p_new;
  }
  return p_1s;
}

point walk_2p(point p_2p,double step,int step_type)
{
  point p_new;
  if(step_type == 1) //normal distribution
    p_new = normal_step(p_2p,0,step);
  else // uniform distribution
    p_new = uniform_step(p_2p,step);

  //metropolis check
  if (rnd.Rannyu() <= metropolis(psi_2p(p_new), psi_2p(p_2p)))
  {
    p_2p = p_new;
  }
  return p_2p;
}

double rate_1s(int throws,double step,int step_type)
{
  int accept=0;
  point p_1s={0,0,0};
  point p_new;

  for(int i = 0; i < throws; i++)
  {
    if(step_type == 1)
      p_new = normal_step(p_1s,0,step);
    else
      p_new = uniform_step(p_1s,step);
    //cout << psi_1s(p_new) << " " << psi_1s(p_1s) << endl;
    if (rnd.Rannyu() <= metropolis(psi_1s(p_new), psi_1s(p_1s)))
    {
      p_1s = p_new;
      accept++;
    }
  }
  return accept*1./throws;
}

double rate_2p(int throws,double step, int step_type)
{
  int accept=0;
  point p_2p={0,0,0};
  point p_new;

  for(int i = 0; i < throws; i++)
  {
    if(step_type == 1)
      p_new = normal_step(p_2p,0,step);
    else
      p_new = uniform_step(p_2p,step);
    if (rnd.Rannyu() <= metropolis(psi_2p(p_new), psi_2p(p_2p)))
    {
      p_2p = p_new;
      accept++;
    }
  }
  return accept*1./throws;
}

double find_step_1s(double s_min, double s_max, double accept_rate, double tolerance,int step_type)
{
    double test_step = 0.5*(s_max+s_min);
    double throws = 1E6;
    double test_rate = rate_1s(throws,test_step,step_type);

    do{
      test_step = 0.5*(s_max+s_min);
      test_rate = rate_1s(throws,test_step,step_type);
      //cout << test_step << " " << test_rate << endl;
      if(test_rate < accept_rate)
      {
        s_max=test_step;
      }
      else
      {
        s_min = test_step;
      }
    }while(abs(test_rate-accept_rate)>tolerance);
    return test_step;
}

double find_step_2p(double s_min, double s_max, double accept_rate, double tolerance, int step_type)
{
    double test_step = 0.5*(s_max+s_min);
    double throws = 1E6;
    double test_rate = rate_2p(throws,test_step,step_type);

    do{
      test_step = 0.5*(s_max+s_min);
      test_rate = rate_2p(throws,test_step,step_type);

      if(test_rate < accept_rate)
      {
        s_max=test_step;
      }
      else
      {
        s_min = test_step;
      }
    }while(abs(test_rate-accept_rate)>tolerance);
    return test_step;
}

point uniform_step(point p, double step)
{
  return p + rnd.Walk_XYZ(step);
}

point normal_step(point p, double mu, double sigma)
{
  return p + rnd.WalkGauss(mu,sigma);
}

double psi_1s(point p)
{
  double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  //return pow(a0, -3) * M_1_PI * exp(-r *2. / a0);
  return M_1_PI * exp(-r *2. / a0);
}

/*****************************************/

double psi_2p(point p)
{
  double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  double c_theta = p.z / sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  //return pow(a0, -5.) * M_1_PI * pow(r,2) * exp(-r / a0) * pow(c_theta,2)/ 32.;
  return M_1_PI * pow(r,2) * exp(-r / a0) * pow(c_theta,2)/ 32.;
}

/***********************************************/

double metropolis(double p_new, double p_old) {
  return min(1., p_new / p_old);
}

/**********************************************/

double error(double ave, double ave2, int n)
{
  if (n == 1)
    return 0;
  else
    return sqrt((ave2/n-pow(ave/n,2))/(n-1));
}

/************************************************/

Random rng_load()
{
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
   if (input.is_open())
   {
      while ( !input.eof() )
      {
         input >> property;
         if( property == "RANDOMSEED" )
         {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   }
   else cerr << "PROBLEM: Unable to open seed.in" << endl;

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
