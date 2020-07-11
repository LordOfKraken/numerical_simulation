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
#include <stdlib.h>     /* atof */
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "rng/random.h"
#include "Monte_Carlo_NVT.h"

using namespace std;

Random rnd = rng_load();

int main(int argc, char *argv[])
{
  // Reading parameters from command line, if present
  if(argc == 4)
  {
    x = atof(argv[1]);
    mu = atof(argv[2]);
    sigma = atof(argv[3]);
  }
  if(argc == 3)
  {
    mu = atof(argv[1]);
    sigma = atof(argv[2]);
  }
  Input(argc); //Inizialization

  //Finding optimal step, ~50% acceptance rate
  // if the optimal step is not reached, use the step provided in the input file
  double opt_step = find_step(0.5,50,0.5,0.02);

  cout << "Metropolis step set at " << opt_step << endl;

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(opt_step);
      if(istep%10 == 0) // Measure every 10 steps
        Measure();
    }
    Averages(iblk);   //Print results for current block
  }
  rnd.SaveSeed();

  return 0;
}
/********************************************************/
/*               Function implementation                */
/********************************************************/


double Psi(double x, double mu, double sigma)
{
  double em = - 0.5 * pow(x-mu,2) /pow(sigma,2);
  double ep = - 0.5 * pow(x+mu,2) /pow(sigma,2);
  return exp(em) + exp(ep);
}

double metropolis(double p_new, double p_old) {
  return min(1., p_new / p_old);
}

double rate(int throws,double step)
{
  double good=0;
  double rate_p_old, rate_p_new;
  double rate_xold, rate_xnew;
  double test_x = x;

  for(int i = 0; i < throws; i++)
  {
    rate_xold = test_x;
    rate_xnew = test_x + step*(rnd.Rannyu() - 0.5) ;
    rate_p_old = pow(Psi(rate_xold, mu, sigma),2);
    rate_p_new = pow(Psi(rate_xnew, mu, sigma),2);

    if(rnd.Rannyu() <= metropolis(rate_p_new,rate_p_old))
    {
          test_x = rate_xnew;
          good++;
    }
  }
  return good*1./throws;
}

double find_step(double s_min, double s_max, double accept_rate, double tolerance)
{
    double test_step = 0.5*(s_max+s_min);
    double throws = 1E4;
    double test_rate = rate(throws,test_step);
    int max_retry = 100, counter=0;

    do{
      test_step = 0.5*(s_max+s_min);
      test_rate = rate(throws,test_step);
      if(test_rate < accept_rate)
      {
        s_max = test_step;
      }
      else
      {
        s_min = test_step;
      }
      //cout << "test step = " << test_step << " --> " << test_rate<<endl;
      counter++;
      if(counter == max_retry)
      {
        cout << "Max retry to find optimal step reached" <<endl;
        cout << "Step = " << test_step <<endl << "Acceptance rate = " << test_rate <<endl;
        cout << "Using the step provided in the input file" << endl;
        test_step = delta;
      }
    }while(counter < max_retry && abs(test_rate-accept_rate)>tolerance);
    return test_step;
}

double V_ext(double x)
{
  return pow(x,4) - 2.5*pow(x,2);
}

double H(double x, double mu, double sigma)
{
    double p2_m = exp(-0.5*pow(x-mu,2)/pow(sigma,2)) * (mu*mu - sigma*sigma + x*x - 2.*mu*x)/pow(sigma,4);
    double p2_p = exp(-0.5*pow(x+mu,2)/pow(sigma,2)) * (mu*mu - sigma*sigma + x*x + 2.*mu*x)/pow(sigma,4);
    return -0.5 * (p2_m + p2_p) + V_ext(x);
}

double Error(double sum, double sum2, int iblk) {
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Input(int args){

    ifstream ReadInput;

    cout << "1D Monte Carlo simulation of a single particle" << endl;
    cout << "External potential V(x) = x^4 + 5/2 * x^2  " << endl << endl;

    ReadInput.open("input.dat");

    ReadInput >> step;
    ReadInput >> nblk;
    ReadInput >> nstep;
    if(args == 4)
    {
      cout <<"Reading x,mu,sigma from command line" <<endl;
    }
    else if(args == 3)
    {
      ReadInput >> x;
      cout <<"Reading mu,sigma from command line" <<endl;
    }
    else
    {
      ReadInput >> x;
      ReadInput >> mu;
      ReadInput >> sigma;
    }
    //hist=new Histogram(-2.5, 2.5, 100, true);

    cout << "This program perform Metropolis moves, with a 50\% acceptance rate" << endl;
    cout << "Default step length = " << step << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Steps per block = " << nstep << endl;
    cout << "Starting point = " << x << endl;
    cout << "mu = " << mu << endl;
    cout << "sigma = " << sigma << endl << endl;

    ih = 0;
    n_props = 1;

    nbins = 100;
    x_min = -3;
    x_max = 3;
    //n_props = n_props + nbins;
    hist_norm = 0;
    bin_size = (x_max-x_min)*1./nbins;

    ReadInput.close();
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

void Move(double step)
{
  double p_old, p_new;
  double xold, xnew;

  xold = x;
  xnew = x + step*(rnd.Rannyu() - 0.5) ;
  p_old = pow(Psi(xold, mu, sigma),2);
  p_new = pow(Psi(xnew, mu, sigma),2);

  if(rnd.Rannyu() <= metropolis(p_new,p_old))
  {
        x = xnew;
        accepted++;
  }
  attempted++;

  // ofstream prova;
  // prova.open("test.dat",ios::app);
  // prova << x << endl;

}

void Measure()
{
  stima_h = H(x,mu,sigma); //Potential energy per particle
  blk_av[ih] += stima_h;
  blk_norm = blk_norm + 1.0;

  //update the histogram
  if(x >= x_min && x <= x_max)
  {
    int bin = (int)((x-x_min)/bin_size);
    hist[bin] += 1.0;
    hist_norm += 1.0;
  }
}

void Averages(int iblk) //Print results for current block
{
    ofstream AveH, Variation;
    ofstream Histo;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << (double)accepted/attempted << endl;

    glob_av[ih]+=blk_av[ih]/blk_norm;
    glob_av2[ih]+=pow(blk_av[ih]/blk_norm,2);

    //string filepath = "data/ave_h";
    //string parameters = to_string(mu) + "_" + to_string(sigma);
    //string filename = filepath + "_" + parameters + ".out";

    AveH.open("data/ave_h.out",ios::app); // print multiple files, with the blocking parameters
    Variation.open("minimize.out",ios::app);
    Histo.open("data/histo.out",ios::app);

    AveH << iblk << " " << glob_av[ih]/(double)iblk << " " << Error(glob_av[ih], glob_av2[ih], iblk) << endl;

    if(iblk==nblk) //save the last averages in the same file, with the parameters
    {
      Variation << mu << " " << sigma << " " << glob_av[ih]/(double)iblk << " " << Error(glob_av[ih], glob_av2[ih], iblk) << endl;
      for(int i = 0; i < nbins; i++)
      {
          double x_bin = x_min + bin_size * (0.5 + i) ;
          Histo << x_bin << " " << hist[i]*1./(hist_norm*bin_size) << endl;
      }
    }
    AveH.close();
    Variation.close();
    Histo.close();
}

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
