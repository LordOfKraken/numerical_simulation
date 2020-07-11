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
#include "md.h"

using namespace std;

int main(int argc, char *argv[])
{
  string input_file="input.dat",config_path="config/";

  if(argc > 1 )
  {
    input_file = argv[1];
    if( argc == 3)
    {
      config_path = argv[2];
      if(config_path.length() > 0  && config_path.at(config_path.length()-1) != '/')
      {
        config_path.append("/");
      }
    }
  }

  cout << "Config path = " << config_path << endl;
  cout << "Input file = " << input_file << endl;

  Input(input_file,config_path);             //Inizialization
  int nconf = 1;

  cout << "Equilibrating for " << eqstep << " steps" << endl;
  for(int i=0; i < eqstep; i++)
      Move();

  cout << "Starting the simulation" << endl;
  for(int iblk=1; iblk<=nblock; ++iblk){
      cout << "Block number : " << iblk << endl;
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep){
          Move();           //Move particles with Verlet algorithm
          if(istep%10 == 0){
              Measure();     //Properties measurement and update block averages
              //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
              nconf += 1;
          }
      }
      Averages(iblk);   //Print results for current block
  }
  ConfFinal(config_path);         //Write final configuration to restart
  return 0;
}


void Input(string input_file, string config_path){ //Prepare all stuff for the simulation
    ifstream ReadInput,ReadConf,ReadOldConf;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    seed = 1;    //Set seed for random numbers
    srand(seed); //Initialize random number generator

    ReadInput.open(input_file); //Read input

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
    ReadInput >> nblock;
    ReadInput >> eqstep;
    ReadInput >> iprint;
    ReadInput >> old;
    ReadInput.close();

    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps per block= " << nstep << endl;
    cout << "Number of blocks = " << nblock << endl << endl;
    cout << "The system configuration files are saved in the corresponding \'config\' folder"<<endl;
    cout << "The output data is saved in the \'data\' folder" << endl;
    cout << "Make sure these two folders exists, for the program to be able to run properly "<<endl<<endl;

    //Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    ip = 4; //Pressure
    n_props = 5; //Number of observables

    //measurement of g(r)
    igofr = 5;
    nbins = 100;
    n_props = n_props + nbins;
    bin_size = (box/(2.0*(double)nbins));

    //Read initial configuration
    string start_conf = config_path;
    string old_conf = config_path;

    start_conf.append("config.0");
    old_conf.append("old.0");
    cout << "Read initial configuration from file "<< start_conf << endl << endl;

    ReadConf.open(start_conf);
    for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
    }
    ReadConf.close();

    //Read old configuration
    if(old == 1)
    {
      cout << "Read old configuration from file " << old_conf << endl;
      ReadOldConf.open(old_conf);
      for (int i=0; i<npart; ++i){
          ReadOldConf >> xold[i] >> yold[i] >> zold[i];
          xold[i] = xold[i] * box;
          yold[i] = yold[i] * box;
          zold[i] = zold[i] * box;
      }
      ReadOldConf.close();

      cout << "Create velocities with Verlet" << endl << endl;

      Move();
    }
    else
    {
      cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
      for (int i=0; i<npart; ++i)
      {
          vx[i] = rand()/double(RAND_MAX) - 0.5;
          vy[i] = rand()/double(RAND_MAX) - 0.5;
          vz[i] = rand()/double(RAND_MAX) - 0.5;
      }
    }

    cout << "Rescale velocities according to desired temperature" << endl;
    double sumv[3] = {0.0, 0.0, 0.0};

    for (int i=0; i<npart; ++i)
    {
        sumv[0] += vx[i];
        sumv[1] += vy[i];
        sumv[2] += vz[i];
    }

      for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
      double sumv2 = 0.0, fs;
      for (int i=0; i<npart; ++i)
      {
          vx[i] = vx[i] - sumv[0];
          vy[i] = vy[i] - sumv[1];
          vz[i] = vz[i] - sumv[2];
          sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
      }
      sumv2 /= (double)npart;

      fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
      for (int i=0; i<npart; ++i)
      {
          vx[i] *= fs;
          vy[i] *= fs;
          vz[i] *= fs;

          xold[i] = Pbc(x[i] - vx[i] * delta);
          yold[i] = Pbc(y[i] - vy[i] * delta);
          zold[i] = Pbc(z[i] - vz[i] * delta);
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

void Measure(){ //Properties measurement
    int bin;
    double v, p, t, vij, pij;
    double dx, dy, dz, dr;
    ofstream Epot, Ekin, Etot, Temp, Pres;

    Epot.open("data/epot.out",ios::app);
    Ekin.open("data/ekin.out",ios::app);
    Etot.open("data/etot.out",ios::app);
    Temp.open("data/temp.out",ios::app);
    Pres.open("data/pres.out",ios::app);

    v = 0.0; //reset observables
    p = 0.0;
    t = 0.0;

    //cycle over pairs of particles
    for (int i=0; i<npart-1; ++i){
        for (int j=i+1; j<npart; ++j){

            dx = Pbc( x[i] - x[j] );
            dy = Pbc( y[i] - y[j] );
            dz = Pbc( z[i] - z[j] );
            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);

            //update of the histogram of g(r)
            bin = igofr + (int)(dr/bin_size);
            if(bin < igofr + nbins)
                blk_av[bin] += 1.0;

            if(dr < rcut){
                vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
                v += vij;

                //Potential energy and virial

                pij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
                p += pij;
            }
        }
    }
    //Kinetic energy
    for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_epot = v/(double)npart; //Potential energy per particle
    stima_ekin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_pres = rho * stima_temp + 48.0 * p / (3.0 * vol); //Pressure

    blk_av[iv]+=stima_epot;
    blk_av[ik]+=stima_ekin;
    blk_av[ie]+=stima_etot;
    blk_av[it]+=stima_temp;
    blk_av[ip]+=stima_pres;

    blk_norm = blk_norm + 1.0;

    Epot << stima_epot  << endl;
    Ekin << stima_ekin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();

    return;
}


void ConfFinal(string config_path){ //Write final configuration
    ofstream WriteConf;
    ofstream WriteOldConf;

    cout << "Print final configuration to file config.final and old.final" << endl << endl;
    string conf_final = config_path;
    string conf_final_old = config_path;

    conf_final.append("config.final");
    conf_final_old.append("old.final");
    WriteConf.open(conf_final);
    WriteOldConf.open(conf_final_old);

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

double Error(double sum, double sum2, int iblk)  //Algorithm for periodic boundary conditions with side L=box
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Averages(int iblk) //Print results for current block
{
    double r, gdir;
    ofstream AveEpot, AveEkin, AveTemp, AveEtot, AvePres;
    ofstream Gofr, Gave;

    glob_av[iv]+=blk_av[iv]/blk_norm;
    glob_av[ik]+=blk_av[ik]/blk_norm;
    glob_av[ie]+=blk_av[ie]/blk_norm;
    glob_av[it]+=blk_av[it]/blk_norm;
    glob_av[ip]+=blk_av[ip]/blk_norm;

    glob_av2[iv]+=pow(blk_av[iv]/blk_norm,2);
    glob_av2[ik]+=pow(blk_av[ik]/blk_norm,2);
    glob_av2[ie]+=pow(blk_av[ie]/blk_norm,2);
    glob_av2[it]+=pow(blk_av[it]/blk_norm,2);
    glob_av2[ip]+=pow(blk_av[ip]/blk_norm,2);

    AveEpot.open("data/ave_epot.out",ios::app);
    AveEkin.open("data/ave_ekin.out",ios::app);
    AveEtot.open("data/ave_etot.out",ios::app);
    AveTemp.open("data/ave_temp.out",ios::app);
    AvePres.open("data/ave_pres.out",ios::app);
    Gofr.open("data/gr.out",ios::app);
    Gave.open("data/ave_gr.out",ios::app);

    AveEpot<< iblk*nstep <<"   "<< glob_av[iv]/(double)iblk <<"   "<< Error(glob_av[iv], glob_av2[iv], iblk) << endl;
    AveEkin<< iblk*nstep <<"   "<< glob_av[ik]/(double)iblk <<"   "<< Error(glob_av[ik], glob_av2[ik], iblk) << endl;
    AveEtot<< iblk*nstep <<"   "<< glob_av[ie]/(double)iblk <<"   "<< Error(glob_av[ie], glob_av2[ie], iblk) << endl;
    AveTemp<< iblk*nstep <<"   "<< glob_av[it]/(double)iblk <<"   "<< Error(glob_av[it], glob_av2[it], iblk) << endl;
    AvePres<< iblk*nstep <<"   "<< glob_av[ip]/(double)iblk <<"   "<< Error(glob_av[ip], glob_av2[ip], iblk) << endl;

    //g(r)
    for (int i=igofr; i<nbins+igofr; i++){

        r = (i-igofr) * bin_size;
        gdir = 1./(rho*(double)npart*4./3.*M_PI*(pow(r+bin_size, 3)-pow(r, 3)) )*blk_av[i]/blk_norm;

        glob_av[i] += gdir;
        glob_av2[i] += gdir*gdir;

        Gofr << r + bin_size/2.0 << " " << gdir << " " << endl;

        if(iblk == nblock) // print only the final result
        {
          Gave << r + bin_size/2.0 << " " << glob_av[i]/(double)iblk << " " << Error(glob_av[i],glob_av2[i],iblk) << endl;
        }
    }

    AveEpot.close();
    AveEkin.close();
    AveTemp.close();
    AveEtot.close();
    AvePres.close();
    Gofr.close();
    Gave.close();
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

    for(int i=0; i<n_props; ++i){
        blk_av[i] = 0;
    }
    blk_norm = 0;
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
