#include <cmath>
#include <stdlib.h>     /* atof */
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include "rng/random.h"
#include "tsp.h"

using namespace std;

Random rnd = rng_load();

int main(int argc, char *argv[])
{
  Input(argc); //Inizialization

  // create and save cities positions
  Cities = GenerateCities(n_cities,world_size,city_type);
  string cities_conf = "data/cities.dat";
  cout << "Cities coordinates are saved in " << cities_conf << endl;
  SaveCities(Cities,cities_conf);

  // Create the initial population
  // - An individual is identified by a vector of the city visit order,
  //   and two double, that are its L(1) and L(2) values

  // - A population is a vector of individuals
  vector<individual> starting_population, new_population;
  for(int i =0; i< pop_size; i++)
  {
    individual dave = GenerateIndividual(Cities); // the visit order is purely random
    starting_population.push_back(dave);
  }

  /*
  // debug purposes, can be uncommented without harming the code
  vector<int> ordered = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
  individual kevin(ordered,L1(Cities,ordered),L2(Cities,ordered));
  starting_population.pop_back();
  starting_population.push_back(kevin);
  */

  //PrintPopulation(starting_population,1);
  PopulationSort(starting_population,1);
  SavePopulation(starting_population,"data/start_population.dat");

  new_population=starting_population;

  vector<individual> xmen; //mutated individuals
  vector<individual> new_generation; // crossover result

  for(int i = 0; i < n_gen; i++)
  {
    //single individual mutation
    xmen=Mutate(new_population,n_mutations);

    //crossover
    for( int j = 0; j < n_crossover; j++)
    {
      int id_m = Select(pick_power,pop_size);
      int id_f = id_m;
      while(id_f == id_m) // avoid self-crossover
      {
         id_f = Select(pick_power,pop_size);
      }
      vector<individual> offspring = Crossover(new_population[id_m],new_population[id_f],p_crossover);
      if(offspring[0]!=new_population[id_m] && offspring[0]!=new_population[id_f])
      {
        new_generation.push_back(offspring[0]);
        new_generation.push_back(offspring[1]);
      }
    }

    if(i%iprint==0)
    {
      cout << endl;
      cout << "Generation : " << i+1 << "/" << n_gen << endl;
      cout << "Average fitness = " << Fitness(new_population,1) << endl;
      cout << "Best fitness = " << new_population[0].L1 << endl;
      cout << "Mutations occurred = " << xmen.size() << endl;
      cout << "Crossover occurred = " << new_generation.size() << endl;
    }

    // join the individuals and double check for errors
    new_population.reserve(new_population.size() + xmen.size() + new_generation.size()); // preallocate memory
    new_population.insert(new_population.end(),xmen.begin(), xmen.end());
    new_population.insert(new_population.end(),new_generation.begin(), new_generation.end());
    // clean memory
    xmen.clear();
    new_generation.clear();

    int errors = 0;
    for(unsigned ii = 0; ii< new_population.size();ii++)
    {
      int e = CheckIndividual(new_population[ii]);
      errors+=e;
      if(e!=0)
      {
        cout << "Anomalous individual : " << ii << "Error = " << e << endl;
        cout << new_population[ii].L1 << " " << new_population[ii].L2 << " ";
        for(unsigned c = 0; c < new_population[ii].city_order.size();c++)
        {
          cout <<new_population[ii].city_order[c] << " " << endl;
        }
      }
    }
    if(errors != 0)
      cout <<"Total errors  =" << errors << endl << endl;

    PopulationSort(new_population,1);
    // mantain the population at the same size, keeping only the fittest.
    while(int(new_population.size()) > pop_size)
    {
      new_population.pop_back();
    }

    // save best L and average L (on the best half) for the current generation
    ofstream l_file;
    l_file.open("fitness.out",ios::app);
    l_file << i << " " << new_population[0].L1 << " " << Fitness(new_population,1) << endl;
    l_file.close();
  }

  // print the best path in the last generation
  ofstream c_file;
  c_file.open("best_path.out");
  for(unsigned i=0; i< Cities.size(); i++)
  {
    c_file << new_population[0].city_order[i] << " " << Cities[new_population[0].city_order[i]].x << " " << Cities[new_population[0].city_order[i]].y<<endl;
  }
  // add first city at the end
  c_file << new_population[0].city_order[0] << " " << Cities[new_population[0].city_order[0]].x << " " << Cities[new_population[0].city_order[0]].y << endl;
  c_file.close();

  return 0;
}

/*************************************************/
/*********** Input/output funtions ***************/
/*************************************************/

void Input(int args){

    ifstream ReadInput;

    cout << "Trade Salesman Problem" << endl;

    ReadInput.open("input.dat");

    ReadInput >> n_cities;
    ReadInput >> city_type;
    ReadInput >> world_size;
    ReadInput >> n_gen;
    ReadInput >> pop_size;
    ReadInput >> iprint;
    ReadInput >> n_mutations;
    ReadInput >> n_crossover;

    // mutation and crossover attempt match the population size
    if(n_mutations < 0)
      n_mutations = pop_size;
    if(n_crossover < 0)
      n_crossover = pop_size;
    ReadInput >> pick_power;
    ReadInput >> p_permutation;
    ReadInput >> p_permutation_chunk;
    ReadInput >> p_shift;
    ReadInput >> p_shift_chunk;
    ReadInput >> p_inversion_chunk;
    ReadInput >> p_crossover;

    cout << "Cities to visit = " << n_cities <<endl;
    if(city_type == 1)
      cout << "Cities are placed on the edge of a circumference of radius = " << world_size <<endl;
    else
      cout << "Cities are placed inside a square of side = " << world_size << endl;
    cout << "Number of generations = " << n_gen << endl;
    cout << "Population size = " << pop_size << endl;
    cout << "Mutation attempts per generation = " << pop_size << endl;
    cout << "Crossover attempts per generation = " << pop_size << endl;
}

void PrintPopulation(vector<individual> population, int v)
{
  if(v==1) // print the fitness and the city order, each individual per line
  {
    for(unsigned i = 0; i< population.size() ;i++)
    {
      cout << population[i].L1;
      for(unsigned j = 0; j<population[i].city_order.size(); j++)
      {
        cout << " " << population[i].city_order[j];
      }
      cout << endl;
    }
  }
  if(v==0) // print only fitness, all in the same line
  {
    for(unsigned i = 0; i< population.size() ;i++)
    {
      cout << population[i].L1 << " ";
    }
    cout << endl;
  }
}

void SaveCities(vector<point> cities, string filename)
{
  ofstream c_file;
  c_file.open(filename);
  for(unsigned i=0; i< cities.size(); i++)
  {
    c_file << cities[i].x << " " << cities[i].y <<endl;
  }
  c_file.close();
}

void SavePopulation(vector<individual> population, string filename)
{
  ofstream p_file;
  p_file.open(filename);
  p_file << "L1 " << "L2 " << "City order"<< endl;
  for(unsigned i = 0; i<population.size(); i++)
  {
    p_file << population[i].L1 << " " << population[i].L2;

    for(unsigned j =0; j< population[i].city_order.size();j++)
    {
      p_file << " " << population[i].city_order[j];
    }
    p_file << endl;

  }
  p_file.close();
}

/*****************************************************************/
/************************* MUTATIONS *****************************/
/*****************************************************************/

vector<individual> Mutate(vector<individual> population, int mutation_attempts)
{
  vector<individual> xmen;

  for( int j = 0; j < mutation_attempts; j++)
  {
    // pick a (maybe) different individual for each mutation.
    int id_permutation = Select(pick_power,pop_size);
    int id_permutation_chunk = Select(pick_power,pop_size);
    int id_shift = Select(pick_power,pop_size);
    int id_shift_chunk = Select(pick_power,pop_size);
    int id_inversion_chunk = Select(pick_power,pop_size);

    // each mutation try to generate a new individual
    // if the mutation happens, the individual it's added to the population

    xmen.push_back(Permutation(population[id_permutation],p_permutation));
    if(xmen.back()==population[id_permutation])
      xmen.pop_back();

    xmen.push_back(PermutationChunk(population[id_permutation_chunk],p_permutation_chunk));
    if(xmen.back()==population[id_permutation_chunk])
      xmen.pop_back();

    xmen.push_back(Shift(population[id_shift],p_shift));
    if(xmen.back()==population[id_shift])
      xmen.pop_back();

    xmen.push_back(ShiftChunk(population[id_shift_chunk],p_shift_chunk));
    if(xmen.back()==population[id_shift_chunk])
      xmen.pop_back();

    xmen.push_back(InversionChunk(population[id_inversion_chunk],p_inversion_chunk));
    if(xmen.back()==population[id_inversion_chunk])
      xmen.pop_back();
  }

  return xmen;
}

individual Permutation(individual a, double p)
{
  individual new_a = a;
  if(rnd.Rannyu() <= p)
  {
    //cout << "Permutation "<< endl;
    int j=rnd.Rannyu(0, a.city_order.size());
    int k=rnd.Rannyu(0, a.city_order.size());
    int tmp = new_a.city_order[j];
    new_a.city_order[j] = new_a.city_order[k];
    new_a.city_order[k] = tmp;
    //cout << a.city_order.size() << ":" << new_a.city_order.size() << endl;
  }
  int e = CheckIndividual(new_a);
  if(e!=0)
    cout << "Permutation error: " << e << endl;
  new_a.L1=L1(Cities,new_a.city_order);
  new_a.L2=L2(Cities,new_a.city_order);

  return new_a;
}

individual PermutationChunk(individual a, double p)
{
  individual new_a = a;
  int g1,g2,m;
  if(rnd.Rannyu() <= p)
  {
    //cout << "ShiftChunk "<< endl;
    g1=rnd.Rannyu(0, a.city_order.size()); // start of the first
    m=rnd.Rannyu(0, a.city_order.size()*0.5); // length of the chunks
    g2=rnd.Rannyu(g1+m+1, g1 - m + a.city_order.size()); // start of the second chunk
    g2 = g2%a.city_order.size();

    for(int i = 0; i < m; i++)
    {
      new_a.city_order[(g1+i)%a.city_order.size()] = a.city_order[(g2+i)%a.city_order.size()];
      new_a.city_order[(g2+i)%a.city_order.size()] = a.city_order[(g1+i)%a.city_order.size()];
    }
  }

  int e = CheckIndividual(new_a);
  if(e!=0)
  {
    cout << "PermutationChunk error: " << e << endl;
    cout << g1 << " " << g2 << " " << m << endl;
    for(unsigned i = 0; i < a.city_order.size(); i++)
    {
      cout << a.city_order[i] << " ";
    }
    cout << endl;
    for(unsigned i = 0; i < new_a.city_order.size(); i++)
    {
      cout << new_a.city_order[i] << " ";
    }
    cout << endl;

  }
  new_a.L1=L1(Cities,new_a.city_order);
  new_a.L2=L2(Cities,new_a.city_order);
  return new_a;
}

individual Shift(individual a, double p )
{
  individual new_a=a;
  if(rnd.Rannyu() <= p)
  {
    //cout << "Shift "<< endl;
    int j=rnd.Rannyu(0, a.city_order.size());

    for(unsigned i = 0; i < a.city_order.size(); i++)
    {
      new_a.city_order[i]=a.city_order[(i+j)%a.city_order.size()];
    }
    //cout << "shift" << a.city_order.size() << ":" << new_a.city_order.size() << endl;

  }
  int e = CheckIndividual(new_a);
  if(e!=0)
    cout << "Shift error: " << e << endl;
  new_a.L1=L1(Cities,new_a.city_order);
  new_a.L2=L2(Cities,new_a.city_order);
  return new_a;
}

individual ShiftChunk(individual a, double p)
{
  individual new_a = a;
  int g,m,n;
  if(rnd.Rannyu() <= p)
  {
    //cout << "ShiftChunk "<< endl;
    g=rnd.Rannyu(0, a.city_order.size()); // start of the chunk to shift
    m=rnd.Rannyu(0, a.city_order.size()-g); // width of the chunk to shift
    n=rnd.Rannyu(0, a.city_order.size()-g-m); // shift lenght

    if(n>0)
    {
      for(int i = 0; i < m; i++)
      {
        new_a.city_order[(g+n+i)%a.city_order.size()] = a.city_order[(g+i)%a.city_order.size()];
      }
      for(int i = 0; i < n ; i++)
      {
        new_a.city_order[(g+i)%a.city_order.size()] = a.city_order[(g+m+i)%a.city_order.size()];
      }
    }
  }
  int e = CheckIndividual(new_a);
  if(e!=0)
    cout << "ShiftChunk error: " << e << endl;
  new_a.L1=L1(Cities,new_a.city_order);
  new_a.L2=L2(Cities,new_a.city_order);
  return new_a;
}


individual InversionChunk(individual a, double p)
{
  individual new_a = a;
  int g,m;
  if(rnd.Rannyu() <= p)
  {
    //cout << "InversionChunk "<< endl;
    g=rnd.Rannyu(0, a.city_order.size());
    m=rnd.Rannyu(0, a.city_order.size());

    for(int i = 0; i <= m; i++)
    {
      new_a.city_order[(g+i)%a.city_order.size()] = a.city_order[(g-i+m)%a.city_order.size()];
    }

    int e = CheckIndividual(new_a);
    if(e!=0)
    {
      cout << "InversionChunk error: " << e << endl;
      cout << g << " " << m << endl;
      for(unsigned i = 0; i < a.city_order.size(); i++)
      {
        cout << a.city_order[i] << " ";
      }
      cout << endl;
      for(unsigned i = 0; i < new_a.city_order.size(); i++)
      {
        cout << new_a.city_order[i] << " ";
      }
      cout << endl;

    }
    new_a.L1=L1(Cities,new_a.city_order);
    new_a.L2=L2(Cities,new_a.city_order);
  }

  return new_a;
}

vector<individual> Crossover(individual m, individual f, double p)
{
  vector<individual> s={m,f};
  vector<int> tail0, tail1;

  if(rnd.Rannyu() <= p)
  {
    int g = rnd.Rannyu(0, m.city_order.size());

    tail0 = s[0].city_order;
    tail1 = s[1].city_order;

    //keep the first g elements
    s[0].city_order.resize(g);
    s[1].city_order.resize(g);

    // erase the first g elements
    tail0.erase(tail0.begin(),tail0.begin()+g);
    tail1.erase(tail1.begin(),tail1.begin()+g);

    for(unsigned j=0; j < m.city_order.size(); j++)
    {
      for(unsigned i = 0 ; i < tail0.size(); i++)
      {
        if(f.city_order[j] == tail0[i])
          s[0].city_order.push_back(f.city_order[j]);
        if(m.city_order[j] == tail1[i])
          s[1].city_order.push_back(m.city_order[j]);
      }
    }
  }

  int e = CheckIndividual(s);
  if(e!=0)
    cout << "Crossover error: " << e << endl;

  s[0].L1=L1(Cities,s[0].city_order);
  s[1].L1=L1(Cities,s[1].city_order);
  s[0].L2=L2(Cities,s[0].city_order);
  s[1].L2=L2(Cities,s[1].city_order);

  return s;
}


/***************************************************************************/
/******************** Parameters evaluation ********************************/
/***************************************************************************/

int Select(double p, int population)
{
  return pow(rnd.Rannyu(), p)*population;
}

double L1(vector<point> cities,vector<int> visit_order)
{
  double L=0;
  point d;

  for(unsigned i=0; i < visit_order.size();i++)
  {
    d = cities[visit_order[i]] - cities[visit_order[(i+1)%visit_order.size()]];
    L += sqrt(d.x*d.x + d.y*d.y);
  }
  return L;
}

double L2(vector<point> cities,vector<int> visit_order)
{
  double L=0;
  point d;

  for(unsigned i=0; i < visit_order.size();i++)
  {
    d = cities[visit_order[i]] - cities[visit_order[(i+1)%visit_order.size()]];
    L += d.x*d.x + d.y*d.y;
  }
  return L;
}

void RegenerateL1(individual &a)
{
  a.L1 = L1(Cities,a.city_order);
}

void RegenerateL1(vector<individual> &population)
{
  for(unsigned i = 0; i< population.size(); i++)
  {
    RegenerateL1(population[i]);
  }
}

void RegenerateL2(individual &a)
{
  a.L2= L2(Cities,a.city_order);
}

void RegenerateL2(vector<individual> &population)
{
  for(unsigned i = 0; i< population.size(); i++)
  {
    RegenerateL2(population[i]);
  }
}

double Fitness(vector<individual> population, int mode)
{
  double l_ave = 0;
  int dim = 0.5*population.size();
  if(mode == 1)
    for(int i = 0; i < dim; i++)
      l_ave += population[i].L1;
  if(mode == 1)
    for(int i = 0; i< dim; i++)
      l_ave += population[i].L2;

  return l_ave*1./dim;
}

int CheckIndividual(individual a)
{
  vector<int> count = {};
  int error =0;

  if(int(a.city_order.size()) != n_cities)
  {
    cout << "ERROR: Number of cities visited doesn't match the total number of cities:";
    cout << a.city_order.size() <<" "<< n_cities << endl;
    error += pop_size*1000;
  }
  else
  {
    for(unsigned i = 0; i< a.city_order.size(); i++)
    {
      for(unsigned j = 0; j<count.size(); j++)
      {
        if(a.city_order[i] == count[j])
        {
          cout <<"ERROR: Duplicated city " << a.city_order[i] << endl;
          error += 1;
        }
      }
      count.push_back(a.city_order[i]);
    }
  }
  return error;
}

int CheckIndividual(vector<individual> a)
{
  int error = 0;
  for(unsigned i = 0; i< a.size(); i++)
  {
    error += CheckIndividual(a[i]);
  }
  return error;
}


void PopulationSort(vector<individual> &population, int mode)
{
  if(mode==1)
    quicksort_L1(population,0,population.size()-1);
  if(mode ==2)
    quicksort_L2(population,0,population.size()-1);
}

void quicksort_L1(vector<individual> &population,int primo,int ultimo)
{
   int i, j, pivot;
   individual temp;

   if(primo<ultimo)
   {
     pivot=primo;
     i=primo;
     j=ultimo;

     while(i<j)
     {
        while(population[i].L1<=population[pivot].L1&&i<ultimo)
          i++;
        while(population[j].L1>population[pivot].L1)
          j--;
        if(i<j)
        {
          temp=population[i];
          population[i]=population[j];
          population[j]=temp;
        }
    }

    temp=population[pivot];
    population[pivot]=population[j];
    population[j]=temp;
    quicksort_L1(population,primo,j-1);
    quicksort_L1(population,j+1,ultimo);
  }
}

void quicksort_L2(vector<individual> &population,int primo,int ultimo)
{
   int i, j, pivot;
   individual temp;

   if(primo<ultimo)
   {
     pivot=primo;
     i=primo;
     j=ultimo;

     while(i<j)
     {
        while(population[i].L1<=population[pivot].L1&&i<ultimo)
          i++;
        while(population[j].L1>population[pivot].L1)
          j--;
        if(i<j)
        {

          temp=population[i];
          population[i]=population[j];
          population[j]=temp;
        }
    }

    temp=population[pivot];
    population[pivot]=population[j];
    population[j]=temp;
    quicksort_L2(population,primo,j-1);
    quicksort_L2(population,j+1,ultimo);
  }
}

/*************************************************************************/
/************************* Initial generation ****************************/
/*************************************************************************/

individual GenerateIndividual(vector<point> cities)
{
  vector<int> city_id, visit_list;

  for( unsigned i = 0; i< cities.size();i++) // int array of cities id
  {
    city_id.push_back(i);
  }
  // sample without dulicates the city id, in random order
  while(visit_list.size() < cities.size())
  {
    int id = rnd.Rannyu(0,city_id.size());
    visit_list.push_back(city_id[id]);
    city_id.erase(city_id.begin() + id);
  }

  individual jeff(visit_list,L1(cities,visit_list),L2(cities,visit_list));
  return jeff;
}


vector<point> GenerateCities(int n, double range, int generation_type)
{
  vector<point> cities;

  if(generation_type == 1) //generate on the circle
  {
    double R = range;
    double theta;
    for(int i = 0; i< n ; i++)
    {
      theta = rnd.Rannyu(0,2.*M_PI);
      point p(R*cos(theta),R*sin(theta),0);
      cities.push_back(p);
    }

  }
  else if(generation_type == -1) // ordered cities for debug
  {
    double R = range;
    double theta;
    for(int i = 0; i< n ; i++)
    {
      theta = M_PI *i*2./n;
      point p(R*cos(theta),R*sin(theta),0);
      cities.push_back(p);
    }
  }
  else //generation inside a square
  {
    double side = range;
    for(int i = 0; i< n ; i++)
    {
      double x = rnd.Rannyu(-0.5*side,0.5*side);
      double y = rnd.Rannyu(-0.5*side,0.5*side);
      point p(x,y,0);
      cities.push_back(p);
    }

  }
  return cities;
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
