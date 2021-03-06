#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include "rng/random.h"
#include "tsp.h"
#include <mpi.h>

using namespace std;


int main()
{
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

  Input();
  //PrintParameters();
  // Create the initial individual, with a random path between the cities

  individual man = GenerateIndividual(Cities);

  individual xman; //mutated individual

  for(int t=0; t <= t_steps; t++)
  {
    accept=0;

    temperature=t_max-t*(t_max-t_min)*1./t_steps;
    beta=1./temperature;

    xman=Mutate(man,n_mutations);

    if(Rank== 0  && t%iprint==0)
    {
      cout << endl;
      cout << " Temperature : " << temperature  << endl;
      // l1 and the accept rate are per-node
      //cout << "L1 = " << xman.L1 << endl;
      //cout << "Mutations accept rate = " << accept*1./n_mutations << endl << endl;
    }
    man=xman;
    SaveL("distance_" + to_string(Rank) + ".out",t,man.L1);
  }

  // print the best path in the last generation
  SaveBestPath("best_path_mpi.out",man,Cities);

  MPI_Finalize();
  return 0;
}

/*************************************************/
/*********** Input/output funtions ***************/
/*************************************************/

void Input()
{
  // random generator loading
  int seed[4];
  int *p1=new int[size];
  int *p2=new int[size];

  if (Rank==0)
  {
    ifstream Primes("rng/Primes");
    if (Primes.is_open())
    {
      for(int i=0; i<size; i++) //change the primes for the rnd ofr each node
      {
        Primes >> p1[i] >> p2[i] ;
      }
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
          }
       }
       input.close();
    }
    else cerr << "PROBLEM: Unable to open seed.in" << endl;
    rnd.SaveSeed();
  }

  MPI_Bcast(p1, size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(p2, size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(seed, 4, MPI_INT, 0, MPI_COMM_WORLD);

  // each node has its own set of primes
  rnd.SetRandom(seed, p1[Rank], p2[Rank]);

  // read parameters from input.dat
  if(Rank==0)
  {
    ifstream ReadInput;

    cout << "Trade Salesman Problem - solved with Simulated Annealing" << endl;
    cout << "Solving with MPI, using " << size << " nodes."<< endl;
    cout << "Node " << Rank << " loading the parameters" << endl;

    ReadInput.open("input.dat");

    ReadInput >> n_cities;
    ReadInput >> city_type;
    ReadInput >> world_size;
    ReadInput >> iprint;
    ReadInput >> t_min;
    ReadInput >> t_max;
    ReadInput >> t_steps;
    ReadInput >> n_mutations;

    cout << "Cities to visit = " << n_cities <<endl;
    if(city_type == 1)
      cout << "Cities are placed on the edge of a circumference of radius = " << world_size <<endl;
    else
      cout << "Cities are placed inside a square of side = " << world_size << endl;
    cout << "Temperature range = ( " << t_min << " , " << t_max << " )" << endl;
    cout << "Number of temperature step = " << t_steps << endl;
    cout << "Mutation attempts per step = " << n_mutations << endl;
  }

  //PrintParameters();

  MPI_Bcast(&iprint, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&t_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&t_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&t_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_mutations, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_cities, 1, MPI_INT, 0, MPI_COMM_WORLD);

  city_x=new double[n_cities];
  city_y=new double[n_cities];
  Cities.resize(n_cities);

  if(Rank==0)
  {
    cout << "Node " << Rank << " generating cities positions" << endl;
    // create and save cities positions to file
    Cities = GenerateCities(n_cities,world_size,city_type);
    string cities_conf = "data/cities.dat";
    cout << "Cities coordinates saved in " << cities_conf << endl;
    SaveCities(Cities,cities_conf);

    for( int i = 0; i < n_cities; i++)
    {
      city_x[i]=Cities[i].x;
      city_y[i]=Cities[i].y;
    }
  }

  MPI_Bcast(city_x, n_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(city_y, n_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //PrintParameters();

  //load the coordinates on the other nodes
  for(int i=0; i<n_cities; i++)
  {
      Cities[i].x=city_x[i];
      Cities[i].y=city_y[i];
  }

  delete [] city_x;
  delete [] city_y;
  delete [] p1;
  delete [] p2;

}

void PrintPath(individual jeff)
{
  cout << jeff.L1;
  for(unsigned j = 0; j<jeff.city_order.size(); j++)
  {
    cout << " " << jeff.city_order[j];
  }
  cout << endl;
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

void SaveL(string filename, int step, double L)
{
  ofstream l_file;
  l_file.open(filename,ios::app);
  l_file << step << " " << L << endl;
}

void SaveBestPath(string filename, individual man, vector<point> city_list)
{
  // gather the L1 from all the nodes, and find the shortest
  double final_l = man.L1;
  double *best_l=new double[size];

  MPI_Gather(&final_l, 1, MPI_DOUBLE, best_l, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int node;
  if(Rank==0)
  {
    double l_min;
    cout << endl << "Comparing path of different nodes:" << endl << endl;
    for(int i=0; i<size; i++)
    {
      cout << "L for node " << i << " = " << best_l[i] << endl;
      if(i==0 || best_l[i] < l_min)
      {
        node=i;
        l_min=best_l[i];
      }
    }
    cout<<"Best path found on node "<<node<<endl;
  }
  MPI_Bcast(&node, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);

  // the node with the best solution saves it
  if(Rank == node)
  {
    ofstream c_file;
    c_file.open(filename);
    for(unsigned i=0; i< city_list.size(); i++)
    {
      c_file << man.city_order[i] << " " << city_list[man.city_order[i]].x << " " << city_list[man.city_order[i]].y<<endl;
    }
    // add first city at the end
    c_file << man.city_order[0] << " " << city_list[man.city_order[0]].x << " " << city_list[man.city_order[0]].y << endl;
    c_file.close();
  }
  delete [] best_l;
}

/*****************************************************************/
/************************* MUTATIONS *****************************/
/*****************************************************************/

individual Mutate(individual old_path, int mutation_attempts)
{
  individual new_path;

  for( int j = 0; j < mutation_attempts; j++)
  {
    // pick a mutation
    int mutation_id = rnd.Rannyu(0,5);
    switch (mutation_id)
    {
      case 0:
        new_path = Shift(old_path);
        //cout << "Shift" << endl;
        break;
      case 1:
        new_path = ShiftChunk(old_path);
        //cout << "ShiftChunk" << endl;
        break;
      case 2:
        new_path = Permutation(old_path);
        //cout << "Permutation" << endl;
        break;
      case 3:
        new_path = PermutationChunk(old_path);
        //cout << "PermutationChunk" << endl;
        break;
      case 4:
        new_path = InversionChunk(old_path);
        //cout << "InversionChunk" << endl;
        break;
      default:
        cerr << "Unrecognized mutation" << endl;
        break;
    }
    double p = exp(-beta*(new_path.L1-old_path.L1));
    if(rnd.Rannyu() < p)
    {

        old_path=new_path;
        accept++;
    }
  }
  return old_path;
}

individual Permutation(individual a)
{
  individual new_a = a;

  //cout << "Permutation "<< endl;
  int j=rnd.Rannyu(0, a.city_order.size());
  int k=rnd.Rannyu(0, a.city_order.size());
  int tmp = new_a.city_order[j];
  new_a.city_order[j] = new_a.city_order[k];
  new_a.city_order[k] = tmp;
  //cout << a.city_order.size() << ":" << new_a.city_order.size() << endl;

  int e = CheckIndividual(new_a);
  if(e!=0)
    cout << "Permutation error: " << e << endl;
  new_a.L1=L1(Cities,new_a.city_order);
  new_a.L2=L2(Cities,new_a.city_order);

  return new_a;
}

individual PermutationChunk(individual a)
{
  individual new_a = a;
  int g1,g2,m;

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

individual Shift(individual a)
{
  individual new_a=a;

  //cout << "Shift "<< endl;
  int j=rnd.Rannyu(0, a.city_order.size());
  for(unsigned i = 0; i < a.city_order.size(); i++)
  {
    new_a.city_order[i]=a.city_order[(i+j)%a.city_order.size()];
  }
  //cout << "shift" << a.city_order.size() << ":" << new_a.city_order.size() << endl;

  int e = CheckIndividual(new_a);
  if(e!=0)
    cout << "Shift error: " << e << endl;
  new_a.L1=L1(Cities,new_a.city_order);
  new_a.L2=L2(Cities,new_a.city_order);
  return new_a;
}

individual ShiftChunk(individual a)
{
  individual new_a = a;
  int g,m,n;

  //cout << "ShiftChunk "<< endl;
  g=rnd.Rannyu(0, a.city_order.size()); // start of the chunk to shift
  m=rnd.Rannyu(1, a.city_order.size()-g); // width of the chunk to shift
  n=rnd.Rannyu(1, a.city_order.size()-g-m); // shift lenght
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

  int e = CheckIndividual(new_a);
  if(e!=0)
  {
    cout << "ShiftChunk error: " << e << endl;
    cout << g << " " << m  << " " << n << endl;
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


individual InversionChunk(individual a)
{
  individual new_a = a;
  int g,m;

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
  return new_a;
}


/***************************************************************************/
/******************** Parameters evaluation ********************************/
/***************************************************************************/

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

int CheckIndividual(individual a)
{
  vector<int> count = {};
  int error =0;

  if(int(a.city_order.size()) != n_cities)
  {
    cout << "ERROR: Number of cities visited doesn't match the total number of cities:";
    cout << a.city_order.size() <<" "<< n_cities << endl;
    error += n_mutations*1000;
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

void PrintParameters()
{
  cout << "Parameters for node " << Rank << " :"<< endl;
  cout << t_min << " " << t_max << " " << t_steps << endl;
  cout << n_mutations << endl;
  cout << iprint << endl;
}
