#ifndef __tsp__
#define __tsp__

#include <vector>

using namespace std;

struct individual
{
  vector<int> city_order;
  double L1;
  double L2;

  individual(const individual&) = default ;
  individual(vector<int> a={-1}, double b=-1, double c=-1)
  : city_order(a), L1(b), L2(c)
  {}

  individual& operator=(const individual& a)
  {
    city_order=a.city_order;
    L1=a.L1;
    L2=a.L2;
    return *this;
  }

  bool operator==(const individual& a)
  {
    return(city_order == a.city_order);
  }
  bool operator!=(const individual& a)
  {
    return( !(*this == a));
  }
};

struct population
{
  vector<individual> orders_list;

  population(vector<individual> a): orders_list(a) {};
};

int n_cities, city_type, world_size; //how many cities has to be created and how
int n_mutations;
double t_min, t_max, temperature, beta;
int t_steps;
int iprint;

int attempt, accept;

vector<point> Cities;
double * city_x;
double * city_y;

int size;
int Rank;
int loaded = 0,loader;

Random rnd;

double L1(vector<point>,vector<int>);
double L2(vector<point>,vector<int>);
vector<point> GenerateCities(int , double, int );

void Input();
void SaveCities(vector<point>, string);
individual GenerateIndividual(vector<point>);

void SaveL(string,int,double);
void SaveBestPath(string , individual, vector<point>);
int CheckIndividual(individual);
int CheckIndividual(vector<individual>);
void PrintPath(individual);
void PrintParameters();

individual Mutate(individual ,int);
// possible mutations, every one returns one individual
individual Permutation(individual);
individual PermutationChunk(individual);
individual Shift(individual);
individual ShiftChunk(individual);
individual InversionChunk(individual);


#endif
