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
int n_gen, pop_size; // number of generation and size of the population in each generation
int n_mutations, n_crossover;

vector<point> Cities;

double pick_power;

double p_shift, p_shift_chunk, p_inversion_chunk, p_permutation;
double p_crossover;

int iprint;

double L1(vector<point>,vector<int>);
double L2(vector<point>,vector<int>);
vector<point> GenerateCities(int , double, int );

Random rng_load();

void Input(int);
void SaveCities(vector<point>, string);
individual GenerateIndividual(vector<point>);
void PopulationSort(vector<individual> &, int );
void RegenerateL1(individual &);
void RegenerateL2(individual &);
void RegenerateL1(vector<individual> &);
void RegenerateL2(vector<individual> &);

void SavePopulation(vector<individual> , string );
void PrintPopulation(vector<individual>,int);
int CheckIndividual(individual);
int CheckIndividual(vector<individual>);
int Select(double , int );
double Fitness(vector<individual> , int );


void quicksort_L1(vector<individual> &,int ,int );
void quicksort_L2(vector<individual> &,int ,int );

vector<individual> Mutate(vector<individual>,int);
// possible mutations, every one returns one individual
individual Permutation(individual,double);
individual Shift(individual,double);
individual ShiftChunk(individual,double);
individual InversionChunk(individual,double);

vector<individual> Crossover(individual, individual,double);


#endif
