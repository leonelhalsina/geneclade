#ifndef dynamic_model
#define dynamic_model

#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace Rcpp;
using namespace std;

// // Declaration specification for exported functions
// #if defined(_WIN32) || defined(__WINDOWS__)
// #define REXPORTS extern "C" __declspec(dllexport)
// #else
// #define REXPORTS extern "C"
// #endif

//ifstream mapp;
//std::string name1 = "map_medium.txt";

//default_random_engine &generator;
struct probabilities_based_traits {
  vector <double> total_probability_species;
  double gammas_total;
  double lambdas_total;
  double mus_total;
  double popchange_rate_total;
  double geneflow_rate_total;
  double mutation_rate_total;
  double total_waiting_times_rates;

};

struct yx {
  int y;
  int x;

};

class landscape {
public:
  int habitat;
  int Nspp;
  double sum_trait;
  double average_trait;
  int k_patch; //this k comes from the map
  int total_abundance_cell;
  int temperature;
  int elevation;
private:
};

class population_structure {
public:
  int pop_size;
  int allelic_a;
  int allelic_b;
  int allelic_c;
  int allelic_d;
  int allelic_e;

  int allelic_z;
  int allelic_y;
  int allelic_x;
  int allelic_w;
  int allelic_v;
  void happening_population_popchange(landscape **map1, double , vector<int> , yx  );
  void happening_mutation();

private:
};



class species {
public:
  bool alive;
  int id;
  int elevation_origin;
  int parent;
  int range = 1;
  int x_coordinate_last_event;
  int y_coordinate_last_event;
  int range_highlands = 0;
  int range_intermediate1 = 0;
  int range_intermediate2 = 0;
  int range_lowlands = 0;
  int succesful_geneflow_events = 0;
  int southernmost;
  int northernmost;
  vector <int> change_northernmost;
  vector <int> change_southernmost;
  vector <double> time_change_northernmost;
  vector <double> time_change_southernmost;
  yx birthplace;
  double total_rate;
  vector <double> temperature_optimum;
  vector <double> computed_rate_based_on_temperature;
  double birth;
  double saturation_grid_birth;
  double death = 0;
  double trait_state;
  int total_pop_size;
  vector <yx> presence;

  vector <population_structure> populations_this_species;
  void happening_population_popchange_this_species(double , landscape **map1, vector<int>  );

  void happening_gene_flow( double,landscape **map1);

  void happening_mutation_this_species() ;
  void happening_trait_evolution(double, double);
  vector <yx> available_neigh_to_colonize_K(int ,int ,int , landscape **map1);
  vector <yx> available_neigh_to_colonize_trait(int ,int,int, double,landscape **map1);
  void happening_expansion(int , int , bool , double , landscape **map1, std::string , vector<int>, double );
  void happening_contraction( double ,landscape **map1,std::string );
  void classify_elevation(bool,landscape **map1);
  void classify_elevation_origin(bool);
  void initial_position(yx);
  vector <yx> find_neighbor(int cell);
private:
};



double calculate_variance(vector<int>);
double give_me_random_normal(double, double);
int give_me_random_uniform(int , int );
vector <int> give_me_random_wallenius (int, vector <int>, vector <int> , int );

void show_all_species_data(vector<species> );

double link_fitnesslike_mu_gamma(yx , population_structure , vector<int> , landscape **map1) ;
void happening_speciation( vector<species>& , vector<int> , int , double , double,landscape **map1);

void set_landscape(IntegerVector,IntegerVector ,IntegerVector,int , int,landscape **map1);
void change_temperature_map(int,int,IntegerVector ,IntegerVector ,landscape **map1);
vector <double> extract_temperature_func ( species);
vector <species> get_species_intocpp(vector<species> , IntegerVector ,IntegerVector, IntegerVector , IntegerVector , IntegerVector , IntegerVector , IntegerVector , NumericVector , NumericVector , NumericVector , IntegerVector , int , landscape **map1, vector<int> , double , double , double , std::string  );
probabilities_based_traits calculate_probabilities_using_traitstate(vector<species> , landscape **map1, std::string ,double, double, double , double , double , double , vector<int> , double  );
void populate_landscape(vector <species> all_species,landscape **map1);
void to_show_richness_map(std::string,vector<species>,landscape **map1);
List extract_species_data(vector <species> process_all_species);
List get_me_output(vector <species>,double );
vector <int> extract_yx_func(species, std::string );
bool final_check(vector<species>,landscape **map1);
#endif
