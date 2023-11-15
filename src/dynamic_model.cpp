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
#include "dynamic_model.h"
#include "stocc.h"
using namespace Rcpp;
using namespace std;

// ifstream mapp;
// std::string name1 = "map_medium.txt";

const int x_max = 22;
const int y_max = 22;
// vector <species> all_species;


void species::update_latitudinal_borders(double t, bool was_it_expansion){


  vector <int> all_ys;
  for(int ij = 0; ij < presence.size();ij ++){

    all_ys.push_back(presence[ij].y);
  }

  int fromhere_southmost;
  fromhere_southmost = *max_element(all_ys.begin(),all_ys.end());
  int fromhere_northmost;
  fromhere_northmost = *min_element(all_ys.begin(),all_ys.end());


  if(was_it_expansion){
    if(fromhere_southmost > southernmost){
      change_southernmost.push_back(1);
      time_change_southernmost.push_back(t);
      southernmost = fromhere_southmost;
    }
    if(fromhere_northmost < northernmost){
      change_northernmost.push_back(1);
      time_change_northernmost.push_back(t);
      northernmost = fromhere_northmost;

    }

  } else {
    // cout << "contract" << endl;
    // cout << fromhere_southmost << "  " << southernmost << endl;
    // cout << fromhere_northmost << "  " <<northernmost << endl;
    if(fromhere_southmost < southernmost){
      change_southernmost.push_back(-1);
      time_change_southernmost.push_back(t);
      southernmost = fromhere_southmost;
      //cout << "contract south" << endl;
    }


    if(fromhere_northmost > northernmost){
      change_northernmost.push_back(-1);
      time_change_northernmost.push_back(t);
      northernmost = fromhere_northmost;
     // cout << "contract north" << endl;
    }
  }

}

double calculate_variance(vector<int> vector_values){

  double the_sum;
  the_sum = std::accumulate(vector_values.begin(),vector_values.end(),
                            decltype(vector_values)::value_type(0));
  double mean;
  mean = the_sum/vector_values.size();

  double variance;
  variance = 0;
  for (int i = 0; i < vector_values.size(); ++i) {
    variance = variance + ((vector_values[i] - mean) * (vector_values[i] - mean));
  }
  variance = variance/(vector_values.size() - 1);
  return variance;
}

void change_temperature_map(int x_max, int y_max, IntegerVector map_temperature_vector,IntegerVector ice_age_change,landscape **map1){

  vector <int> temperature_bands1;
  for(int ij = 0; ij < map_temperature_vector.size();ij ++)
  {
    if(map_temperature_vector[ij] > 0)
    {
      temperature_bands1.push_back(map_temperature_vector[ij]);
    }
  }

  vector<int>::iterator ip;
  // Sorting the array
  std::sort(temperature_bands1.begin(), temperature_bands1.end());
  // Now vector becomes 1 1 2 2 3 3 3 3 7 7 8 10
  // Using std::unique
  ip = std::unique(temperature_bands1.begin(), temperature_bands1.begin() + temperature_bands1.size());
  // Now vector becomes {1 2 3 7 8 10 * * * * * *}
  // * means undefined
  // Resizing the vector so as to remove the undefined terms
  temperature_bands1.resize(std::distance(temperature_bands1.begin(), ip));

  if(ice_age_change.size() != temperature_bands1.size() )
  {
    stop("ice_age_change needs to be same length than the number of temperature bands");
  }

  for (int i = 0; i < x_max; i++)
  {
    for (int j = 0; j < y_max; j++)
    {

      int this_temperature;
      this_temperature = map1[i][j].temperature;
      if(this_temperature > 0){
        int identi_band;
        for(int jj = 0;  jj < temperature_bands1.size(); jj++)
        {
          if(this_temperature == temperature_bands1[jj])
          {
            identi_band = jj;
          }
        }
        map1[i][j].temperature = (this_temperature + ice_age_change[identi_band]);

      }


    }
  }
}




bool final_check(vector<species> all_species,landscape **map1)
{
  cout << "doing all sort of checks" << endl;
  bool no_failure;
  no_failure = true;
  bool message_final_check_to_do;
  message_final_check_to_do = true;


  // create a map

  landscape **map_for_check = new landscape *[y_max];

  // dynamically allocate memory of size N for each row
  for (int i = 0; i < y_max; ++i)
  {
    map_for_check[i] = new landscape[x_max];
  }


  for (int i = 0; i < y_max; i++)
  {
    for (int j = 0; j < x_max; j++)
    {
      map_for_check[i][j].total_abundance_cell = 0;
    }
  }



  for(int iji = 0; iji < all_species.size(); ++iji)
  {
    species do_this_species;
    do_this_species = all_species[iji];
    if(do_this_species.alive)
    {

      int total_pop_size_from_individual_populations;
      total_pop_size_from_individual_populations = 0;
      for(int iij = 0; iij < do_this_species.presence.size(); ++iij)
      {

        int current_abundance_this_cell_from_map;

        current_abundance_this_cell_from_map = map_for_check[do_this_species.presence[iij].y - 1 ][do_this_species.presence[iij].x - 1].total_abundance_cell;

        map_for_check[do_this_species.presence[iij].y - 1 ][do_this_species.presence[iij].x - 1].total_abundance_cell = current_abundance_this_cell_from_map + do_this_species.populations_this_species[iij].pop_size;




        population_structure do_this_population;
        do_this_population = do_this_species.populations_this_species[iij];

        if(do_this_species.presence.size() != do_this_species.range)
        {
          no_failure = false;
          cout << " ++++++++++++++++++Something went wrong: issue 1 +++++++++++" << endl;
        }
        int sum_allelic;
        int sum_allelic_neutral;
        sum_allelic = do_this_population.allelic_a + do_this_population.allelic_b + do_this_population.allelic_c + do_this_population.allelic_d + do_this_population.allelic_e;
        sum_allelic_neutral = do_this_population.allelic_v + do_this_population.allelic_w + do_this_population.allelic_x + do_this_population.allelic_y + do_this_population.allelic_z;

        if(message_final_check_to_do && sum_allelic > 0)
        {
          cout << "doing final check: " << endl;
          message_final_check_to_do = false;
        }

        if(sum_allelic != sum_allelic_neutral)
        {
          no_failure = false;
          cout << " ++++++++++++++++++Something went wrong: issue 2 +++++++++++" << endl;
        }
        if(sum_allelic != do_this_population.pop_size)
        {
          no_failure = false;
          cout << " ++++++++++++++++++Something went wrong: issue 3 +++++++++++" << endl;
        }

        total_pop_size_from_individual_populations = total_pop_size_from_individual_populations + do_this_population.pop_size;
      }

      if(total_pop_size_from_individual_populations != do_this_species.total_pop_size)
      {
        no_failure = false;
        cout << " ++++++++++++++++++Something went wrong: issue 3 +++++++++++" << endl;
      }


    }

  }

  cout << "Checking spatial abundance" << endl;
  for (int i = 0; i < y_max; i++)
  {
    for (int j = 0; j < x_max; j++)
    {
      if(map_for_check[i][j].total_abundance_cell != map1[i][j].total_abundance_cell)
        cout << " ++++++++++++++++++Something went wrong: issue 4 +++++++++++" << endl;
    }
  }



  return no_failure;
}



void population_structure::happening_mutation()
{
  // decrease frequency of an allele by one and increase the frecuency of another one by one.
  // cout << "here at happening mutation from population structure " << endl;
  // cout << "allelic_a mutation " << allelic_a << endl;
  // cout << "allelic_b mutation " << allelic_b << endl;
  // cout << "allelic_c mutation " << allelic_c << endl;
  // cout << "allelic_d mutation " << allelic_d << endl;
  // cout << "allelic_e mutation " << allelic_e << endl;


  //  cout << endl;
  // cout << "______________after mutation: " << endl;
  vector <int> position_alleles;
  if(allelic_a >= 1){
    position_alleles.push_back(0);
  }
  if(allelic_b >= 1){
    position_alleles.push_back(1);
  }
  if(allelic_c >= 1){
    position_alleles.push_back(2);
  }
  if(allelic_d >= 1){
    position_alleles.push_back(3);
  }
  if(allelic_e >= 1){
    position_alleles.push_back(4);
  }

  //  cout << "here position_alleles.size()" << position_alleles.size() << endl;
  int which_allele_from;
  int position_allele_from;
  position_allele_from = give_me_random_uniform(0, (position_alleles.size()-1)); // the ones where we can take alleles from

  which_allele_from = position_alleles[position_allele_from];

  if(which_allele_from == 0)
  {
    allelic_a = allelic_a - 1;
  }
  if(which_allele_from == 1)
  {
    allelic_b = allelic_b - 1;
  }
  if(which_allele_from == 2)
  {
    allelic_c = allelic_c - 1;
  }
  if(which_allele_from == 3)
  {
    allelic_d = allelic_d - 1;
  }
  if(which_allele_from == 4)
  {
    allelic_e = allelic_e - 1;
  }

  vector <int> allelic_positions_all;
  allelic_positions_all.push_back(0);
  allelic_positions_all.push_back(1);
  allelic_positions_all.push_back(2);
  allelic_positions_all.push_back(3);
  allelic_positions_all.push_back(4);

  allelic_positions_all.erase(allelic_positions_all.begin() + which_allele_from);

  //  cout << "there position_alleles.size()" << position_alleles.size() << endl;
  int which_allele_to;
  int position_allele_to;
  position_allele_to = give_me_random_uniform(0, 3);

  which_allele_to = allelic_positions_all[position_allele_to];
  if(which_allele_to == 0)
  {
    allelic_a = allelic_a + 1;
  }
  if(which_allele_to == 1)
  {
    allelic_b = allelic_b + 1;
  }
  if(which_allele_to == 2)
  {
    allelic_c = allelic_c + 1;
  }
  if(which_allele_to == 3)
  {
    allelic_d = allelic_d + 1;
  }
  if(which_allele_to == 4)
  {
    allelic_e = allelic_e + 1;
  }
  //
  //   cout << "allelic_a mutation " << allelic_a << endl;
  //   cout << "allelic_b mutation " << allelic_b << endl;
  //   cout << "allelic_c mutation " << allelic_c << endl;
  //   cout << "allelic_d mutation " << allelic_d << endl;
  //   cout << "allelic_e mutation " << allelic_e << endl;

  // neutral now


  //cout << endl;
  // cout << "______________after mutation: " << endl;
  vector <int> position_alleles_neutral;
  if(allelic_v >= 1){
    position_alleles_neutral.push_back(0);
  }
  if(allelic_w >= 1){
    position_alleles_neutral.push_back(1);
  }
  if(allelic_x >= 1){
    position_alleles_neutral.push_back(2);
  }
  if(allelic_y >= 1){
    position_alleles_neutral.push_back(3);
  }
  if(allelic_z >= 1){
    position_alleles_neutral.push_back(4);
  }

  // cout << "here position_alleles.size()" << position_alleles_neutral.size() << endl;
  int which_allele_from_neutral;
  int position_allele_from_neutral;
  position_allele_from_neutral = give_me_random_uniform(0, (position_alleles_neutral.size()-1)); // the onces where we can take alleles from

  which_allele_from_neutral = position_alleles_neutral[position_allele_from_neutral];

  if(which_allele_from_neutral == 0)
  {
    allelic_v = allelic_v - 1;
  }
  if(which_allele_from_neutral == 1)
  {
    allelic_w = allelic_w - 1;
  }
  if(which_allele_from_neutral == 2)
  {
    allelic_x = allelic_x - 1;
  }
  if(which_allele_from_neutral == 3)
  {
    allelic_y = allelic_y - 1;
  }
  if(which_allele_from_neutral == 4)
  {
    allelic_z = allelic_z - 1;
  }

  vector <int> allelic_positions_all_neutral;
  allelic_positions_all_neutral.push_back(0);
  allelic_positions_all_neutral.push_back(1);
  allelic_positions_all_neutral.push_back(2);
  allelic_positions_all_neutral.push_back(3);
  allelic_positions_all_neutral.push_back(4);

  allelic_positions_all_neutral.erase(allelic_positions_all_neutral.begin() + which_allele_from_neutral);

  // cout << "there allelic_positions_all_neutral.size()" << allelic_positions_all_neutral.size() << endl;
  int which_allele_to_neutral;
  int position_allele_to_neutral;
  position_allele_to_neutral = give_me_random_uniform(0, 3);

  which_allele_to_neutral = allelic_positions_all_neutral[position_allele_to_neutral];
  if(which_allele_to_neutral == 0)
  {
    allelic_v = allelic_v + 1;
  }
  if(which_allele_to_neutral == 1)
  {
    allelic_w = allelic_w + 1;
  }
  if(which_allele_to_neutral == 2)
  {
    allelic_x = allelic_x + 1;
  }
  if(which_allele_to_neutral == 3)
  {
    allelic_y = allelic_y + 1;
  }
  if(which_allele_to_neutral == 4)
  {
    allelic_z = allelic_z + 1;
  }

  // cout << "allelic_v mutation " << allelic_v << endl;
  // cout << "allelic_w mutation " << allelic_w << endl;
  // cout << "allelic_x mutation " << allelic_x << endl;
  // cout << "allelic_y mutation " << allelic_y << endl;
  // cout << "allelic_z mutation " << allelic_z << endl;
}



void species::happening_mutation_this_species()
{
  //  cout << "here starting happening mutation " << endl;
  int random_population_to_mutate;
  random_population_to_mutate = give_me_random_uniform(1, presence.size());
  random_population_to_mutate = random_population_to_mutate - 1;
  //  cout << "random_number " << random_population_to_mutate << endl;
  populations_this_species[random_population_to_mutate].happening_mutation();
}

/*
 // gene flow with equal number of indivuals exchanged
 void species::happening_gene_flow(double percentage_flow)
 {
 population_structure allelic_pool_neighb;

 int focal_cell;


 focal_cell = give_me_random_uniform(0, (presence.size() - 1));

 vector<yx> neighbors_focal;
 neighbors_focal = find_neighbor(focal_cell);
 vector <int> id_neigh_cells;

 for(int ii = 0; ii < neighbors_focal.size(); ++ii)
 {
 for(int jj = 0; jj < presence.size(); ++jj)
 {
 if(neighbors_focal[ii].x == presence[jj].x && neighbors_focal[ii].y == presence[jj].y)
 {
 if(focal_cell != jj) // to exclude the focal cell itself
 {
 id_neigh_cells.push_back(jj);
 }

 }
 }
 }
 //  cout << "total id_neigh_cells " << id_neigh_cells.size() << endl;
 if(id_neigh_cells.size() > 0)// it is otherwise  an isolated pop
 {
 vector <int> weights_for_wallenius;
 weights_for_wallenius.push_back(1);
 weights_for_wallenius.push_back(1);
 weights_for_wallenius.push_back(1);
 weights_for_wallenius.push_back(1);
 weights_for_wallenius.push_back(1);

 int random_id_selector;
 random_id_selector = id_neigh_cells[give_me_random_uniform(0, (id_neigh_cells.size() - 1))];
 // cout << "focal_cell " << focal_cell << endl;
 // cout <<  "random_id_selector " << random_id_selector << endl;
 // cout << "populations_this_species[focal_cell].pop_size: " << populations_this_species[focal_cell].pop_size  << endl;
 // cout << "total pops this species here: " << populations_this_species.size() << endl;
 // cout << "populations_this_species[random_id_selector].pop_size: " << populations_this_species[random_id_selector].pop_size  << endl;
 // cout << "     pre-GENEflow" << endl;
 // cout<<"focal_cell.allelic_a: " << populations_this_species[focal_cell].allelic_a << endl;
 // cout<<"focal_cell.allelic_b: " << populations_this_species[focal_cell].allelic_b << endl;
 // cout<<"focal_cell.allelic_c: " << populations_this_species[focal_cell].allelic_c << endl;
 // cout<<"focal_cell.allelic_d: " << populations_this_species[focal_cell].allelic_d << endl;
 // cout<<"focal_cell.allelic_e: " << populations_this_species[focal_cell].allelic_e << endl;
 // cout << "           " << endl;
 // cout<<"pop_exchange_1.allelic_a: " << populations_this_species[random_id_selector].allelic_a << endl;
 // cout<<"pop_exchange_1.allelic_b: " << populations_this_species[random_id_selector].allelic_b << endl;
 // cout<<"pop_exchange_1.allelic_c: " << populations_this_species[random_id_selector].allelic_c << endl;
 // cout<<"pop_exchange_1.allelic_d: " << populations_this_species[random_id_selector].allelic_d << endl;
 // cout<<"pop_exchange_1.allelic_e: " << populations_this_species[random_id_selector].allelic_e << endl;

 double number_alleles_to_exchage;

 vector <double> percentage_per_each_pop;

 percentage_per_each_pop.push_back(populations_this_species[focal_cell].pop_size * (percentage_flow/100.0));
 percentage_per_each_pop.push_back(populations_this_species[random_id_selector].pop_size * (percentage_flow/100.0));

 // select the smaller amount out of both populations. If a population is much bigger than the other one, a 20% might
 // be more than the total population of the smaller population
 number_alleles_to_exchage = *min_element(percentage_per_each_pop.begin(),percentage_per_each_pop.end());
 number_alleles_to_exchage = floor(number_alleles_to_exchage);
 if(number_alleles_to_exchage == 0)

 {
 number_alleles_to_exchage = 1;
 }

 vector <int> contribution_allelic_frequency_focal;
 vector <int> current_allelic_frequency_focal;
 current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_a);
 current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_b);
 current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_c);
 current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_d);
 current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_e);
 contribution_allelic_frequency_focal = give_me_random_wallenius(weights_for_wallenius.size(),current_allelic_frequency_focal, weights_for_wallenius,number_alleles_to_exchage );

 vector <int> contribution_allelic_frequency_other;
 vector <int> current_allelic_frequency_other;
 current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_a);
 current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_b);
 current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_c);
 current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_d);
 current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_e);
 contribution_allelic_frequency_other = give_me_random_wallenius(weights_for_wallenius.size(),current_allelic_frequency_other, weights_for_wallenius,number_alleles_to_exchage );

 vector <int> pooled_alleles;
 pooled_alleles.push_back(contribution_allelic_frequency_focal[0] + contribution_allelic_frequency_other[0]);
 pooled_alleles.push_back(contribution_allelic_frequency_focal[1] + contribution_allelic_frequency_other[1]);
 pooled_alleles.push_back(contribution_allelic_frequency_focal[2] + contribution_allelic_frequency_other[2]);
 pooled_alleles.push_back(contribution_allelic_frequency_focal[3] + contribution_allelic_frequency_other[3]);
 pooled_alleles.push_back(contribution_allelic_frequency_focal[4] + contribution_allelic_frequency_other[4]);

 vector <int>  sampled_for_focal;

 sampled_for_focal = give_me_random_wallenius(weights_for_wallenius.size(),pooled_alleles, weights_for_wallenius,number_alleles_to_exchage );
 vector <int> sampled_for_other;
 sampled_for_other.push_back(pooled_alleles[0] - sampled_for_focal[0]);
 sampled_for_other.push_back(pooled_alleles[1] - sampled_for_focal[1]);
 sampled_for_other.push_back(pooled_alleles[2] - sampled_for_focal[2]);
 sampled_for_other.push_back(pooled_alleles[3] - sampled_for_focal[3]);
 sampled_for_other.push_back(pooled_alleles[4] - sampled_for_focal[4]);




 populations_this_species[focal_cell].allelic_a = populations_this_species[focal_cell].allelic_a - (contribution_allelic_frequency_focal[0] - sampled_for_focal[0]);
 populations_this_species[focal_cell].allelic_b = populations_this_species[focal_cell].allelic_b - (contribution_allelic_frequency_focal[1] - sampled_for_focal[1]);
 populations_this_species[focal_cell].allelic_c = populations_this_species[focal_cell].allelic_c - (contribution_allelic_frequency_focal[2] - sampled_for_focal[2]);
 populations_this_species[focal_cell].allelic_d = populations_this_species[focal_cell].allelic_d - (contribution_allelic_frequency_focal[3] - sampled_for_focal[3]);
 populations_this_species[focal_cell].allelic_e = populations_this_species[focal_cell].allelic_e - (contribution_allelic_frequency_focal[4] - sampled_for_focal[4]);

 populations_this_species[random_id_selector].allelic_a = populations_this_species[random_id_selector].allelic_a - (contribution_allelic_frequency_other[0] - sampled_for_other[0]);
 populations_this_species[random_id_selector].allelic_b = populations_this_species[random_id_selector].allelic_b - (contribution_allelic_frequency_other[1] - sampled_for_other[1]);
 populations_this_species[random_id_selector].allelic_c = populations_this_species[random_id_selector].allelic_c - (contribution_allelic_frequency_other[2] - sampled_for_other[2]);
 populations_this_species[random_id_selector].allelic_d = populations_this_species[random_id_selector].allelic_d - (contribution_allelic_frequency_other[3] - sampled_for_other[3]);
 populations_this_species[random_id_selector].allelic_e = populations_this_species[random_id_selector].allelic_e - (contribution_allelic_frequency_other[4] - sampled_for_other[4]);


 vector <int> contribution_allelic_frequency_focal_neutral;
 vector <int> current_allelic_frequency_focal_neutral;
 current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_v);
 current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_w);
 current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_x);
 current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_y);
 current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_z);
 contribution_allelic_frequency_focal_neutral = give_me_random_wallenius(weights_for_wallenius.size(),current_allelic_frequency_focal_neutral, weights_for_wallenius,number_alleles_to_exchage );

 vector <int> contribution_allelic_frequency_other_neutral;
 vector <int> current_allelic_frequency_other_neutral;
 current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_v);
 current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_w);
 current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_x);
 current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_y);
 current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_z);
 contribution_allelic_frequency_other_neutral = give_me_random_wallenius(weights_for_wallenius.size(),current_allelic_frequency_other_neutral, weights_for_wallenius,number_alleles_to_exchage );

 vector <int> pooled_alleles_neutral;
 pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[0] + contribution_allelic_frequency_other_neutral[0]);
 pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[1] + contribution_allelic_frequency_other_neutral[1]);
 pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[2] + contribution_allelic_frequency_other_neutral[2]);
 pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[3] + contribution_allelic_frequency_other_neutral[3]);
 pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[4] + contribution_allelic_frequency_other_neutral[4]);

 vector <int>  sampled_for_focal_neutral;

 sampled_for_focal_neutral = give_me_random_wallenius(weights_for_wallenius.size(),pooled_alleles_neutral, weights_for_wallenius,number_alleles_to_exchage );
 vector <int> sampled_for_other_neutral;
 sampled_for_other_neutral.push_back(pooled_alleles_neutral[0] - sampled_for_focal_neutral[0]);
 sampled_for_other_neutral.push_back(pooled_alleles_neutral[1] - sampled_for_focal_neutral[1]);
 sampled_for_other_neutral.push_back(pooled_alleles_neutral[2] - sampled_for_focal_neutral[2]);
 sampled_for_other_neutral.push_back(pooled_alleles_neutral[3] - sampled_for_focal_neutral[3]);
 sampled_for_other_neutral.push_back(pooled_alleles_neutral[4] - sampled_for_focal_neutral[4]);

 populations_this_species[focal_cell].allelic_v = populations_this_species[focal_cell].allelic_v - (contribution_allelic_frequency_focal_neutral[0] - sampled_for_focal_neutral[0]);
 populations_this_species[focal_cell].allelic_w = populations_this_species[focal_cell].allelic_w - (contribution_allelic_frequency_focal_neutral[1] - sampled_for_focal_neutral[1]);
 populations_this_species[focal_cell].allelic_x = populations_this_species[focal_cell].allelic_x - (contribution_allelic_frequency_focal_neutral[2] - sampled_for_focal_neutral[2]);
 populations_this_species[focal_cell].allelic_y = populations_this_species[focal_cell].allelic_y - (contribution_allelic_frequency_focal_neutral[3] - sampled_for_focal_neutral[3]);
 populations_this_species[focal_cell].allelic_z = populations_this_species[focal_cell].allelic_z - (contribution_allelic_frequency_focal_neutral[4] - sampled_for_focal_neutral[4]);

 populations_this_species[random_id_selector].allelic_v = populations_this_species[random_id_selector].allelic_v - (contribution_allelic_frequency_other_neutral[0] - sampled_for_other_neutral[0]);
 populations_this_species[random_id_selector].allelic_w = populations_this_species[random_id_selector].allelic_w - (contribution_allelic_frequency_other_neutral[1] - sampled_for_other_neutral[1]);
 populations_this_species[random_id_selector].allelic_x = populations_this_species[random_id_selector].allelic_x - (contribution_allelic_frequency_other_neutral[2] - sampled_for_other_neutral[2]);
 populations_this_species[random_id_selector].allelic_y = populations_this_species[random_id_selector].allelic_y - (contribution_allelic_frequency_other_neutral[3] - sampled_for_other_neutral[3]);
 populations_this_species[random_id_selector].allelic_z = populations_this_species[random_id_selector].allelic_z - (contribution_allelic_frequency_other_neutral[4] - sampled_for_other_neutral[4]);


 }
 else
 {
 cout <<" gene flow attempted but no population was adjacent to focal one" << endl;
 }



 } //end of it
 */




// happening_gene_flow with assymetrical exchange
void species::happening_gene_flow(double percentage_flow,landscape **map1)
{
  int focal_cell;
  //
  // cout << "species presence " << endl;
  //   for(int i = 0; i < presence.size(); ++i)
  //   {
  //   cout << presence[i].x << "_"<< presence[i].y << "map k: " << map1 [presence[i].y - 1][presence[i].x - 1].temperature << endl;
  //   }


  focal_cell = give_me_random_uniform(0, (presence.size() - 1) );

  //cout << "focal_cell is " <<  focal_cell << " : " <<presence[focal_cell].x << " : "<< presence[focal_cell].y << endl;
  vector<yx> neighbors_focal;
  neighbors_focal = find_neighbor(focal_cell);


  // cout << "neighbors_focal" << endl;
  // for(int i = 0; i < neighbors_focal.size(); ++i)
  // {
  //   cout << neighbors_focal[i].x << "_"<< neighbors_focal[i].y << endl;
  // }




  vector <int> id_neigh_cells;

  // find the cells occupied by the species that neighbor the focal cell
  for(int ii = 0; ii < presence.size(); ++ii)
  {
    for(int jj = 0; jj < neighbors_focal.size(); ++jj)
    {
      if(presence[ii].x == neighbors_focal[jj].x  && presence[ii].y == neighbors_focal[jj].y )
      {
        if( ii != focal_cell) // to exclude the focal cell itself
        {
          // cout << "id_this: " << ii << endl;
          id_neigh_cells.push_back(ii);
        }

      }
    }
  }
  //  cout << "total id_neigh_cells " << id_neigh_cells.size() << endl;
  if(id_neigh_cells.size() > 0)// it is otherwise  an isolated pop
  {
    vector <int> weights_for_wallenius;
    weights_for_wallenius.push_back(1);
    weights_for_wallenius.push_back(1);
    weights_for_wallenius.push_back(1);
    weights_for_wallenius.push_back(1);
    weights_for_wallenius.push_back(1);

    int random_id_selector;
    random_id_selector = id_neigh_cells[give_me_random_uniform(0, (id_neigh_cells.size() - 1))];
    // cout << "focal_cell " << focal_cell << endl;
    // cout <<  "random_id_selector " << random_id_selector << endl;
    // cout << "populations_this_species[focal_cell].pop_size: " << populations_this_species[focal_cell].pop_size  << endl;
    // cout << "total pops this species here: " << populations_this_species.size() << endl;
    // cout << "populations_this_species[random_id_selector].pop_size: " << populations_this_species[random_id_selector].pop_size  << endl;
    // cout << "     pre-GENEflow" << endl;
    // cout<<"focal_cell.allelic_a: " << populations_this_species[focal_cell].allelic_a << endl;
    // cout<<"focal_cell.allelic_b: " << populations_this_species[focal_cell].allelic_b << endl;
    // cout<<"focal_cell.allelic_c: " << populations_this_species[focal_cell].allelic_c << endl;
    // cout<<"focal_cell.allelic_d: " << populations_this_species[focal_cell].allelic_d << endl;
    // cout<<"focal_cell.allelic_e: " << populations_this_species[focal_cell].allelic_e << endl;
    // cout << "           " << endl;
    // cout<<"pop_exchange_1.allelic_a: " << populations_this_species[random_id_selector].allelic_a << endl;
    // cout<<"pop_exchange_1.allelic_b: " << populations_this_species[random_id_selector].allelic_b << endl;
    // cout<<"pop_exchange_1.allelic_c: " << populations_this_species[random_id_selector].allelic_c << endl;
    // cout<<"pop_exchange_1.allelic_d: " << populations_this_species[random_id_selector].allelic_d << endl;
    // cout<<"pop_exchange_1.allelic_e: " << populations_this_species[random_id_selector].allelic_e << endl;
    // cout << "random_id_selector: " << random_id_selector << endl;
    //     cout << "neigh to go: " << presence[random_id_selector].x << "_" <<presence[random_id_selector].y << endl;

    int number_alleles_to_giveaway_focal;
    int number_alleles_to_giveaway_neighbour;

    int k_at_focal;
    int k_at_neighbour;

    k_at_focal = map1[presence[focal_cell].y - 1 ][presence[focal_cell ].x - 1].k_patch;
    k_at_neighbour = map1[presence[random_id_selector].y - 1][presence[random_id_selector].x - 1].k_patch;
    int pre_genflow_cell_abundance_focal;
    int pre_genflow_cell_abundance_neighbour;

    pre_genflow_cell_abundance_focal = map1[presence[focal_cell ].y - 1][presence[focal_cell ].x - 1].total_abundance_cell;
    pre_genflow_cell_abundance_neighbour = map1[presence[random_id_selector].y - 1][presence[random_id_selector].x - 1].total_abundance_cell;
    // int t_at_focal;
    // int t_at_neighbour;
    // t_at_focal = map1[presence[focal_cell ].y - 1 ][presence[focal_cell ].x - 1].temperature;
    // t_at_neighbour = map1[presence[random_id_selector].y - 1][presence[random_id_selector].x - 1].temperature;

    // cout << "k_at_focal " << k_at_focal << "k_at_neighbour " << k_at_neighbour << endl;
    // cout << "t_at_focal " << t_at_focal << "t_at_neighbour " << t_at_neighbour << endl;

    if (k_at_focal == -9)
    {
      stop("invalid k value at  k_at_focal ");
    }
    if (k_at_neighbour == -9)
    {
      stop("invalid k value at  k_at_neighbour ");
    }

    // cout << " original populations_this_species[focal_cell].pop_size: " << populations_this_species[focal_cell].pop_size << endl;
    // cout << "original populations_this_species[random_id_selector].pop_size: " << populations_this_species[random_id_selector].pop_size << endl;
    // number_alleles_to_giveaway_focal = floor((populations_this_species[focal_cell].pop_size * (percentage_flow/100.0)));
    number_alleles_to_giveaway_neighbour = floor((populations_this_species[random_id_selector].pop_size * (percentage_flow/100.0)));

    if(number_alleles_to_giveaway_focal == 0)
    {
      number_alleles_to_giveaway_focal = 1;
    }
    if(number_alleles_to_giveaway_neighbour == 0)
    {
      number_alleles_to_giveaway_neighbour = 1;
    }
    // increase in the number of individuals cannot go beyond k


    // cout << "     considering k " << endl;
    //     cout << "populations_this_species[focal_cell].pop_size: " << populations_this_species[focal_cell].pop_size  << endl;
    //     cout << "populations_this_species[random_id_selector].pop_size: " << populations_this_species[random_id_selector].pop_size  << endl;
    //     cout << "drawn number_alleles_to_giveaway_focal: " << number_alleles_to_giveaway_focal << endl;
    //     cout << "drawn number_alleles_to_giveaway_neighbour: " << number_alleles_to_giveaway_neighbour << endl;
    int went_beyond = 0;
    if((pre_genflow_cell_abundance_focal - number_alleles_to_giveaway_focal  + number_alleles_to_giveaway_neighbour) > k_at_focal )
    {
      // cout << "k_at_neighbour " << k_at_neighbour << endl;
      // cout << "from if statement " << (pre_genflow_cell_abundance_focal - number_alleles_to_giveaway_focal  + number_alleles_to_giveaway_neighbour) << endl;
      // cout << "number_alleles_to_giveaway_focal " <<  number_alleles_to_giveaway_focal <<endl;
      // cout << "number_alleles_to_giveaway_neighbour " << number_alleles_to_giveaway_neighbour << endl;

      went_beyond = (pre_genflow_cell_abundance_focal - number_alleles_to_giveaway_focal + number_alleles_to_giveaway_neighbour)  - k_at_focal;
      //cout << "went_beyond1 : " << went_beyond << endl;
      number_alleles_to_giveaway_neighbour =  number_alleles_to_giveaway_neighbour - went_beyond;

    }

    if((pre_genflow_cell_abundance_neighbour - number_alleles_to_giveaway_neighbour  + number_alleles_to_giveaway_focal) > k_at_neighbour )
    {
      // cout << "k_at_neighbour " << k_at_neighbour << endl;
      // cout << "from if statement " << (pre_genflow_cell_abundance_neighbour - number_alleles_to_giveaway_neighbour  + number_alleles_to_giveaway_focal) << endl;
      // cout << "number_alleles_to_giveaway_focal " <<  number_alleles_to_giveaway_focal <<endl;
      // cout << "number_alleles_to_giveaway_neighbour " << number_alleles_to_giveaway_neighbour << endl;

      went_beyond = (pre_genflow_cell_abundance_neighbour - number_alleles_to_giveaway_neighbour  + number_alleles_to_giveaway_focal)  - k_at_neighbour;
      //cout << "went_beyond2 : " << went_beyond << endl;
      number_alleles_to_giveaway_focal =  number_alleles_to_giveaway_focal - went_beyond;
    }
    // cout << "final number_alleles_to_giveaway_focal: " << number_alleles_to_giveaway_focal << endl;
    // cout << "final number_alleles_to_giveaway_neighbour: " << number_alleles_to_giveaway_neighbour << endl;
    if(went_beyond != 0){
      // cout << "this is how the problem was solved " << "went_beyond1 " << number_alleles_to_giveaway_neighbour << "-went_beyond1 " <<number_alleles_to_giveaway_focal << endl;
      // cout << populations_this_species[focal_cell].pop_size - number_alleles_to_giveaway_focal  + number_alleles_to_giveaway_neighbour  << endl;
      //
      // cout << populations_this_species[random_id_selector].pop_size - number_alleles_to_giveaway_neighbour  + number_alleles_to_giveaway_focal << endl;
      //
      // cout << "number_alleles_to_giveaway_neighbour " << number_alleles_to_giveaway_neighbour << endl;
      // cout << "number_alleles_to_giveaway_focal " << number_alleles_to_giveaway_focal << endl;
      // cout << "focal " << map1[presence[focal_cell].y - 1 ][presence[focal_cell ].x - 1].total_abundance_cell << endl;
      // cout << "neighgh " << map1[presence[random_id_selector].y - 1][presence[random_id_selector].x - 1].total_abundance_cell << endl;
      //

      //stop("it went beyond so I stop");
    }

    // Neutral part

    //  cout<<"focal_cell.allelic_v: " << populations_this_species[focal_cell].allelic_v << endl;
    // cout<<"focal_cell.allelic_w: " << populations_this_species[focal_cell].allelic_w << endl;
    // cout<<"focal_cell.allelic_x: " << populations_this_species[focal_cell].allelic_x << endl;
    // cout<<"focal_cell.allelic_y: " << populations_this_species[focal_cell].allelic_y << endl;
    // cout<<"focal_cell.allelic_z: " << populations_this_species[focal_cell].allelic_z << endl;
    // cout << "           " << endl;
    //
    // cout<<"pop_exchange_1.allelic_v: " << populations_this_species[random_id_selector].allelic_v << endl;
    // cout<<"pop_exchange_1.allelic_w: " << populations_this_species[random_id_selector].allelic_w << endl;
    // cout<<"pop_exchange_1.allelic_x: " << populations_this_species[random_id_selector].allelic_x << endl;
    // cout<<"pop_exchange_1.allelic_y: " << populations_this_species[random_id_selector].allelic_y << endl;
    // cout<<"pop_exchange_1.allelic_z: " << populations_this_species[random_id_selector].allelic_z << endl;


    vector <int> contribution_allelic_frequency_focal;
    vector <int> current_allelic_frequency_focal;
    current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_a);
    current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_b);
    current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_c);
    current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_d);
    current_allelic_frequency_focal.push_back(populations_this_species[focal_cell].allelic_e);
    contribution_allelic_frequency_focal = give_me_random_wallenius(weights_for_wallenius.size(),current_allelic_frequency_focal, weights_for_wallenius,number_alleles_to_giveaway_focal );

    vector <int> contribution_allelic_frequency_other;
    vector <int> current_allelic_frequency_other;
    current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_a);
    current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_b);
    current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_c);
    current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_d);
    current_allelic_frequency_other.push_back(populations_this_species[random_id_selector].allelic_e);
    contribution_allelic_frequency_other = give_me_random_wallenius(weights_for_wallenius.size(),current_allelic_frequency_other, weights_for_wallenius,number_alleles_to_giveaway_neighbour );

    vector <int> pooled_alleles;
    pooled_alleles.push_back(contribution_allelic_frequency_focal[0] + contribution_allelic_frequency_other[0]);
    pooled_alleles.push_back(contribution_allelic_frequency_focal[1] + contribution_allelic_frequency_other[1]);
    pooled_alleles.push_back(contribution_allelic_frequency_focal[2] + contribution_allelic_frequency_other[2]);
    pooled_alleles.push_back(contribution_allelic_frequency_focal[3] + contribution_allelic_frequency_other[3]);
    pooled_alleles.push_back(contribution_allelic_frequency_focal[4] + contribution_allelic_frequency_other[4]);

    vector <int>  sampled_for_focal;

    sampled_for_focal = give_me_random_wallenius(weights_for_wallenius.size(),pooled_alleles, weights_for_wallenius,number_alleles_to_giveaway_neighbour );
    vector <int> sampled_for_other;
    sampled_for_other.push_back(pooled_alleles[0] - sampled_for_focal[0]);
    sampled_for_other.push_back(pooled_alleles[1] - sampled_for_focal[1]);
    sampled_for_other.push_back(pooled_alleles[2] - sampled_for_focal[2]);
    sampled_for_other.push_back(pooled_alleles[3] - sampled_for_focal[3]);
    sampled_for_other.push_back(pooled_alleles[4] - sampled_for_focal[4]);

    // cout << "pooled_alleles[0] " <<pooled_alleles[0] << endl;
    // cout << "sampled_for_focal[0] " <<sampled_for_focal[0] << endl;
    // cout << "sampled_for_other[0] " <<sampled_for_other[0] << endl;


    populations_this_species[focal_cell].allelic_a = populations_this_species[focal_cell].allelic_a - (contribution_allelic_frequency_focal[0] - sampled_for_focal[0]);
    populations_this_species[focal_cell].allelic_b = populations_this_species[focal_cell].allelic_b - (contribution_allelic_frequency_focal[1] - sampled_for_focal[1]);
    populations_this_species[focal_cell].allelic_c = populations_this_species[focal_cell].allelic_c - (contribution_allelic_frequency_focal[2] - sampled_for_focal[2]);
    populations_this_species[focal_cell].allelic_d = populations_this_species[focal_cell].allelic_d - (contribution_allelic_frequency_focal[3] - sampled_for_focal[3]);
    populations_this_species[focal_cell].allelic_e = populations_this_species[focal_cell].allelic_e - (contribution_allelic_frequency_focal[4] - sampled_for_focal[4]);

    populations_this_species[random_id_selector].allelic_a = populations_this_species[random_id_selector].allelic_a - (contribution_allelic_frequency_other[0] - sampled_for_other[0]);
    populations_this_species[random_id_selector].allelic_b = populations_this_species[random_id_selector].allelic_b - (contribution_allelic_frequency_other[1] - sampled_for_other[1]);
    populations_this_species[random_id_selector].allelic_c = populations_this_species[random_id_selector].allelic_c - (contribution_allelic_frequency_other[2] - sampled_for_other[2]);
    populations_this_species[random_id_selector].allelic_d = populations_this_species[random_id_selector].allelic_d - (contribution_allelic_frequency_other[3] - sampled_for_other[3]);
    populations_this_species[random_id_selector].allelic_e = populations_this_species[random_id_selector].allelic_e - (contribution_allelic_frequency_other[4] - sampled_for_other[4]);

    // cout << "populations_this_species[focal_cell].pop_size: " << populations_this_species[focal_cell].pop_size  << endl;
    // cout << "populations_this_species[random_id_selector].pop_size: " << populations_this_species[random_id_selector].pop_size  << endl;
    // cout << "real number_alleles_to_exchage: " << number_alleles_to_exchage << endl;
    // cout << "     post-GENEflow" << endl;
    // cout<<"focal_cell.allelic_a: " << populations_this_species[focal_cell].allelic_a << endl;
    // cout<<"focal_cell.allelic_b: " << populations_this_species[focal_cell].allelic_b << endl;
    // cout<<"focal_cell.allelic_c: " << populations_this_species[focal_cell].allelic_c << endl;
    // cout<<"focal_cell.allelic_d: " << populations_this_species[focal_cell].allelic_d << endl;
    // cout<<"focal_cell.allelic_e: " << populations_this_species[focal_cell].allelic_e << endl;
    // cout << "           " << endl;
    // cout<<"pop_exchange_1.allelic_a: " << populations_this_species[random_id_selector].allelic_a << endl;
    // cout<<"pop_exchange_1.allelic_b: " << populations_this_species[random_id_selector].allelic_b << endl;
    // cout<<"pop_exchange_1.allelic_c: " << populations_this_species[random_id_selector].allelic_c << endl;
    // cout<<"pop_exchange_1.allelic_d: " << populations_this_species[random_id_selector].allelic_d << endl;
    // cout<<"pop_exchange_1.allelic_e: " << populations_this_species[random_id_selector].allelic_e << endl;


    vector <int> contribution_allelic_frequency_focal_neutral;
    vector <int> current_allelic_frequency_focal_neutral;
    current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_v);
    current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_w);
    current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_x);
    current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_y);
    current_allelic_frequency_focal_neutral.push_back(populations_this_species[focal_cell].allelic_z);
    contribution_allelic_frequency_focal_neutral = give_me_random_wallenius(weights_for_wallenius.size(),current_allelic_frequency_focal_neutral, weights_for_wallenius,number_alleles_to_giveaway_focal );

    vector <int> contribution_allelic_frequency_other_neutral;
    vector <int> current_allelic_frequency_other_neutral;
    current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_v);
    current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_w);
    current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_x);
    current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_y);
    current_allelic_frequency_other_neutral.push_back(populations_this_species[random_id_selector].allelic_z);
    contribution_allelic_frequency_other_neutral = give_me_random_wallenius(weights_for_wallenius.size(),current_allelic_frequency_other_neutral, weights_for_wallenius,number_alleles_to_giveaway_neighbour );

    vector <int> pooled_alleles_neutral;
    pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[0] + contribution_allelic_frequency_other_neutral[0]);
    pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[1] + contribution_allelic_frequency_other_neutral[1]);
    pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[2] + contribution_allelic_frequency_other_neutral[2]);
    pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[3] + contribution_allelic_frequency_other_neutral[3]);
    pooled_alleles_neutral.push_back(contribution_allelic_frequency_focal_neutral[4] + contribution_allelic_frequency_other_neutral[4]);

    vector <int>  sampled_for_focal_neutral;

    sampled_for_focal_neutral = give_me_random_wallenius(weights_for_wallenius.size(),pooled_alleles_neutral, weights_for_wallenius,number_alleles_to_giveaway_neighbour );
    vector <int> sampled_for_other_neutral;
    sampled_for_other_neutral.push_back(pooled_alleles_neutral[0] - sampled_for_focal_neutral[0]);
    sampled_for_other_neutral.push_back(pooled_alleles_neutral[1] - sampled_for_focal_neutral[1]);
    sampled_for_other_neutral.push_back(pooled_alleles_neutral[2] - sampled_for_focal_neutral[2]);
    sampled_for_other_neutral.push_back(pooled_alleles_neutral[3] - sampled_for_focal_neutral[3]);
    sampled_for_other_neutral.push_back(pooled_alleles_neutral[4] - sampled_for_focal_neutral[4]);

    // cout << "pooled_alleles[0] " <<pooled_alleles_neutral[0] << endl;
    // cout << "sampled_for_focal[0] " <<sampled_for_focal_neutral[0] << endl;
    // cout << "sampled_for_other_neutral[0] " <<sampled_for_focal_neutral[0] << endl;


    populations_this_species[focal_cell].allelic_v = populations_this_species[focal_cell].allelic_v - (contribution_allelic_frequency_focal_neutral[0] - sampled_for_focal_neutral[0]);
    populations_this_species[focal_cell].allelic_w = populations_this_species[focal_cell].allelic_w - (contribution_allelic_frequency_focal_neutral[1] - sampled_for_focal_neutral[1]);
    populations_this_species[focal_cell].allelic_x = populations_this_species[focal_cell].allelic_x - (contribution_allelic_frequency_focal_neutral[2] - sampled_for_focal_neutral[2]);
    populations_this_species[focal_cell].allelic_y = populations_this_species[focal_cell].allelic_y - (contribution_allelic_frequency_focal_neutral[3] - sampled_for_focal_neutral[3]);
    populations_this_species[focal_cell].allelic_z = populations_this_species[focal_cell].allelic_z - (contribution_allelic_frequency_focal_neutral[4] - sampled_for_focal_neutral[4]);

    populations_this_species[random_id_selector].allelic_v = populations_this_species[random_id_selector].allelic_v - (contribution_allelic_frequency_other_neutral[0] - sampled_for_other_neutral[0]);
    populations_this_species[random_id_selector].allelic_w = populations_this_species[random_id_selector].allelic_w - (contribution_allelic_frequency_other_neutral[1] - sampled_for_other_neutral[1]);
    populations_this_species[random_id_selector].allelic_x = populations_this_species[random_id_selector].allelic_x - (contribution_allelic_frequency_other_neutral[2] - sampled_for_other_neutral[2]);
    populations_this_species[random_id_selector].allelic_y = populations_this_species[random_id_selector].allelic_y - (contribution_allelic_frequency_other_neutral[3] - sampled_for_other_neutral[3]);
    populations_this_species[random_id_selector].allelic_z = populations_this_species[random_id_selector].allelic_z - (contribution_allelic_frequency_other_neutral[4] - sampled_for_other_neutral[4]);



    // cout << "     post-GENEflow _____________NEUTRAL" << endl;
    // cout<<"focal_cell.allelic_v: " << populations_this_species[focal_cell].allelic_v << endl;
    // cout<<"focal_cell.allelic_w: " << populations_this_species[focal_cell].allelic_w << endl;
    // cout<<"focal_cell.allelic_x: " << populations_this_species[focal_cell].allelic_x << endl;
    // cout<<"focal_cell.allelic_y: " << populations_this_species[focal_cell].allelic_y << endl;
    // cout<<"focal_cell.allelic_z: " << populations_this_species[focal_cell].allelic_z << endl;
    // cout << "           " << endl;
    //
    // cout<<"pop_exchange_1.allelic_v: " << populations_this_species[random_id_selector].allelic_v << endl;
    // cout<<"pop_exchange_1.allelic_w: " << populations_this_species[random_id_selector].allelic_w << endl;
    // cout<<"pop_exchange_1.allelic_x: " << populations_this_species[random_id_selector].allelic_x << endl;
    // cout<<"pop_exchange_1.allelic_y: " << populations_this_species[random_id_selector].allelic_y << endl;
    // cout<<"pop_exchange_1.allelic_z: " << populations_this_species[random_id_selector].allelic_z << endl;

    populations_this_species[focal_cell].pop_size = populations_this_species[focal_cell].pop_size -  number_alleles_to_giveaway_focal  + number_alleles_to_giveaway_neighbour;
    populations_this_species[random_id_selector].pop_size = populations_this_species[random_id_selector].pop_size - number_alleles_to_giveaway_neighbour + number_alleles_to_giveaway_focal;
    // cout << "populations_this_species[focal_cell].pop_size: " << populations_this_species[focal_cell].pop_size  << endl;
    // cout << "populations_this_species[random_id_selector].pop_size: " << populations_this_species[random_id_selector].pop_size  << endl;


    //
    //     cout << "pre_genflow_cell_abundance_focal " << pre_genflow_cell_abundance_focal << endl;
    //     cout << "pre_genflow_cell_abundance_neighbour " << pre_genflow_cell_abundance_neighbour << endl;
    //
    //     cout << "number_alleles_to_giveaway_focal " << number_alleles_to_giveaway_focal << endl;
    //     cout << "number_alleles_to_giveaway_neighbour " << number_alleles_to_giveaway_neighbour << endl;


    map1[presence[focal_cell  ].y - 1][presence[focal_cell ].x - 1].total_abundance_cell =  pre_genflow_cell_abundance_focal - number_alleles_to_giveaway_focal + number_alleles_to_giveaway_neighbour ;
    map1[presence[random_id_selector].y - 1][presence[random_id_selector].x - 1].total_abundance_cell =  pre_genflow_cell_abundance_neighbour - number_alleles_to_giveaway_neighbour + number_alleles_to_giveaway_focal;

    //  cout <<     "map1[presence[focal_cell  ].y - 1][presence[focal_cell ].x - 1].total_abundance_cell   "  <<     map1[presence[focal_cell  ].y - 1][presence[focal_cell ].x - 1].total_abundance_cell << endl;
    // cout << "map1[presence[random_id_selector].y - 1][presence[random_id_selector].x - 1].total_abundance_cell  " <<  map1[presence[random_id_selector].y - 1][presence[random_id_selector].x - 1].total_abundance_cell << endl;
    //  cout << "cordinates map for focal and neighbour" << endl;
    //  cout <<     "focal   "  <<     presence[focal_cell  ].y - 1 << "y: " << presence[focal_cell ].x - 1 << endl;
    //  cout <<     "random_id_selector   "  <<     presence[random_id_selector  ].y - 1 << "y: " << presence[random_id_selector ].x - 1 << endl;

    succesful_geneflow_events = succesful_geneflow_events + 1;
  }
  else
  {
    //cout <<" gene flow attempted but no population was adjacent to focal one" << endl;
  }
}

void show_all_species_data(vector <species> all_species)
{
  cout << "richness in show_all_species_data function: " << all_species.size() << endl;
  for(int ii = 0; ii < all_species.size(); ++ii)
  {
    species do_this_species;
    do_this_species = all_species[ii];
    cout << "id: " << do_this_species.id << endl;
    cout << "range: " << do_this_species.range<< endl;
    cout << "total pop: " << do_this_species.total_pop_size<< endl;
    for(int ij = 0; ij < do_this_species.presence.size(); ++ij)
    {
      population_structure do_this_population;
      do_this_population = do_this_species.populations_this_species[ij];
      cout<< "                  pop size this population: " << do_this_population.pop_size<< endl;
      cout<< "                  fitness like: " << do_this_species.computed_rate_based_on_temperature[ij]<< endl;
      cout<< "            x y coordinate: " << do_this_species.presence[ij].x<<" , " << do_this_species.presence[ij].y << endl;
      cout<< "      allele a: " << do_this_population.allelic_a
          << " allele b: " << do_this_population.allelic_b
          << " allele c: " << do_this_population.allelic_c
          << " allele d: " << do_this_population.allelic_d
          << " allele e: " << do_this_population.allelic_e<< endl;

      cout<< "      allele v: " << do_this_population.allelic_v
          << " allele w: " << do_this_population.allelic_w
          << " allele x: " << do_this_population.allelic_x
          << " allele y: " << do_this_population.allelic_y
          << " allele z: " << do_this_population.allelic_z<< endl;
    }

  }

}
void species::happening_population_popchange_this_species(double sd_normal_distribution, landscape **map1, vector<int> alleles_adaptation_coef)
{

  int random_population_to_popchange;
  random_population_to_popchange = give_me_random_uniform(1, presence.size());
  yx this_cell;
  int original_size_this_pop;
  this_cell = presence[random_population_to_popchange - 1];
  original_size_this_pop = populations_this_species[random_population_to_popchange - 1].pop_size;
  populations_this_species[random_population_to_popchange - 1].happening_population_popchange(map1, sd_normal_distribution, alleles_adaptation_coef, this_cell); // -1 as it is index
  total_pop_size = total_pop_size + (populations_this_species[random_population_to_popchange - 1].pop_size - original_size_this_pop);
  computed_rate_based_on_temperature[random_population_to_popchange - 1] = link_fitnesslike_mu_gamma(this_cell, populations_this_species[random_population_to_popchange - 1], alleles_adaptation_coef, map1);
}

void population_structure::happening_population_popchange(landscape **map1,double sd_normal_distribution, vector<int> alleles_adaptation_coef, yx this_cell)
{
  int cell_temperature;
  cell_temperature = map1[this_cell.y - 1][this_cell.x - 1].temperature;
  int total_abundance_cell;
  int k_this_cell;
  total_abundance_cell = map1[this_cell.y - 1][this_cell.x - 1].total_abundance_cell;
  k_this_cell = map1[this_cell.y - 1][this_cell.x - 1].k_patch;
  int change_in_population;
  int new_pop_size;
  change_in_population = round(give_me_random_normal(0, sd_normal_distribution));
  if(change_in_population == 0){
    int random_value_to_decide;
    random_value_to_decide = give_me_random_normal(0,1);
    if(random_value_to_decide == 0)
    {
      change_in_population = 1;
    }
    else
    {
      change_in_population = -1;
    }

  }
  //change_in_population = 5; // to DELETE


  // if the population grows, it cannot exceed K
  if((change_in_population + total_abundance_cell) > k_this_cell)
  {
    change_in_population = k_this_cell - total_abundance_cell;
  }

  new_pop_size = pop_size + change_in_population ;
  // cout << "pop_size before "<< pop_size << endl;
  // cout << "change_in_population: " << change_in_population << endl;

  if(new_pop_size <= 0) // I do not let extirpation happening by the means of population change
  {
    //cout << "change_in_population stopped as it tried to delete this population" << endl;
  }
  else
  {
    vector<int> weights_for_wallenius;
    vector<int> weights_for_wallenius_neutral;

    for (int ii = 0; ii < alleles_adaptation_coef.size(); ++ii)
    {
      weights_for_wallenius.push_back(abs(alleles_adaptation_coef[ii] - cell_temperature));
      weights_for_wallenius_neutral.push_back(1);
    }
    if (change_in_population < 0) // population shrinks
    {
      //cout << "population shrinks" << endl;
      vector<int> current_allelic_frequency;
      current_allelic_frequency.push_back(allelic_a);
      current_allelic_frequency.push_back(allelic_b);
      current_allelic_frequency.push_back(allelic_c);
      current_allelic_frequency.push_back(allelic_d);
      current_allelic_frequency.push_back(allelic_e);

      vector<int> current_allelic_frequency_neutral;
      current_allelic_frequency_neutral.push_back(allelic_v);
      current_allelic_frequency_neutral.push_back(allelic_w);
      current_allelic_frequency_neutral.push_back(allelic_x);
      current_allelic_frequency_neutral.push_back(allelic_y);
      current_allelic_frequency_neutral.push_back(allelic_z);

      int number_items_sample;                         //
      number_items_sample = abs(change_in_population); //
      vector <int> change_allelic_frequency;
      vector <int> change_allelic_frequency_neutral;
      int number_alelles;
      number_alelles = current_allelic_frequency.size();
      // cout << "allelic_freq: "<<current_allelic_frequency[0] << " "<< current_allelic_frequency[1] << " "<< current_allelic_frequency[2] << " "<<
      //   current_allelic_frequency[3] <<" "<< current_allelic_frequency[4]<< " "<< endl;
      // cout << "number_items_sample " << number_items_sample << endl;
      //
      // cout << "allelic_freqNeutral: "<<current_allelic_frequency_neutral[0] << " "<< current_allelic_frequency_neutral[1] << " "<< current_allelic_frequency_neutral[2] << " "<<
      //   current_allelic_frequency_neutral[3] <<" "<< current_allelic_frequency_neutral[4]<< " "<< endl;
      // we will remove those individuals that are less well adapted, so weights_for_wallenius works fine
      change_allelic_frequency = give_me_random_wallenius(number_alelles,current_allelic_frequency, weights_for_wallenius, number_items_sample);
      //  cout << "now NEUTRAL " << endl;
      change_allelic_frequency_neutral = give_me_random_wallenius(number_alelles,current_allelic_frequency_neutral, weights_for_wallenius_neutral, number_items_sample);


      //   cout << "now NEUTRAL " << endl;
      //  cout << "allelic_freq"<<current_allelic_frequency_neutral[0] <<current_allelic_frequency_neutral[1] << current_allelic_frequency_neutral[2] << current_allelic_frequency_neutral[3] << current_allelic_frequency_neutral[4] << endl;


      // for (int ii = 0; ii < change_allelic_frequency.size(); ++ii)
      // {
      //   cout << "show this sampling from wallenius " << change_allelic_frequency[ii] << endl;
      // }


      allelic_a = allelic_a - change_allelic_frequency[0];
      allelic_b = allelic_b - change_allelic_frequency[1];
      allelic_c = allelic_c - change_allelic_frequency[2];
      allelic_d = allelic_d - change_allelic_frequency[3];
      allelic_e = allelic_e - change_allelic_frequency[4];

      allelic_v = allelic_v - change_allelic_frequency_neutral[0];
      allelic_w = allelic_w - change_allelic_frequency_neutral[1];
      allelic_x = allelic_x - change_allelic_frequency_neutral[2];
      allelic_y = allelic_y - change_allelic_frequency_neutral[3];
      allelic_z = allelic_z - change_allelic_frequency_neutral[4];


    }
    else // population grows
    {
      //cout << "population grows" << endl;

      // I will create new individuals, and need to know what alleles they will carry. This is function of the cell temperature. Individuals carrying the fittest
      // allele will produce more offspring. I need to define a gigantic pool of individuals where the new ones will be taking from based on the weights. But the
      // weights_for_wallenius needs to be treated differently, otherwise the sampling will only take be biased to individuals with higher weight i.e., far from temp.
      // If an allele does not exist at the
      // population, the new individuals cannot have such allele.

      vector <int> weights_for_wallenius_more_adapted;
      for(int ii = 0; ii < weights_for_wallenius.size(); ++ii)
      {
        if(weights_for_wallenius[ii] == 0){ // perfectly adapted
          weights_for_wallenius_more_adapted.push_back(1.0/0.5);
        } else{
          weights_for_wallenius_more_adapted.push_back(round((1.0/weights_for_wallenius[ii]) * 100));
        }

        // cout << "this weights: " << round((1.0/weights_for_wallenius[ii]) * 100) << endl;
      }

      vector<int> current_allelic_frequency;
      if(allelic_a > 0){
        current_allelic_frequency.push_back(1900); // any large number will do
      }
      else
      {
        current_allelic_frequency.push_back(0);
      }

      if(allelic_b > 0){
        current_allelic_frequency.push_back(1800); // any large number will do
      }
      else
      {
        current_allelic_frequency.push_back(0);
      }

      if(allelic_c > 0){
        current_allelic_frequency.push_back(1800); // any large number will do
      }
      else
      {
        current_allelic_frequency.push_back(0);
      }
      if(allelic_d > 0){
        current_allelic_frequency.push_back(1800); // any large number will do
      }
      else
      {
        current_allelic_frequency.push_back(0);
      }
      if(allelic_e > 0){
        current_allelic_frequency.push_back(1800); // any large number will do
      }
      else
      {
        current_allelic_frequency.push_back(0);
      }

      vector<int> current_allelic_frequency_neutral;
      if(allelic_v > 0){
        current_allelic_frequency_neutral.push_back(1800);
      }
      else
      {
        current_allelic_frequency_neutral.push_back(0);
      }
      if(allelic_w > 0){
        current_allelic_frequency_neutral.push_back(1800);
      }
      else
      {
        current_allelic_frequency_neutral.push_back(0);
      }
      if(allelic_x > 0){
        current_allelic_frequency_neutral.push_back(1800);
      }
      else
      {
        current_allelic_frequency_neutral.push_back(0);
      }
      if(allelic_y > 0){
        current_allelic_frequency_neutral.push_back(1800);
      }
      else
      {
        current_allelic_frequency_neutral.push_back(0);
      }
      if(allelic_z > 0){
        current_allelic_frequency_neutral.push_back(1800);
      }
      else
      {
        current_allelic_frequency_neutral.push_back(0);
      }
      int number_items_sample;                         //
      number_items_sample = abs(change_in_population);
      int number_alelles;
      number_alelles = current_allelic_frequency.size();
      vector <int> change_allelic_frequency;
      vector <int> change_allelic_frequency_neutral;


      change_allelic_frequency = give_me_random_wallenius(number_alelles,current_allelic_frequency, weights_for_wallenius_more_adapted, number_items_sample);
      // cout << "now NEUTRAL " << endl;
      //  cout << "allelic_freq: "<<current_allelic_frequency_neutral[0] << " "<< current_allelic_frequency_neutral[1] << " "<< current_allelic_frequency_neutral[2] << " "<<
      //    current_allelic_frequency_neutral[3] <<" "<< current_allelic_frequency_neutral[4]<< " "<< endl;
      change_allelic_frequency_neutral = give_me_random_wallenius(number_alelles,current_allelic_frequency_neutral, weights_for_wallenius_neutral, number_items_sample);


      allelic_a = allelic_a + change_allelic_frequency[0];
      allelic_b = allelic_b + change_allelic_frequency[1];
      allelic_c = allelic_c + change_allelic_frequency[2];
      allelic_d = allelic_d + change_allelic_frequency[3];
      allelic_e = allelic_e + change_allelic_frequency[4];

      allelic_v = allelic_v + change_allelic_frequency_neutral[0];
      allelic_w = allelic_w + change_allelic_frequency_neutral[1];
      allelic_x = allelic_x + change_allelic_frequency_neutral[2];
      allelic_y = allelic_y + change_allelic_frequency_neutral[3];
      allelic_z = allelic_z + change_allelic_frequency_neutral[4];
    }

    // cout << "allelic_a after:  "<< allelic_a << endl;
    // cout << "allelic_b after:  "<< allelic_b << endl;
    // cout << "allelic_c after:  "<< allelic_c << endl;
    // cout << "allelic_d after:  "<< allelic_d << endl;
    // cout << "allelic_e after:  "<< allelic_e << endl;
    //
    // cout << "new_pop_size :  " << new_pop_size << endl;


    // update map
    int current_pop_map;
    current_pop_map = map1[(this_cell.y - 1)][(this_cell.x - 1)].total_abundance_cell;
    map1[(this_cell.y - 1)][(this_cell.x - 1)].total_abundance_cell = current_pop_map - (pop_size - (new_pop_size));
    if(map1[(this_cell.y - 1)][(this_cell.x - 1)].total_abundance_cell == 0)
    {
      stop("this thing came down to zero");
    }
    pop_size = new_pop_size;
    // cout << "pop_size from species after: " << pop_size <<endl;

  }
  //new_pop_size = pop_size + abs(new_pop_size - pop_size); // to DELETE
}

void species::classify_elevation(bool add, landscape **map1)
{

  int elevation_last_event;
  elevation_last_event = map1[y_coordinate_last_event - 1][x_coordinate_last_event - 1].elevation;

  if (elevation_last_event == 4)
  {
    // cout <<  "  high:  " << x_coordinate_last_event << endl;
    if (add)
    {
      range_highlands = range_highlands + 1;
    }
    else
    {
      range_highlands = range_highlands - 1;
    }
  }
  if (elevation_last_event == 3)
  {
    // cout <<  "  high:  " << x_coordinate_last_event << endl;
    if (add)
    {
      range_intermediate1 = range_intermediate1 + 1;
    }
    else
    {
      range_intermediate1 = range_intermediate1 - 1;
    }
  }
  if (elevation_last_event == 2)
  {
    // cout <<  "  high:  " << x_coordinate_last_event << endl;
    if (add)
    {
      range_intermediate2 = range_intermediate2 + 1;
    }
    else
    {
      range_intermediate2 = range_intermediate2 - 1;
    }
  }
  if (elevation_last_event == 1)
  {
    // cout <<  "  high:  " << x_coordinate_last_event << endl;
    if (add)
    {
      range_lowlands = range_lowlands + 1;
    }
    else
    {
      range_lowlands = range_lowlands - 1;
    }
  }
}

void species::classify_elevation_origin(bool add)
{
  if (range_highlands == 1)
  {
    elevation_origin = 4;
  }

  if (range_intermediate1 == 1)
  {
    elevation_origin = 3;
  }

  if (range_intermediate2 == 1)
  {
    elevation_origin = 2;
  }

  if (range_lowlands == 1)
  {
    elevation_origin = 1;
  }
}

double link_fitnesslike_mu_gamma(yx this_location, population_structure this_pop, vector<int> alleles_adaptation_coef, landscape **map1)

{
  double unfitnesslike_thispopulation;
  double fitnesslike_thispopulation;
  int cell_temperature;
  cell_temperature = map1[this_location.y - 1][this_location.x - 1].temperature;
  vector <int> fitnesslike_per_allele;
  for (int i = 0; i < alleles_adaptation_coef.size(); ++i)
  {
    fitnesslike_per_allele.push_back(abs(alleles_adaptation_coef[i] - cell_temperature));
  }
  unfitnesslike_thispopulation = 0;
  unfitnesslike_thispopulation = unfitnesslike_thispopulation + (this_pop.allelic_a * fitnesslike_per_allele[0]);
  unfitnesslike_thispopulation = unfitnesslike_thispopulation + (this_pop.allelic_b * fitnesslike_per_allele[1]);
  unfitnesslike_thispopulation = unfitnesslike_thispopulation + (this_pop.allelic_c * fitnesslike_per_allele[2]);
  unfitnesslike_thispopulation = unfitnesslike_thispopulation + (this_pop.allelic_d * fitnesslike_per_allele[3]);
  unfitnesslike_thispopulation = unfitnesslike_thispopulation + (this_pop.allelic_e * fitnesslike_per_allele[4]);


  int max_unfit;
  max_unfit = *max_element(fitnesslike_per_allele.begin(),fitnesslike_per_allele.end());
  max_unfit = max_unfit + 1;
  double average_unfitness;

  average_unfitness = unfitnesslike_thispopulation/this_pop.pop_size;

  // to have a sort of index that goes from 0 to 1 where 1 is perfectly adapted to a given cell
  fitnesslike_thispopulation = (max_unfit - average_unfitness)/max_unfit;

  // cout << "      +++++++ fitnesslike_thispopulation :" << fitnesslike_thispopulation << endl;
  return fitnesslike_thispopulation;
  // vector <int> all_populations_fitnesslike;
  //
  // for (int ii = 0; ii < use_this_species.populations_this_species.size(); ++ii){
  //   int cell_temperature;
  //   cell_temperature = map1[use_this_species.presence[ii].y - 1][use_this_species.presence[ii].x - 1].temperature;
  //    vector <int> fitnesslike_per_allele;
  //
  //
  //
  //   vector <int> fitnesslike_thispopulation;
  //   fitnesslike_thispopulation = 0;
  //
  //
  //   fitnesslike_thispopulation = fitnesslike_thispopulation + (use_this_species.populations_this_species[ii].allelic_a * fitnesslike_per_allele[0]);
  //   fitnesslike_thispopulation = fitnesslike_thispopulation + (use_this_species.populations_this_species[ii].allelic_b * fitnesslike_per_allele[1]);
  //   fitnesslike_thispopulation = fitnesslike_thispopulation + (use_this_species.populations_this_species[ii].allelic_c * fitnesslike_per_allele[2]);
  //   fitnesslike_thispopulation = fitnesslike_thispopulation + (use_this_species.populations_this_species[ii].allelic_d * fitnesslike_per_allele[3]);
  //   fitnesslike_thispopulation = fitnesslike_thispopulation + (use_this_species.populations_this_species[ii].allelic_e * fitnesslike_per_allele[4]);
  //
  //   all_populations_fitnesslike.push_back(fitnesslike_thispopulation)
  // }
  //
  //
  //
  // //  //  old bit for trait evolution
  // // //cout << "linkin function: " << gamma << endl;
  // // // take a value of mu
  // // double calculated_fitness;
  // //
  // // calculated_fitness = exp(-((species_optimum_temperature-cell_temperature)*(species_optimum_temperature-cell_temperature))/(2*v) );
  // // //  calculated_fitness = abs(cell_temperature - species_optimum_temperature) ;
  // return all_populations_fitnesslike;
}

void species::initial_position(yx initial)
{
  presence.push_back(initial);
}

void species::happening_trait_evolution(double mean_normal_distribution, double sd_normal_distribution)
{
  double evol_rate;
  evol_rate = give_me_random_normal(mean_normal_distribution, sd_normal_distribution);
  trait_state = trait_state + evol_rate;
}

vector<yx> species::find_neighbor(int cell)
{

  yx neighbor_a;
  yx neighbor_b;
  yx neighbor_c;
  yx neighbor_d;
  yx neighbor_e;
  yx neighbor_f;
  yx neighbor_g;
  yx neighbor_h;

  yx cell_sending_colonizer;
  vector<yx> adjacent_cells;
  cell_sending_colonizer = presence[cell];
  // cout << "cell_sending_colonizer xy" << cell_sending_colonizer.x <<" "<< cell_sending_colonizer.y  <<endl;

  neighbor_a.x = cell_sending_colonizer.x - 1;
  neighbor_b.x = cell_sending_colonizer.x - 1;
  neighbor_c.x = cell_sending_colonizer.x - 1;
  neighbor_d.x = cell_sending_colonizer.x;
  neighbor_e.x = cell_sending_colonizer.x;
  neighbor_f.x = cell_sending_colonizer.x + 1;
  neighbor_g.x = cell_sending_colonizer.x + 1;
  neighbor_h.x = cell_sending_colonizer.x + 1;

  neighbor_a.y = cell_sending_colonizer.y - 1;
  neighbor_b.y = cell_sending_colonizer.y;
  neighbor_c.y = cell_sending_colonizer.y + 1;
  neighbor_d.y = cell_sending_colonizer.y - 1;
  neighbor_e.y = cell_sending_colonizer.y + 1;
  neighbor_f.y = cell_sending_colonizer.y - 1;
  neighbor_g.y = cell_sending_colonizer.y;
  neighbor_h.y = cell_sending_colonizer.y + 1;

  // cout << "focal cell:" << cell_sending_colonizer.x << "," << cell_sending_colonizer.y << endl;

  if ((neighbor_a.x > 0) && (neighbor_a.x < (x_max + 1)) && (neighbor_a.y > 0) && (neighbor_a.x < (y_max + 1)))
  {
    adjacent_cells.push_back(neighbor_a);
  }
  if ((neighbor_b.x > 0) && (neighbor_b.x < (x_max + 1)) && (neighbor_b.y > 0) && (neighbor_b.x < (y_max + 1)))
  {
    adjacent_cells.push_back(neighbor_b);
  }
  if ((neighbor_c.x > 0) && (neighbor_c.x < (x_max + 1)) && (neighbor_c.y > 0) && (neighbor_c.x < (y_max + 1)))
  {
    adjacent_cells.push_back(neighbor_c);
  }
  if ((neighbor_d.x > 0) && (neighbor_d.x < (x_max + 1)) && (neighbor_d.y > 0) && (neighbor_d.x < (y_max + 1)))
  {
    adjacent_cells.push_back(neighbor_d);
  }
  if ((neighbor_e.x > 0) && (neighbor_e.x < (x_max + 1)) && (neighbor_e.y > 0) && (neighbor_e.x < (y_max + 1)))
  {
    adjacent_cells.push_back(neighbor_e);
  }
  if ((neighbor_f.x > 0) && (neighbor_f.x < (x_max + 1)) && (neighbor_f.y > 0) && (neighbor_f.x < (y_max + 1)))
  {
    adjacent_cells.push_back(neighbor_f);
  }
  if ((neighbor_g.x > 0) && (neighbor_g.x < (x_max + 1)) && (neighbor_g.y > 0) && (neighbor_g.x < (y_max + 1)))
  {
    adjacent_cells.push_back(neighbor_g);
  }
  if ((neighbor_h.x > 0) && (neighbor_h.x < (x_max + 1)) && (neighbor_h.y > 0) && (neighbor_h.x < (y_max + 1)))
  {
    adjacent_cells.push_back(neighbor_h);
  }
  //cout << "number of neighbor: " << adjacent_cells.size() << endl;
  return adjacent_cells;
}

probabilities_based_traits calculate_probabilities_using_traitstate(vector<species> all_species, landscape **map1, std::string temperature_influencing, double mutation_rate, double geneflow_rate, double popchange_rate, double lambda, double gamma, double mu, vector<int> id_alive_species, double v)
{
  // int species_to_do;
  // species_to_do = id_alive_species[give_me_random_uniform(0, (id_alive_species.size()-1))];

  double gammas_total = 0;
  double lambdas_total = 0;
  double mus_total = 0;
  double popchange_rate_total = 0;
  double geneflow_rate_total = 0;
  double mutation_rate_total = 0;
  double differences_temperature = 0;
  int total_populations = 0;
  vector<double> total_probability_species;
  for (int i = 0; i < all_species.size(); ++i)
  {
    species work_this_species;
    work_this_species = all_species[i];
    double total_computed_rate_based_on_temperature;
    double this_species_lambda;
    double this_species_gamma;
    double this_species_mu;
    double this_species_popchange;
    double this_species_geneflow;
    double this_species_mutation;
    total_computed_rate_based_on_temperature = 0;
    if (work_this_species.alive)
    {
      // work_this_species.computed_rate_based_on_temperature.clear();
      double total_fitnesslike;
      total_fitnesslike = std::accumulate(work_this_species.computed_rate_based_on_temperature.begin(), work_this_species.computed_rate_based_on_temperature.end(),
                                          decltype(work_this_species.computed_rate_based_on_temperature)::value_type(0));
      if (temperature_influencing == "mu_depending")
      {
        this_species_mu = total_fitnesslike * mu;
        this_species_gamma = work_this_species.range * gamma;
      }
      if (temperature_influencing == "gamma_depending" || temperature_influencing == "gamma_depending_mu_onPopsize")
      {
        this_species_gamma = (total_fitnesslike/work_this_species.range) * gamma;

        // being poorly adapted increases the chances that extirpation happens, so this event is more likely to happen
        // I use the range size of the species as that value would be equivalent to populations with the highest
        // fitness_like score
        this_species_mu = ((work_this_species.range - total_fitnesslike)/work_this_species.range) * mu * work_this_species.range;
      }
      this_species_geneflow = work_this_species.range * geneflow_rate;
      this_species_lambda = work_this_species.range * lambda;
      this_species_popchange = work_this_species.range * popchange_rate;

      this_species_mutation = work_this_species.range * mutation_rate;
      // cout << "probs this species.  mu: " << this_species_mu
      //
      // << " gamma: " << this_species_gamma
      //
      // << " lambda: " << this_species_lambda
      //
      // << " popchange: " << this_species_popchange <<
      // " mutation: " << this_species_mutation << " total fitness: " << total_fitnesslike << endl;

      mus_total = mus_total + this_species_mu;
      gammas_total = gammas_total + this_species_gamma;
      lambdas_total = lambdas_total + this_species_lambda;
      popchange_rate_total = popchange_rate_total + this_species_popchange;
      geneflow_rate_total = geneflow_rate_total + this_species_geneflow;
      mutation_rate_total = mutation_rate_total + this_species_mutation;
      work_this_species.total_rate = this_species_mu + this_species_gamma + this_species_lambda + this_species_popchange + this_species_geneflow + this_species_mutation;
      total_probability_species.push_back(work_this_species.total_rate);

      all_species[i] = work_this_species;
    }
  }

  double total_waiting_times_rates;
  total_waiting_times_rates = gammas_total + lambdas_total + mus_total + popchange_rate_total + geneflow_rate_total + mutation_rate_total;
  probabilities_based_traits calculation_probabilities;
  calculation_probabilities.total_probability_species = total_probability_species;
  calculation_probabilities.gammas_total = gammas_total;
  calculation_probabilities.lambdas_total = lambdas_total;
  calculation_probabilities.mus_total = mus_total;
  calculation_probabilities.popchange_rate_total = popchange_rate_total;
  calculation_probabilities.geneflow_rate_total = geneflow_rate_total;
  calculation_probabilities.mutation_rate_total = mutation_rate_total;
  calculation_probabilities.total_waiting_times_rates = total_waiting_times_rates;
  return calculation_probabilities;
}

void to_show_richness_map(std::string show_richness_map,vector<species>all_species, landscape **map1)
{
  // create_dynamic_map ( x_max, y_max);
  //populate_landscape(all_species, map1);
  cout << endl;
  // it can be trait, richness or none
  if (show_richness_map == "both")
  {
    cout << endl;
    cout << "this is total_abundance map" << endl;
    for (int i = 0; i < y_max; i++)
    {
      for (int j = 0; j < x_max; j++)
      {
        cout << map1[i][j].total_abundance_cell << "\t";
      }
      cout << endl;
    }
    cout << endl;
    cout << "this is richness map" << endl;
    for (int i = 0; i < y_max; i++)
    {
      for (int j = 0; j < x_max; j++)
      {
        cout << map1[i][j].Nspp << "\t";
      }
      cout << endl;
    }
  }
  if (show_richness_map == "total_abundance")
  {
    cout << endl;
    cout << "this is total_abundance map" << endl;
    for (int i = 0; i < y_max; i++)
    {
      for (int j = 0; j < x_max; j++)
      {
        cout << map1[i][j].total_abundance_cell << "\t";
      }
      cout << endl;
    }
  }
  if (show_richness_map == "richness")
  {
    cout << endl;
    cout << "this is richness map" << endl;
    for (int i = 0; i < y_max; i++)
    {
      for (int j = 0; j < x_max; j++)
      {
        cout << map1[i][j].Nspp << "\t";
      }
      cout << endl;
    }
  }
  if (show_richness_map == "richness")
  {
    cout << endl;
  }
  if (show_richness_map == "none")
  {
    cout << "no map is requested to be shown" << endl;
  }
}

void populate_landscape(vector<species> all_species, landscape **map1)
{
  for (int ii = 0; ii < all_species.size(); ii++)
  {
    // cout << "plotting the species" << ii << endl;
    species focal;
    focal = all_species[ii];
    // richness
    vector<yx> distribution_focal = focal.presence;
    for (int i = 0; i < distribution_focal.size(); i++)
    {
      // I substract 1 as it is and index
      int current_richness;
      double per_cell_trait;
      int current_abundance;
      current_richness = map1[(distribution_focal[i].y - 1)][(distribution_focal[i].x - 1)].Nspp;
      map1[(distribution_focal[i].y - 1)][(distribution_focal[i].x - 1)].Nspp = current_richness + 1;
      // average trait
      per_cell_trait = map1[(distribution_focal[i].y - 1)][(distribution_focal[i].x - 1)].sum_trait;
      map1[(distribution_focal[i].y - 1)][(distribution_focal[i].x - 1)].sum_trait = per_cell_trait + focal.trait_state;

      current_abundance = map1[(distribution_focal[i].y - 1)][(distribution_focal[i].x - 1)].total_abundance_cell;
      map1[(distribution_focal[i].y - 1)][(distribution_focal[i].x - 1)].total_abundance_cell = current_abundance + focal.populations_this_species[i].pop_size;
    }
  }

  for (int i = 0; i < y_max; i++)
  {
    for (int j = 0; j < x_max; j++)
    {
      if (map1[i][j].Nspp > 0)
      {
        map1[i][j].average_trait = map1[i][j].sum_trait / map1[i][j].Nspp;
      }
    }
  }
}
// function to make a subset of neighbouring cells that can be colonized

vector<yx> species::available_neigh_to_colonize_K(int x_max, int y_max, int focal_cell, landscape **map1)
{
  vector<yx> neighbors_focal;
  neighbors_focal = find_neighbor(focal_cell);
  // to get rid of the -9 from map
  vector<yx> neighbors_focal_good_habitat;
  for (int i = 0; i < neighbors_focal.size(); i++)
  {
    int type_habitat;
    // cout << "coordinate x: " << neighbors_focal[i].x << " and y: " << neighbors_focal[i].y << endl;
    type_habitat = map1[neighbors_focal[i].y - 1][neighbors_focal[i].x - 1].habitat;
    if (type_habitat != -9)
    {
      neighbors_focal_good_habitat.push_back(neighbors_focal[i]);
    }
  }
  // to remove cells where focal is present. I will artificially, put -9 to those cells
  // where the focal species inhabits. Then, I will use the code above to remove those cells.
  vector<yx> neighbors_focal_good_habitat2;

  landscape map_artificial[y_max][x_max];
  for (int i = 0; i < y_max; i++)
  {
    for (int j = 0; j < x_max; j++)
    {
      map_artificial[i][j].habitat = 0;
    }
  }
  vector<yx> distribution_focal = presence;
  for (int i = 0; i < distribution_focal.size(); i++)
  {
    // I substract 1 as it is and index
    int current_richness;
    map_artificial[(distribution_focal[i].y - 1)][(distribution_focal[i].x - 1)].habitat = -9;
  }
  // place in neighbors_focal_good_habitat2 only the ones where focal does not live

  for (int i = 0; i < neighbors_focal_good_habitat.size(); i++)
  {
    int type_habitat;
    // cout << "coordinate x: " << neighbors_focal[i].x << " and y: " << neighbors_focal[i].y << endl;
    type_habitat = map_artificial[neighbors_focal_good_habitat[i].y - 1][neighbors_focal_good_habitat[i].x - 1].habitat;
    //cout << "habitat cell: " << type_habitat << endl;
    if (type_habitat != -9)
    {
      neighbors_focal_good_habitat2.push_back(neighbors_focal_good_habitat[i]);
    }
  }
  // to find those cells under K
  vector<yx> available_to_colonize;

  for (int i = 0; i < neighbors_focal_good_habitat2.size(); i++)
  {
    int richness_cell;
    int k_patch;
    int cell_temperature;
    int total_abundance_cell_cell;
    cell_temperature = map1[neighbors_focal_good_habitat2[i].y - 1][neighbors_focal_good_habitat2[i].x - 1].temperature;
    richness_cell = map1[neighbors_focal_good_habitat2[i].y - 1][neighbors_focal_good_habitat2[i].x - 1].Nspp;
    k_patch = map1[neighbors_focal_good_habitat2[i].y - 1][neighbors_focal_good_habitat2[i].x - 1].k_patch;
    total_abundance_cell_cell = map1[neighbors_focal_good_habitat2[i].y - 1][neighbors_focal_good_habitat2[i].x - 1].total_abundance_cell;
    // cout << "this is k patch" << k_patch << "and this is total_abundance_cell_cell: " << total_abundance_cell_cell << endl;
    // cout << "this is k patch" << k_patch << "and this is richness_cell: " << richness_cell << endl;
    // cout << "this is k patch" << k_patch << "and this is cell_temperature: " << cell_temperature << endl;
    if (total_abundance_cell_cell < k_patch)
    {
      //  cout <<  "this one is available and under K" << endl;
      available_to_colonize.push_back(neighbors_focal_good_habitat2[i]);
    }
  }
  //cout << "num cells to colonize:" << available_to_colonize.size() << endl;
  return available_to_colonize;
}

vector<yx> species::available_neigh_to_colonize_trait(int x_max, int y_max, int focal_cell, double trait_dissimilarity_threshold, landscape **map1)
{
  vector<yx> neighbors_focal;
  neighbors_focal = find_neighbor(focal_cell);
  // to get rid of the -9 from map
  vector<yx> neighbors_focal_good_habitat;
  for (int i = 0; i < neighbors_focal.size(); i++)
  {
    int type_habitat;
    // cout << "coordinate x: " << neighbors_focal[i].x << " and y: " << neighbors_focal[i].y << endl;
    type_habitat = map1[neighbors_focal[i].y - 1][neighbors_focal[i].x - 1].habitat;
    // cout << "habitat " << type_habitat << endl;
    if (type_habitat != -9)
    {
      neighbors_focal_good_habitat.push_back(neighbors_focal[i]);
    }
  }
  // to remove cells where focal is present. I will artificially, put -9 to those cells
  // where the focal species inhabits. Then, I will use the code above to remove those cells.
  vector<yx> neighbors_focal_good_habitat2;
  landscape map_artificial[y_max][x_max];
  for (int i = 0; i < y_max; i++)
  {
    for (int j = 0; j < x_max; j++)
    {
      map_artificial[i][j].habitat = 0;
    }
  }
  // richness
  vector<yx> distribution_focal = presence;
  for (int i = 0; i < distribution_focal.size(); i++)
  {
    // I substract 1 as it is and index
    int current_richness;
    map_artificial[(distribution_focal[i].y - 1)][(distribution_focal[i].x - 1)].habitat = -9;
  }
  // place in neighbors_focal_good_habitat2 only the ones where focal does not live
  for (int i = 0; i < neighbors_focal_good_habitat.size(); i++)
  {
    int type_habitat;
    // cout << "coordinate x: " << neighbors_focal[i].x << " and y: " << neighbors_focal[i].y << endl;
    type_habitat = map_artificial[neighbors_focal_good_habitat[i].y - 1][neighbors_focal_good_habitat[i].x - 1].habitat;
    if (type_habitat != 9)
    {
      neighbors_focal_good_habitat2.push_back(neighbors_focal_good_habitat[i]);
    }
  }
  // find cells under dissimilary threshold
  vector<yx> available_to_colonize;
  for (int i = 0; i < neighbors_focal_good_habitat2.size(); i++)
  {
    double per_cell_trait;
    // cout << "coordinate x: " << neighbors_focal[i].x << " and y: " << neighbors_focal[i].y << endl;
    per_cell_trait = map1[neighbors_focal_good_habitat2[i].y - 1][neighbors_focal_good_habitat2[i].x - 1].average_trait;

    double percentage_difference; // Percentage difference formula
    percentage_difference = abs(trait_state - per_cell_trait) / ((trait_state + per_cell_trait) / 2);
    // cout << "the trait is: " << trait_state << " and the average of the community is: " << per_cell_trait << endl;
    // cout << "                           the percentage_difference is " << percentage_difference << " and I will";
    if (percentage_difference >= (trait_dissimilarity_threshold / 100.0))
    {
      available_to_colonize.push_back(neighbors_focal_good_habitat2[i]);
      //   cout << "Include it" << endl;
    }
    else
    {
      //   cout << "NOT include it" << endl;
      //   cout << neighbors_focal_good_habitat[i].x << neighbors_focal_good_habitat[i].y << endl;
    }
  }
  // cout << "number avai cells, under local limit (either k or trait dissimilarity): " << available_to_colonize.size() << endl;
  for (int ii = 0; ii < available_to_colonize.size(); ii++)
  {
    //  cout << "coordinate: " << available_to_colonize[ii].x << "," << available_to_colonize[ii].y << endl;
  }
  return available_to_colonize;
}

vector <int> give_me_random_wallenius(int number_alelles,vector<int> current_allelic_frequency, vector<int> weights_for_wallenius, int number_items_sample)
{

  // int number_alelles;
  // number_alelles = current_allelic_frequency.size();
  int32_t frequencies[number_alelles];         // output vector
  int32_t mlist[number_alelles];// number of balls of each color

  int total_poolsize;
  total_poolsize = 0;
  for(int ii = 0; ii < current_allelic_frequency.size(); ++ii){
    mlist[ii] = current_allelic_frequency[ii];
    total_poolsize = total_poolsize + current_allelic_frequency[ii];
  }
  // double wlist[] ={1, 5, 3, 6};   // weight of each color
  double wlist[number_alelles]; // weight of each color

  for(int ii = 0; ii < weights_for_wallenius.size(); ++ii){
    if(weights_for_wallenius[ii] == 0){
      wlist[ii] = 0.1;// in the case that only that allele is present, it needs to be selected anyway
    } else {
      wlist[ii] = weights_for_wallenius[ii];
    }

  }

  int32_t n = number_items_sample;        // number of balls to pick
  StochasticLib3 randomwallenius((int32_t)time(0));

  if(number_items_sample > total_poolsize){
    for(int ii = 0; ii < weights_for_wallenius.size(); ++ii){

      cout << "current_allelic_frequency: " << current_allelic_frequency[ii] << "and weight: " << weights_for_wallenius[ii]<<endl;
    }
    cout << "and this is the number of items to sample: " << number_items_sample << endl;
  }

  randomwallenius.MultiWalleniusNCHyp(frequencies, mlist, wlist, n, number_alelles);
  vector <int> frequencies2;

  for(int ii = 0; ii < current_allelic_frequency.size(); ++ii){
    frequencies2.push_back(frequencies[ii]);
  }
  return frequencies2;
}




int give_me_random_uniform(int from, int to)
{
  random_device rd;
  default_random_engine generator1(rd());

  int my_random;
  uniform_int_distribution<> dis(from, to);
  my_random = dis(generator1);
  // cout << "random:" << my_random << endl;
  return my_random;
}

double give_me_random_normal(double my_mean, double my_sd)
{
  random_device rd;
  default_random_engine generator2(rd());
  double my_random;
  std::normal_distribution<> dis(my_mean, my_sd);
  my_random = dis(generator2);
  // cout << "random:" << my_random << endl;
  return my_random;
}

void species::happening_contraction(double t, landscape **map1, std::string temperature_influencing)
{
  random_device rd;
  default_random_engine generator3(rd());
  int random_cell_to_remove_population;
  if (temperature_influencing == "gamma_depending")
  {
    random_cell_to_remove_population = give_me_random_uniform(0, presence.size() - 1);
  }
  if (temperature_influencing == "gamma_depending_mu_onPopsize")
  {

    vector<int> all_cells_id;
    vector <int> pop_sizes_vector;
    for (int i = 0; i < presence.size(); ++i)
    {
      all_cells_id.push_back(i);
      // give higher probabilities to those populations better adapted

      //cout << "probab to be selected for contraction " << (1 - computed_rate_based_on_temperature[i]) << endl;
      pop_sizes_vector.push_back((1 - computed_rate_based_on_temperature[i]));
    }

    discrete_distribution<int> cell_probabilities_to_pick(pop_sizes_vector.begin(), pop_sizes_vector.end());
    random_cell_to_remove_population = all_cells_id[cell_probabilities_to_pick(generator3)];
  }

  x_coordinate_last_event = presence[random_cell_to_remove_population].x;
  y_coordinate_last_event = presence[random_cell_to_remove_population].y;
  // classify_elevation( false, map1);

  // update the map
  int current_richness_map;
  current_richness_map = map1[(presence[random_cell_to_remove_population].y - 1)][(presence[random_cell_to_remove_population].x - 1)].Nspp;
  map1[(presence[random_cell_to_remove_population].y - 1)][(presence[random_cell_to_remove_population].x - 1)].Nspp = current_richness_map - 1;

  int current_popsize_map;
  current_popsize_map = map1[(presence[random_cell_to_remove_population].y - 1)][(presence[random_cell_to_remove_population].x - 1)].total_abundance_cell;
  map1[(presence[random_cell_to_remove_population].y - 1)][(presence[random_cell_to_remove_population].x - 1)].total_abundance_cell = current_popsize_map - populations_this_species[random_cell_to_remove_population].pop_size;


  // cout << " ++ cells presence Before contraction:" << presence.size() << endl;
  presence.erase(presence.begin() + random_cell_to_remove_population);
  temperature_optimum.erase(temperature_optimum.begin() + random_cell_to_remove_population);
  computed_rate_based_on_temperature.erase(computed_rate_based_on_temperature.begin() + random_cell_to_remove_population);
  range = range - 1;
  // cout << "total_pop_size from here: "<< total_pop_size << endl;
  // cout << "population size to be removed: " << populations_this_species[random_cell_to_remove_population].pop_size <<endl;
  total_pop_size = total_pop_size - populations_this_species[random_cell_to_remove_population].pop_size;
  // cout << "total_pop_size from here: "<< total_pop_size << endl;
  populations_this_species.erase(populations_this_species.begin() + random_cell_to_remove_population);


  // part to see whether there was a latitudinal contraction

  // if(presence[random_cell_to_remove_population].y == (southernmost + 1))
  // {
  //   change_southernmost.push_back(-1);
  //   time_change_southernmost.push_back(t);
  //   southernmost = southernmost + 1;
  // }
  // if(presence[random_cell_to_remove_population].y == (northernmost - 1))
  // {
  //   change_northernmost.push_back(-1);
  //   time_change_northernmost.push_back(t);
  //   northernmost = northernmost - 1;
  // }



  // cout << " ++ cells presence After contraction:" << presence.size() << endl;
  if (range == 0)
  { // it died by contraction too much
    death = t;
    alive = false;
  }
  else
  {
    update_latitudinal_borders(t,false);
  }
}

void species::happening_expansion(int x_max, int y_max, bool use_k, double restiction_par, landscape **map1, std::string temperature_influencing, vector<int> alleles_adaptation_coef,double t)
{ // restriction par will be either k or trait dissimilarity
  // If there is any cell available for focal species to expand, colonization will take place
  // empty cells that are surrounded by focal species, will be more likely to be colonized
  random_device rd;
  default_random_engine generator4(rd());
  vector<yx> ready_to_colonize2;
  vector<double> probabilities_temperature2;
  vector<int> id_cells2;
  vector<int> all_cells_id;
  for (int i = 0; i < presence.size(); ++i)
  {
    vector<yx> ready_to_colonize;
    double probabilities_temperature;
    int id_cells;
    if (use_k == true)
    {
      ready_to_colonize = available_neigh_to_colonize_K(x_max, y_max, i, map1);
      probabilities_temperature = computed_rate_based_on_temperature[i];
      id_cells = i;
    }
    else
    {
      //  cout << "colonization based on trait similarity" << endl;
      ready_to_colonize = available_neigh_to_colonize_trait(x_max, y_max, i, restiction_par, map1);
    }
    for (int ii = 0; ii < ready_to_colonize.size(); ++ii)
    {
      ready_to_colonize2.push_back(ready_to_colonize[ii]);
      probabilities_temperature2.push_back(probabilities_temperature);
      id_cells2.push_back(id_cells);
      all_cells_id.push_back((all_cells_id.size() + 1)); // all_cells_id is just a sequence of integers of equal length than probabilities_temperature2
    }
  }
  // cout << "to colonize: " << ready_to_colonize2.size() << endl;
  if (ready_to_colonize2.size() > 0)
  { // it has, otherwise, no where to go
    int random_cell_to_go;
    if (temperature_influencing == "mu_depending")
    {
      random_cell_to_go = give_me_random_uniform(1, ready_to_colonize2.size());

    }
    if (temperature_influencing == "gamma_depending" || temperature_influencing == "gamma_depending_mu_onPopsize")
    {
      discrete_distribution<int> cell_probabilities_to_pick(probabilities_temperature2.begin(), probabilities_temperature2.end());
      random_cell_to_go = all_cells_id[cell_probabilities_to_pick(generator4)];

    }
    if (temperature_influencing == "independent")
    {
      random_cell_to_go = give_me_random_uniform(1, ready_to_colonize2.size());

    }
    population_structure current_population_sender;
    current_population_sender = populations_this_species[id_cells2[random_cell_to_go - 1]];
    //cout << " show current_population_sender.pop_size " << current_population_sender.pop_size << endl;
    if(current_population_sender.pop_size >= 2)
    {
      presence.push_back(ready_to_colonize2[random_cell_to_go - 1]); // -1 as it is index

      // log if latitudinal expansion took place
      // cout << "y coor: "<<ready_to_colonize2[random_cell_to_go - 1].y << endl;

      //    if (ready_to_colonize2[random_cell_to_go - 1].y == (northernmost + 1)) // northward expansion
      //    {
      //    //  cout << "here in north" << endl;
      //      change_northernmost.push_back(1); // Ones are expansions whereas -1 are contractions
      //      time_change_northernmost.push_back(t);
      // //  cout << "before " << northernmost << endl;
      //      northernmost = northernmost + 1;
      //   //   cout << "after " << northernmost << endl;
      //    }
      //    if (ready_to_colonize2[random_cell_to_go - 1].y == (southernmost - 1)) // southward expansion
      //    {
      //    //  cout << "here in south" << endl;
      //      change_southernmost.push_back(1); // Ones are expansions whereas -1 are contractions
      //      time_change_southernmost.push_back(t);
      //      //cout << "before " << southernmost << endl;
      //      southernmost = southernmost - 1;
      //     // cout << "after " << southernmost << endl;
      //    }


      // the part where the new population could have evolved a different temperature from parental population
      double evol_rate;
      double current_temperature_sender;
      double temperature_new_population;
      // evol_rate = give_me_random_normal(mean_normal_distribution, sd_normal_distribution);
      current_temperature_sender = temperature_optimum[id_cells2[random_cell_to_go - 1]];
      //    cout << "evol_rate: " << evol_rate << endl;
      temperature_new_population = current_temperature_sender + evol_rate;
      temperature_optimum.push_back(temperature_new_population);
      //        cout << "sender: " << current_temperature_sender + evol_rate << endl;
      //   // computed_rate_based_on_temperature.push_back(probabilities_temperature2[id_cells2[random_cell_to_go]]);
      //
      // part that adds the new population and splits population size between parent and new population. Also separates allelic frequencies

      population_structure new_population;
      int percentage_to_split_pop;
      double new_pop_size;
      percentage_to_split_pop = give_me_random_uniform(10, 80);
      //  cout << "percentage_to_split_pop " << percentage_to_split_pop << endl;
      new_pop_size = (percentage_to_split_pop / 100.0) * current_population_sender.pop_size;
      //   cout << "new_pop_size from int: " << new_pop_size << endl;
      new_pop_size = floor(new_pop_size);
      if(new_pop_size < 1)
      {
        new_pop_size = 1;
      }
      int original_population_sender_pop_size;
      original_population_sender_pop_size = current_population_sender.pop_size;
      //  cout << "current_population_sender.pop_size before: " << current_population_sender.pop_size << endl;

      // the number of individuals going to the cell should not exceed K

      int abundance_individuals_cell_to_go;
      abundance_individuals_cell_to_go = map1[(ready_to_colonize2[random_cell_to_go - 1].y - 1)][(ready_to_colonize2[random_cell_to_go - 1].x - 1)].total_abundance_cell;
      int k_cell;
      k_cell = map1[(ready_to_colonize2[random_cell_to_go - 1].y - 1)][(ready_to_colonize2[random_cell_to_go - 1].x - 1)].k_patch;
      if((new_pop_size + abundance_individuals_cell_to_go) > k_cell )
      {
        new_pop_size = k_cell - abundance_individuals_cell_to_go;
      }


      int cell_temperature;
      cell_temperature = map1[ready_to_colonize2[random_cell_to_go - 1].y - 1][ready_to_colonize2[random_cell_to_go - 1].x - 1].temperature;
      vector<int> weights_for_wallenius;
      vector<int> neutral_weights_for_wallenius;
      vector<int> difference_cell_allele_adapt;
      for (int ii = 0; ii < alleles_adaptation_coef.size(); ++ii)
      {
        difference_cell_allele_adapt.push_back(abs(alleles_adaptation_coef[ii] - cell_temperature));
      }

      int max_unfit;
      max_unfit = *max_element(difference_cell_allele_adapt.begin(),difference_cell_allele_adapt.end());

      for (int ii = 0; ii < alleles_adaptation_coef.size(); ++ii)
      {
        weights_for_wallenius.push_back((max_unfit - abs(alleles_adaptation_coef[ii] - cell_temperature)) + 1);
        neutral_weights_for_wallenius.push_back(1);
      }



      vector<int> current_allelic_frequency;
      current_allelic_frequency.push_back(current_population_sender.allelic_a);
      current_allelic_frequency.push_back(current_population_sender.allelic_b);
      current_allelic_frequency.push_back(current_population_sender.allelic_c);
      current_allelic_frequency.push_back(current_population_sender.allelic_d);
      current_allelic_frequency.push_back(current_population_sender.allelic_e);

      vector<int> current_allelic_frequency_neutral;
      current_allelic_frequency_neutral.push_back(current_population_sender.allelic_v);
      current_allelic_frequency_neutral.push_back(current_population_sender.allelic_w);
      current_allelic_frequency_neutral.push_back(current_population_sender.allelic_x);
      current_allelic_frequency_neutral.push_back(current_population_sender.allelic_y);
      current_allelic_frequency_neutral.push_back(current_population_sender.allelic_z);
      //
      //       cout << "current_allelic_frequency A "<<current_allelic_frequency[0] << endl;
      //       cout << "current_allelic_frequency B "<<current_allelic_frequency[1] << endl;
      //       cout << "current_allelic_frequency C "<<current_allelic_frequency[2] << endl;
      //       cout << "current_allelic_frequency D "<<current_allelic_frequency[3] << endl;
      //       cout << "current_allelic_frequency E "<<current_allelic_frequency[4] << endl;


      // cout << "weights_for_wallenius A "<<weights_for_wallenius[0] << endl;
      // cout << "weights_for_wallenius B "<<weights_for_wallenius[1] << endl;
      // cout << "weights_for_wallenius C "<<weights_for_wallenius[2] << endl;
      // cout << "weights_for_wallenius D "<<weights_for_wallenius[3] << endl;
      // cout << "weights_for_wallenius E "<<weights_for_wallenius[4] << endl;

      int number_alelles;
      number_alelles = weights_for_wallenius.size();
      vector<int> allelic_frequency_new;
      vector<int> allelic_frequency_new_neutral;


      int pop_size_from_sum_alle;
      pop_size_from_sum_alle = std::accumulate(current_allelic_frequency.begin(),current_allelic_frequency.end(),
                                               decltype(current_allelic_frequency)::value_type(0));
      //
      //       cout << "pop_size_from_sum_alle " << pop_size_from_sum_alle << endl;
      //       cout << "new pop size this point : " << new_pop_size << endl;
      //       cout << "number_alelles" << number_alelles << endl;
      // cout << "vector_alleliccurrent "<<current_allelic_frequency[0] << current_allelic_frequency[1] <<current_allelic_frequency[2] << current_allelic_frequency[3] <<current_allelic_frequency[4] << endl;
      allelic_frequency_new = give_me_random_wallenius(number_alelles,current_allelic_frequency, weights_for_wallenius, new_pop_size);
      // cout << "now neutral" << endl;
      // cout << "vector_alleliccurrent NEUTRAL"<<current_allelic_frequency_neutral[0] << current_allelic_frequency_neutral[1] <<current_allelic_frequency_neutral[2] << current_allelic_frequency_neutral[3] <<current_allelic_frequency_neutral[4] << endl;
      allelic_frequency_new_neutral = give_me_random_wallenius(number_alelles,current_allelic_frequency_neutral, neutral_weights_for_wallenius, new_pop_size);


      new_population.allelic_a = allelic_frequency_new[0];
      new_population.allelic_b = allelic_frequency_new[1];
      new_population.allelic_c = allelic_frequency_new[2];
      new_population.allelic_d = allelic_frequency_new[3];
      new_population.allelic_e = allelic_frequency_new[4];

      new_population.allelic_v = allelic_frequency_new_neutral[0];
      new_population.allelic_w = allelic_frequency_new_neutral[1];
      new_population.allelic_x = allelic_frequency_new_neutral[2];
      new_population.allelic_y = allelic_frequency_new_neutral[3];
      new_population.allelic_z = allelic_frequency_new_neutral[4];

      new_population.pop_size = new_pop_size;
      // cout << "new pop size here: " << new_population.pop_size << endl;
      // cout << "new_population a" << new_population.allelic_a << endl;
      // cout << "new_population b" << new_population.allelic_b << endl;
      // cout << "new_population c" << new_population.allelic_c << endl;
      // cout << "new_population d" << new_population.allelic_d << endl;
      // cout << "new_population e" << new_population.allelic_e << endl;

      current_population_sender.allelic_a = current_population_sender.allelic_a - allelic_frequency_new[0];
      current_population_sender.allelic_b = current_population_sender.allelic_b - allelic_frequency_new[1];
      current_population_sender.allelic_c = current_population_sender.allelic_c - allelic_frequency_new[2];
      current_population_sender.allelic_d = current_population_sender.allelic_d - allelic_frequency_new[3];
      current_population_sender.allelic_e = current_population_sender.allelic_e - allelic_frequency_new[4];

      current_population_sender.allelic_v = current_population_sender.allelic_v - allelic_frequency_new_neutral[0];
      current_population_sender.allelic_w = current_population_sender.allelic_w - allelic_frequency_new_neutral[1];
      current_population_sender.allelic_x = current_population_sender.allelic_x - allelic_frequency_new_neutral[2];
      current_population_sender.allelic_y = current_population_sender.allelic_y - allelic_frequency_new_neutral[3];
      current_population_sender.allelic_z = current_population_sender.allelic_z - allelic_frequency_new_neutral[4];
      //
      //       cout << "  parent_population_sender.allelic_a" <<   current_population_sender.allelic_a << endl;
      //       cout << "  parent_population_sender.allelic_b" <<   current_population_sender.allelic_b << endl;
      //       cout << "  parent_population_sender.allelic_c" <<   current_population_sender.allelic_c << endl;
      //       cout << "  parent_population_sender.allelic_d" <<   current_population_sender.allelic_d << endl;
      //       cout << "  parent_population_sender.allelic_e" <<   current_population_sender.allelic_e << endl;
      //

      current_population_sender.pop_size =  current_population_sender.pop_size - new_population.pop_size;
      //      cout << "current_population_sender.pop_size after: " << current_population_sender.pop_size << endl;


      double new_fitnesslike;
      new_fitnesslike = link_fitnesslike_mu_gamma(ready_to_colonize2[random_cell_to_go - 1], new_population, alleles_adaptation_coef, map1);
      computed_rate_based_on_temperature.push_back(new_fitnesslike);
      populations_this_species.push_back(new_population);
      populations_this_species[id_cells2[random_cell_to_go - 1]] = current_population_sender;
      // cout << "total_pop_size at this point " << total_pop_size << endl;
      // cout << "original_population_sender_pop_size " << original_population_sender_pop_size << endl;
      // cout << "current_population_sender.pop_size "<< current_population_sender.pop_size << endl;
      total_pop_size = total_pop_size - (original_population_sender_pop_size - (current_population_sender.pop_size + new_population.pop_size));
      //     cout << "total_pop_size new: " << total_pop_size << endl;
      range = range + 1;
      x_coordinate_last_event = ready_to_colonize2[random_cell_to_go - 1].x;
      y_coordinate_last_event = ready_to_colonize2[random_cell_to_go - 1].y;
      // classify_elevation(true,map1);

      // for(int i = 0; i < computed_rate_based_on_temperature.size(); ++i){
      //   cout << "computed_rate_based_on_temperature: " << computed_rate_based_on_temperature[i] << endl;
      // }

      // update the map
      int current_richness_map;
      current_richness_map = map1[(ready_to_colonize2[random_cell_to_go - 1].y - 1)][(ready_to_colonize2[random_cell_to_go - 1].x - 1)].Nspp;
      map1[(ready_to_colonize2[random_cell_to_go - 1].y - 1)][(ready_to_colonize2[random_cell_to_go - 1].x - 1)].Nspp = current_richness_map + 1;


      // for the new colonized cell
      map1[(ready_to_colonize2[random_cell_to_go - 1].y - 1)][(ready_to_colonize2[random_cell_to_go - 1].x - 1)].total_abundance_cell = abundance_individuals_cell_to_go + new_pop_size;
      //cout << "abundance of new colonized cell:" << map1[(ready_to_colonize2[random_cell_to_go - 1].y - 1)][(ready_to_colonize2[random_cell_to_go - 1].x - 1)].total_abundance_cell << endl;
      // for the senders cell
      int current_pop_map2;
      current_pop_map2 = map1[(presence[id_cells2[random_cell_to_go - 1]].y - 1)][(presence[id_cells2[random_cell_to_go - 1]].x - 1)].total_abundance_cell;
      map1[(presence[id_cells2[random_cell_to_go - 1]].y - 1)][(presence[id_cells2[random_cell_to_go - 1]].x - 1)].total_abundance_cell = (current_pop_map2 - original_population_sender_pop_size ) + (original_population_sender_pop_size - new_pop_size);

      // some checks
      int new_sum_from_allele;
      new_sum_from_allele = new_population.allelic_a + new_population.allelic_b + new_population.allelic_c + new_population.allelic_d + new_population.allelic_e;
      int new_sum_from_allele_neutral;
      new_sum_from_allele_neutral = new_population.allelic_v + new_population.allelic_w + new_population.allelic_x + new_population.allelic_y + new_population.allelic_z;


      int sender_sum_from_allele;
      sender_sum_from_allele = current_population_sender.allelic_a + current_population_sender.allelic_b + current_population_sender.allelic_c + current_population_sender.allelic_d + current_population_sender.allelic_e;
      int sender_sum_from_allele_neutral;
      sender_sum_from_allele_neutral = current_population_sender.allelic_v + current_population_sender.allelic_w + current_population_sender.allelic_x + current_population_sender.allelic_y + current_population_sender.allelic_z;


      if(new_sum_from_allele != new_sum_from_allele_neutral && sender_sum_from_allele!=sender_sum_from_allele_neutral)
      {
        stop(" +++++++++++++ Error: mistmatch between neutral and adaptive loci ++++++++++++++");

      }


      if((new_sum_from_allele + sender_sum_from_allele) != original_population_sender_pop_size)
      {
        stop(" +++++++++++++ Error: original populatio did not split correctly (alleles) ++++++++++++++");
      }

      if((new_population.pop_size + current_population_sender.pop_size) != pop_size_from_sum_alle)
      {
        stop(" +++++++++++++ Error: original populatio did not split correctly (pop_size) ++++++++++++++");
      }

      update_latitudinal_borders(t,true);

    }
    else
    {
      //cout << "this population had less than 3 individuals and expansion did not take place" << endl;
    }


  }
}

void happening_speciation(vector<species>& all_species, vector<int> alleles_adaptation_coef, int species_to_do, double t, double full_saturation_indi, landscape **map1)
{
  random_device rd;
  default_random_engine generator5(rd());
  species focal;
  focal = all_species[species_to_do];
  if(focal.presence.size() > 1 ){
    vector <int> prob_based_on_popsize;
    vector <int> cells_presence;
    for (int i = 0; i < focal.presence.size(); ++i)
    {
      cells_presence.push_back(i);
      prob_based_on_popsize.push_back(focal.populations_this_species[i].pop_size);
    }

    //  cout << "prob:" << endl;
    //  for(int i = 0; i < focal.presence.size(); ++i){
    //
    //    cout <<  prob_based_on_elevation[i] << "- ";
    //  }
    int random_cell_to_mutate;
    discrete_distribution<int> events_probabilities_to_pick(prob_based_on_popsize.begin(), prob_based_on_popsize.end());
    random_cell_to_mutate = cells_presence[events_probabilities_to_pick(generator5)];
    //  cout << "picked "<<random_cell_to_mutate << endl;
    //  cout <<" numb of pop parent before speciation: " << focal.populations_this_species.size() << endl;
    species new_species;
    new_species.initial_position(focal.presence[random_cell_to_mutate]);
    new_species.birthplace.x = focal.presence[random_cell_to_mutate].x;
    new_species.birthplace.y = focal.presence[random_cell_to_mutate].y;
    new_species.southernmost = focal.presence[random_cell_to_mutate].y;
    new_species.northernmost = focal.presence[random_cell_to_mutate].y;
    new_species.parent = focal.id;
    new_species.birth = t;
    new_species.saturation_grid_birth = full_saturation_indi;
    new_species.trait_state = focal.trait_state;
    new_species.alive = true;
    new_species.id = all_species.size() + 1;
    new_species.x_coordinate_last_event = focal.presence[random_cell_to_mutate].x;
    new_species.y_coordinate_last_event = focal.presence[random_cell_to_mutate].y;
    new_species.temperature_optimum.push_back(focal.temperature_optimum[random_cell_to_mutate]);
    // new_species.classify_elevation(true,map1);
    new_species.classify_elevation_origin(true);
    population_structure first_population_newspecies;
    first_population_newspecies.allelic_a = focal.populations_this_species[random_cell_to_mutate].allelic_a;
    first_population_newspecies.allelic_b = focal.populations_this_species[random_cell_to_mutate].allelic_b;
    first_population_newspecies.allelic_c = focal.populations_this_species[random_cell_to_mutate].allelic_c;
    first_population_newspecies.allelic_d = focal.populations_this_species[random_cell_to_mutate].allelic_d;
    first_population_newspecies.allelic_e = focal.populations_this_species[random_cell_to_mutate].allelic_e;

    first_population_newspecies.allelic_v = focal.populations_this_species[random_cell_to_mutate].allelic_v;
    first_population_newspecies.allelic_w = focal.populations_this_species[random_cell_to_mutate].allelic_w;
    first_population_newspecies.allelic_x = focal.populations_this_species[random_cell_to_mutate].allelic_x;
    first_population_newspecies.allelic_y = focal.populations_this_species[random_cell_to_mutate].allelic_y;
    first_population_newspecies.allelic_z = focal.populations_this_species[random_cell_to_mutate].allelic_z;
    first_population_newspecies.pop_size = focal.populations_this_species[random_cell_to_mutate].pop_size;
    new_species.populations_this_species.push_back(first_population_newspecies);
    new_species.total_pop_size = first_population_newspecies.pop_size;
    new_species.computed_rate_based_on_temperature.push_back(link_fitnesslike_mu_gamma(focal.presence[random_cell_to_mutate], first_population_newspecies, alleles_adaptation_coef, map1));
    all_species.push_back(new_species);
    // cout << "position new species:" << new_species.presence[0].x << "," << new_species.presence[0].y << endl;
    // cout << " +++++  cells presence before speciation:" << focal.presence.size() << endl;
    focal.x_coordinate_last_event = focal.presence[random_cell_to_mutate].x;
    focal.y_coordinate_last_event = focal.presence[random_cell_to_mutate].y;




    // focal.classify_elevation(false, map1);
    focal.presence.erase(focal.presence.begin() + random_cell_to_mutate);
    focal.computed_rate_based_on_temperature.erase(focal.computed_rate_based_on_temperature.begin() + random_cell_to_mutate);
    focal.populations_this_species.erase(focal.populations_this_species.begin() + random_cell_to_mutate);
    focal.temperature_optimum.erase(focal.temperature_optimum.begin() + random_cell_to_mutate);
    focal.range = focal.range - 1;

    // latitudinal limits have changed with a population turning into a different species?

    vector <int> all_ys_this_species;
    for(int i = 0; i < focal.presence.size(); i++)
    {
      all_ys_this_species.push_back(focal.presence[i].y);

    }
    focal.northernmost = *min_element(all_ys_this_species.begin(),all_ys_this_species.end());

    focal.southernmost = *max_element(all_ys_this_species.begin(),all_ys_this_species.end());

    //  cout <<" population size_parent before speciation: " << focal.total_pop_size << endl;

    focal.total_pop_size = focal.total_pop_size - first_population_newspecies.pop_size;

    // cout <<" population size_parent after speciation: " << focal.total_pop_size << endl;
    // cout <<" numb of pop parent after speciation: " << focal.populations_this_species.size() << endl;
    // cout <<" population size_child: " << first_population_newspecies.pop_size << endl;
    //
    //   cout << "first_population_newspecies.allelic_a" << first_population_newspecies.allelic_a << endl;
    //   cout << "first_population_newspecies.allelic_b" << first_population_newspecies.allelic_b << endl;
    //   cout << "first_population_newspecies.allelic_c" << first_population_newspecies.allelic_c << endl;
    //   cout << "first_population_newspecies.allelic_d" << first_population_newspecies.allelic_d << endl;
    //   cout << "first_population_newspecies.allelic_e" << first_population_newspecies.allelic_e << endl;

    //  int its_elevation;
    //  its_elevation = focal.presence[random_cell_to_mutate].x;
    //  if(its_elevation <= lower_limit_highlands){
    //    new_species.range_highlands =  1;
    //    focal.range_highlands = focal.range_highlands - 1;
    //  } else {
    //    if(its_elevation <= lower_limit_intermediate1){
    //      new_species.range_intermediate1 =  1;
    //      focal.range_highlands = focal.range_intermediate2 - 1;
    //    } else {
    //      if(its_elevation <= lower_limit_intermediate2){
    //        new_species.range_intermediate2 =  1;
    //        focal.range_intermediate2 = focal.range_intermediate2 - 1;
    //      } else{
    //        new_species.range_lowlands =  1;
    //        focal.range_lowlands = focal.range_lowlands - 1;
    //      }
    //    }
    //  }
    // cout << " +++++  cells presence After speciation:" << focal.presence.size() << endl;
    if (focal.range == 0)
    { // it died giving birth
      focal.death = t;
      focal.alive = false;
    }
    all_species[species_to_do] = focal; // this line updates the all_species vector
    all_species[species_to_do] = all_species[species_to_do]; // this line updates the all_species vector
    // cout << "richness in here: " << all_species.size() << endl;
    //show_all_species_data(all_species);
  }

}

void set_landscape(IntegerVector map_elevation_vector, IntegerVector map_k_vector, IntegerVector map_temperature_vector, int y_max, int x_max, landscape **map1)
{

  // landscape map1[x_max][y_max];
  // mapp.open(name1.c_str());

  for (int i = 0; i < y_max; i++)
  {
    for (int j = 0; j < x_max; j++)
    {
      //  cout << map_k_vector[(i*y_max) + j] << "\t";
      map1[i][j].habitat = map_k_vector[(i * y_max) + j];
      map1[i][j].k_patch = map_k_vector[(i * y_max) + j];
      map1[i][j].temperature = map_temperature_vector[(i * y_max) + j];
      map1[i][j].elevation = map_elevation_vector[(i * y_max) + j];
      map1[i][j].total_abundance_cell = 0;
    }
    // cout << endl;
  }

  /*
   for (int i = 0; i < y_max; i++) {
   for (int j = 0; j < x_max; j++) {
   cout << map1[i][j].habitat << "\t";
   }
   cout << endl;
   }
   */

  // to "clean" the map

  for (int i = 0; i < y_max; i++)
  {
    for (int j = 0; j < x_max; j++)
    {
      map1[i][j].Nspp = 0;
    }
  }

  for (int i = 0; i < y_max; i++)
  {
    for (int j = 0; j < x_max; j++)
    {
      map1[i][j].sum_trait = 0;
    }
  }

  for (int i = 0; i < y_max; i++)
  {
    for (int j = 0; j < x_max; j++)
    {
      map1[i][j].average_trait = 0;
    }
  }

  /*
   mapp.open(name2.c_str());

   for (int i = 0; i < x_max; i++) {
   for (int j = 0; j < y_max; j++) {
   mapp >> map2[i][j].habitat;
   }
   }
   mapp.close();
   */
  // return map1;
}

// functions for Rcpp
vector<int> extract_yx_func(species process_one_species, std::string either_yx)
{
  vector<int> extract_yx;
  for (int i = 0; i < process_one_species.presence.size(); ++i)
  {
    if (either_yx == "x")
    {
      //  cout << "show the x " << process_one_species.presence[i].x << endl;
      extract_yx.push_back(process_one_species.presence[i].x);
    }
    if (either_yx == "y")
    {
      //    cout << "show the y " << process_one_species.presence[i].y << endl;
      extract_yx.push_back(process_one_species.presence[i].y);
    }
  }
  return extract_yx;
}

vector<double> extract_temperature_func(species process_one_species)
{
  vector<double> extract_temperature;
  for (int i = 0; i < process_one_species.presence.size(); ++i)
  {
    extract_temperature.push_back(process_one_species.temperature_optimum[i]);
  }
  return extract_temperature;
}



List get_me_output (vector<species> all_species, double t)
{
  List list_all_species_distribution = List::create();
  List list_all_species_allele_A = List::create();
  List list_all_species_allele_B = List::create();
  List list_all_species_allele_C = List::create();
  List list_all_species_allele_D = List::create();
  List list_all_species_allele_E = List::create();

  List list_all_species_allele_V = List::create();
  List list_all_species_allele_W = List::create();
  List list_all_species_allele_X = List::create();
  List list_all_species_allele_Y = List::create();
  List list_all_species_allele_Z = List::create();


  List all_change_northernmost = List::create();
  List all_change_southernmost = List::create();
  List all_time_change_northernmost = List::create();
  List all_time_change_southernmost = List::create();

  List list_all_species_popsize_perPop = List::create();

  List list_time = List::create();
  for (int iij = 0; iij < all_species.size(); ++iij)
  {
    List xy_one_species = List::create();
    List popsize_one_species = List::create();
    List list_allele_A_one_species = List::create();
    List list_allele_B_one_species = List::create();
    List list_allele_C_one_species = List::create();
    List list_allele_D_one_species = List::create();
    List list_allele_E_one_species = List::create();

    List list_allele_V_one_species = List::create();
    List list_allele_W_one_species = List::create();
    List list_allele_X_one_species = List::create();
    List list_allele_Y_one_species = List::create();
    List list_allele_Z_one_species = List::create();
    List list_popsize_perPop =  List::create();
    // cout << " seeing the output: " << iij << endl;
    species take_one_ouput;
    take_one_ouput = all_species[iij];

    if (take_one_ouput.alive)
    {
      vector <int> allele_A_one_species;
      vector <int> allele_B_one_species;
      vector <int> allele_C_one_species;
      vector <int> allele_D_one_species;
      vector <int> allele_E_one_species;

      vector <int> allele_V_one_species;
      vector <int> allele_W_one_species;
      vector <int> allele_X_one_species;
      vector <int> allele_Y_one_species;
      vector <int> allele_Z_one_species;
      vector <int> popsize_perPop;
      for(int ii = 0; ii < take_one_ouput.populations_this_species.size(); ++ii)
      {


        allele_A_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_a);
        allele_B_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_b);
        allele_C_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_c);
        allele_D_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_d);
        allele_E_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_e);

        allele_V_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_v);
        allele_W_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_w);
        allele_X_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_x);
        allele_Y_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_y);
        allele_Z_one_species.push_back(take_one_ouput.populations_this_species[ii].allelic_z);
        popsize_perPop.push_back(take_one_ouput.populations_this_species[ii].pop_size);
      }

      vector<int> one_species_x;
      vector<int> one_species_y;

      vector<double> vector_popsize;
      vector_popsize = extract_temperature_func(take_one_ouput);

      list_allele_A_one_species.push_back(allele_A_one_species);
      list_allele_B_one_species.push_back(allele_B_one_species);
      list_allele_C_one_species.push_back(allele_C_one_species);
      list_allele_D_one_species.push_back(allele_D_one_species);
      list_allele_E_one_species.push_back(allele_E_one_species);

      list_allele_V_one_species.push_back(allele_V_one_species);
      list_allele_W_one_species.push_back(allele_W_one_species);
      list_allele_X_one_species.push_back(allele_X_one_species);
      list_allele_Y_one_species.push_back(allele_Y_one_species);
      list_allele_Z_one_species.push_back(allele_Z_one_species);

      list_popsize_perPop.push_back(popsize_perPop);

      one_species_x = extract_yx_func(take_one_ouput, "x");
      one_species_y = extract_yx_func(take_one_ouput, "y");
      xy_one_species.push_back(one_species_x);
      xy_one_species.push_back(one_species_y);
      popsize_one_species.push_back(vector_popsize);
    }
    else
    {
      std::string dead_species;
      dead_species = "extinct";
      xy_one_species = dead_species;
      popsize_one_species = dead_species;
      // allele_A_one_species = dead_species;
      // allele_B_one_species = dead_species;
      // popsize_perPop = dead_species;
      // cout << "this went extinct" << endl;
    }
    list_all_species_allele_A.push_back(list_allele_A_one_species);
    list_all_species_allele_B.push_back(list_allele_B_one_species);
    list_all_species_allele_C.push_back(list_allele_C_one_species);
    list_all_species_allele_D.push_back(list_allele_D_one_species);
    list_all_species_allele_E.push_back(list_allele_E_one_species);

    list_all_species_allele_V.push_back(list_allele_V_one_species);
    list_all_species_allele_W.push_back(list_allele_W_one_species);
    list_all_species_allele_X.push_back(list_allele_X_one_species);
    list_all_species_allele_Y.push_back(list_allele_Y_one_species);
    list_all_species_allele_Z.push_back(list_allele_Z_one_species);
    list_all_species_popsize_perPop.push_back(list_popsize_perPop);

    list_all_species_distribution.push_back(xy_one_species);
    all_change_northernmost.push_back(take_one_ouput.change_northernmost);
    all_change_southernmost.push_back(take_one_ouput.change_southernmost);
    all_time_change_northernmost.push_back(take_one_ouput.time_change_northernmost);
    all_time_change_southernmost.push_back(take_one_ouput.time_change_southernmost);

  }
  list_time.push_back(t);
  List list_extract_species_data = List::create();
  list_extract_species_data = extract_species_data(all_species);
  //}

  List model_output = List::create(Named("Distribution") = list_all_species_distribution,
                                   _["popsize_perPop"] = list_all_species_popsize_perPop,
                                   _["all_change_northernmost"] = all_change_northernmost,
                                   _["all_change_southernmost"] = all_change_southernmost,
                                   _["all_time_change_northernmost"] = all_time_change_northernmost,
                                   _["all_time_change_southernmost"] = all_time_change_southernmost,
                                   _["Allele_A"] = list_all_species_allele_A,
                                   _["Allele_B"] = list_all_species_allele_B,
                                   _["Allele_C"] = list_all_species_allele_C,
                                   _["Allele_D"] = list_all_species_allele_D,
                                   _["Allele_E"] = list_all_species_allele_E,
                                   _["Allele_V"] = list_all_species_allele_V,
                                   _["Allele_W"] = list_all_species_allele_W,
                                   _["Allele_X"] = list_all_species_allele_X,
                                   _["Allele_Y"] = list_all_species_allele_Y,
                                   _["Allele_Z"] = list_all_species_allele_Z,
                                   _["Species_Info"] = list_extract_species_data,
                                   _["total_time"] = list_time);
  return model_output;
}







List extract_species_data(vector<species> process_all_species)
{
  vector<int> extract_id;
  vector<int> extract_parent;
  vector<double> extract_birth;
  vector <int> extract_birthplace_x;
  vector <int> extract_birthplace_y;
  vector <double> extract_saturation_grid_birth;
  vector<double> extract_death;
  vector<double> extract_trait;
  vector<int> extract_range;
  vector<int> extract_popsize;
  vector <int> extract_succefulgeneflow;
  vector <int> extract_southernmost;
  vector <int> extract_northernmost;




  for (int i = 0; i < process_all_species.size(); ++i)
  {
    species process_one_species = process_all_species[i];
    extract_id.push_back(process_one_species.id);
    extract_parent.push_back(process_one_species.parent);
    extract_birth.push_back(process_one_species.birth);
    extract_birthplace_x.push_back(process_one_species.birthplace.x);
    extract_birthplace_y.push_back(process_one_species.birthplace.y);
    extract_death.push_back(process_one_species.death);
    extract_trait.push_back(process_one_species.trait_state);
    extract_range.push_back(process_one_species.range);
    extract_popsize.push_back(process_one_species.total_pop_size);
    extract_succefulgeneflow.push_back(process_one_species.succesful_geneflow_events);
    extract_southernmost.push_back(process_one_species.southernmost);
    extract_northernmost.push_back(process_one_species.northernmost);
    // cout << "process_one_species.saturation_grid_birth " << process_one_species.saturation_grid_birth << endl;
    extract_saturation_grid_birth.push_back(process_one_species.saturation_grid_birth);
  }
  List list_extract_species_data = List::create(Named("ID") = extract_id,
                                                _["Parent"] = extract_parent,
                                                _["Birth"] = extract_birth,
                                                _["BirthPlaceX"] = extract_birthplace_x,
                                                _["BirthPlaceY"] = extract_birthplace_y,
                                                _["Death"] = extract_death,
                                                _["TraitValue"] = extract_trait,
                                                _["RangeSize"] = extract_range,
                                                _["populationSize"] = extract_popsize,
                                                _["saturation_grid_birth"] = extract_saturation_grid_birth,
                                                _["geneflow_events"] =  extract_succefulgeneflow,
                                                _["northernnmost_locat"] =   extract_northernmost,
                                                _["southernmost_locat"] =   extract_southernmost);


  return list_extract_species_data;
}
// end for Rcpp

vector<species> get_species_intocpp(vector<species> all_species, IntegerVector all_alleles,IntegerVector all_alleles_neutral, IntegerVector all_popsize, IntegerVector all_x, IntegerVector all_y, IntegerVector all_IDs, IntegerVector all_parents, NumericVector all_births, NumericVector all_deaths, NumericVector all_traits, IntegerVector all_ranges, int number_spp, landscape **map1, vector<int> alleles_adaptation_coef, double v, double gamma, double mu, std::string temperature_influencing)
{

  // vector <species> all_species;
  int counter_species_range = 0;
  for (int i = 0; i < number_spp; ++i)
  {
    species put_one_species;
    put_one_species.id = all_IDs[i];
    put_one_species.parent = all_parents[i];
    put_one_species.range = all_ranges[i];
    put_one_species.birth = all_births[i];
    put_one_species.death = all_deaths[i];
    put_one_species.x_coordinate_last_event = all_x[i];
    put_one_species.y_coordinate_last_event = all_y[i];
    put_one_species.temperature_optimum.push_back(all_traits[i]);
    put_one_species.total_pop_size = all_popsize[i];
    put_one_species.birthplace.x = all_x[i];
    put_one_species.birthplace.y = all_y[i];
    put_one_species.northernmost = all_x[i];
    put_one_species.southernmost = all_x[i];
    put_one_species.saturation_grid_birth = 0;
    if (put_one_species.range == 0)
    {
    }
    else
    {
      for (int ii = 0; ii < put_one_species.range; ++ii)
      {

        yx this_is_one_yx;
        double fitnesslike_thispopulation;

        this_is_one_yx.x = all_x[counter_species_range + ii];
        // cout <<  all_x[counter_species_range + ii ] << endl;
        // cout <<  counter_species_range + ii   << endl;
        this_is_one_yx.y = all_y[counter_species_range + ii];
        put_one_species.presence.push_back(this_is_one_yx);

        population_structure this_pop;
        this_pop.allelic_a = all_alleles[0];
        this_pop.allelic_b = all_alleles[1];
        this_pop.allelic_c = all_alleles[2];
        this_pop.allelic_d = all_alleles[3];
        this_pop.allelic_e = all_alleles[4];

        this_pop.allelic_v = all_alleles_neutral[0];
        this_pop.allelic_w = all_alleles_neutral[1];
        this_pop.allelic_x = all_alleles_neutral[2];
        this_pop.allelic_y = all_alleles_neutral[3];
        this_pop.allelic_z = all_alleles_neutral[4];
        this_pop.pop_size = all_popsize[i];
        fitnesslike_thispopulation = link_fitnesslike_mu_gamma(this_is_one_yx, this_pop, alleles_adaptation_coef, map1);
        put_one_species.populations_this_species.push_back(this_pop);
        put_one_species.computed_rate_based_on_temperature.push_back(fitnesslike_thispopulation);
      }
    }
    put_one_species.trait_state = all_traits[i];
    if (put_one_species.death == 0)
    {
      put_one_species.alive = true;
    }
    else
    {
      put_one_species.alive = false;
    }
    // put_one_species.classify_elevation(true,map1);
    put_one_species.classify_elevation_origin(true);
    counter_species_range = counter_species_range + put_one_species.range;
    all_species.push_back(put_one_species);
  }
  for (int i = 0; i < number_spp; ++i)
  {
    // cout << "birth of " << all_species[i].id << "at " << all_species[i].birth << endl;
    if (all_species[i].alive)
    {
      // for (int ii = 0; ii < all_species[i].range; ++ii)
      // {
      //   //  cout << "x: " << all_species[i].presence[ii].x <<" y: " << all_species[i].presence[ii].y  << endl;
      // }
    }
  }
  return (all_species);
}
