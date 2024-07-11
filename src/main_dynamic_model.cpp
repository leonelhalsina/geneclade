#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include "dynamic_model.h"
// #include "dynamic_model.cpp"

using namespace Rcpp;
using namespace std;

//default_random_engine &generator;
//' @export
 // [[Rcpp::export]]
 List do_simulation(IntegerVector map_elevation_vector, IntegerVector map_k_vector, IntegerVector map_temperature_vector, IntegerVector map_temperature_vector2,bool extirpation_depen_temperature,
                    bool colonization_depen_temperature, int x_max, int y_max, IntegerVector all_x, IntegerVector all_y,
                    IntegerVector all_IDs, IntegerVector all_parents, NumericVector all_births, NumericVector all_deaths,
                    NumericVector all_traits, IntegerVector all_ranges, IntegerVector all_alleles, IntegerVector all_alleles_neutral,
                    IntegerVector all_popsize, int number_spp, int the_seed, double mutation_rate ,double percentage_flow,
                    double geneflow_rate, double popchange_rate,NumericVector the_gammas, NumericVector the_mus,
                    double q, double lambda, bool species_trait_state_gamma,double sd_normal_distribution_traitevol, double mean_normal_distribution_traitevol,
                    double sd_normal_distribution_pop_change,double starting_time,double simulated_time, int max_spp,int maximum_cycles, bool use_k, double restiction_par, std::string show_richness_map,
                    double v, IntegerVector alleles_adaptation_coef2,bool do_change_map_rates, bool vicariant_speciation,bool speciation_rangesize_unlinked, bool colonization_rangesize_unlinked, double time_percent_stop_after_first_equilibrium_and_disturbance)
 {

   if((x_max * y_max ) != map_k_vector.size()){
     stop("map size does not match y_max and/or x_max");
   }

   random_device rd;
   default_random_engine generator(rd());
   double gamma;
   double second_gamma;
   gamma = the_gammas[0];
   second_gamma = the_gammas[1];

   double mu;
   double second_mu;
   mu = the_mus[0];
   second_mu = the_mus[1];

   // if (all_ranges.size() > 1)
   // {
   //
   //   stop("initialization of simulation should be with one population only, or go and work on get_species_intocpp function" );
   // }
   int total_pop_from_allelevector;
   total_pop_from_allelevector = 0;
   for(int ij = 0; ij < all_alleles.size(); ++ij){
     total_pop_from_allelevector = total_pop_from_allelevector + all_alleles[ij];
   }

   if(total_pop_from_allelevector != all_popsize[0]){
     cout << all_popsize[0] << endl;
     cout << total_pop_from_allelevector << endl;
     stop("check your all_alleles and all_popsize, they should match");
   }

   vector<species> all_species;
   landscape **map1 = new landscape *[y_max];

   // dynamically allocate memory of size N for each row
   for (int i = 0; i < y_max; ++i)
   {
     map1[i] = new landscape[x_max];
   }

   set_landscape(map_elevation_vector, map_k_vector, map_temperature_vector, y_max, x_max, map1);

   int absolute_potential_abundance = 0;
   int absolute_cells_to_live = 0;
   for (int ij = 0; ij < map_k_vector.size(); ++ij)
   {
     if(map_k_vector[ij] != -9)
     {
       absolute_potential_abundance = absolute_potential_abundance + map_k_vector[ij];
       absolute_cells_to_live = absolute_cells_to_live + 1;
     }
   }
   // if (the_seed != 0)
   // {
   //   generator.seed(the_seed);
   // }




   vector<int> alleles_adaptation_coef;
   for (int ij = 0; ij < alleles_adaptation_coef2.size(); ++ij)
   {
     alleles_adaptation_coef.push_back(alleles_adaptation_coef2[ij]);
   }


   all_species = get_species_intocpp(all_species,  all_alleles, all_alleles_neutral, all_popsize,  all_x,  all_y,  all_IDs,  all_parents,  all_births,  all_deaths,  all_traits,  all_ranges,  number_spp, map1, alleles_adaptation_coef,  v,  gamma,  mu, extirpation_depen_temperature,  colonization_depen_temperature);

   int total_num_populations = 0;
   for (int ij = 0; ij < all_ranges.size(); ++ij)
   {
     cout << all_ranges[ij] << endl;
     total_num_populations = total_num_populations + all_ranges[ij];
   }

   //  landscape map1[][];
   vector<string> vector_events_tookplace;
   int total_expansion_events = 0;
   int total_contraction_events (0);
   int total_mutation_events = 0;
   int attempted_geneflow_events = 0;
   int total_speciation_events = 0;
   int total_geneflow_events = 0;
   int total_popchange_events = 0;
   int total_traitevolution_events = 0;
   std::string event_to_do;

   int full_saturation_indi;
   int total_indviduals;

   double first_equi_at;
   int richness_at_equilibrium;
   int populations_at_equilibrium;
   int individuals_at_equilibrium;
   double t;
   t = starting_time;
   bool pending_change_in_rates;
   bool equilibrium_achieved;
   equilibrium_achieved = false;
   pending_change_in_rates = true;

   vector <int> accumul_satu_values;
   int previous_satu = 0;
   int count_same_satura = 0;
   vector <int> how_long_took_to_change_percentage;
   vector <int> show_the_percentages_saturation;
   int counting_to_print_variance = 0;

   int cycles = 0;
   // for (int cycles = 0; cycles < maximum_cycles; ++cycles) {
   //
   set_landscape(map_elevation_vector, map_k_vector, map_temperature_vector, y_max, x_max, map1);
   populate_landscape(y_max,x_max,all_species, map1);
   to_show_richness_map(y_max,x_max,show_richness_map,all_species,map1);
   while (t < simulated_time && total_num_populations > 0 && cycles < maximum_cycles)
   {

     cycles = cycles + 1;
     //show_all_species_data(all_species);
     // to pick a species

     // vector<int> population_sizes;

     // int total_num_populations;
     // total_num_populations = 0;

     total_indviduals = 0;
     int total_num_populations = 0;
     vector<int> id_alive_species;
     for (int i = 0; i < all_species.size(); ++i)
     {
       if (all_species[i].alive)
       {
         id_alive_species.push_back(all_species[i].id);
         //population_sizes.push_back(all_species[i].total_pop_size);
         total_num_populations = total_num_populations + all_species[i].range;
         total_indviduals = total_indviduals + all_species[i].total_pop_size;
         //     cout << "range size, species" << all_species[i].id << " is " << all_species[i].range << endl;
       }
     }


     int species_to_do;
     // species_to_do = id_alive_species[give_me_random_uniform(0, (id_alive_species.size()-1))];

     // the trait state determines all the probabilities.
     vector<double> total_probability_species;
     probabilities_based_traits calculation_probabilities;

     calculation_probabilities = calculate_probabilities_using_traitstate(all_species, map1,  extirpation_depen_temperature, colonization_depen_temperature, species_trait_state_gamma, mutation_rate, geneflow_rate, popchange_rate, lambda, gamma, mu, id_alive_species, v);
     total_probability_species = calculation_probabilities.total_probability_species;
     discrete_distribution<int> species_probabilities_to_pick(total_probability_species.begin(), total_probability_species.end());
     species_to_do = id_alive_species[species_probabilities_to_pick(generator)];

     // Gillispie algorithm

     double time_elapsed;
     exponential_distribution<double> waiting_times(calculation_probabilities.total_waiting_times_rates);

     // cout << " this is the TOTAL rate for Gillispie: " << calculation_probabilities.total_waiting_times_rates << endl;

     time_elapsed = waiting_times(generator);
     double t_previous_cycle;
     t_previous_cycle = t;
     t = t + time_elapsed;
     // cout << "total abundance: " << total_num_populations << "..and computed from elevation info:" << (populations_highlands +populations_intermediate1 +populations_intermediate2 + populations_lowlands) << endl;

     full_saturation_indi = round(((total_indviduals/(double)absolute_potential_abundance) * 100.0 ));


     full_saturation_indi = round((((total_num_populations/id_alive_species.size())/(double)absolute_cells_to_live) * 100.0 ));

     // bit that looks for equilibrium and stops simulation when certain conditions are met






     // show_the_percentages_saturation.push_back(full_saturation_indi);
     //
     // if(event_to_do == "speciation" ||event_to_do == "contraction" || event_to_do == "expansion" ){
     //   if(full_saturation_indi == previous_satu)
     //   {
     //     count_same_satura = count_same_satura + 1;
     //   }
     //   else
     //
     //   {
     //     how_long_took_to_change_percentage.push_back(count_same_satura);
     //     count_same_satura = 0;
     //   }
     //   //      if(full_saturation_indi > previous_satu)
     //   //      {
     //   //
     //   // count_increase_satura = count_increase_satura + 1;
     //   //      }
     //   //      if(full_saturation_indi < previous_satu)
     //   //      {
     //   //
     //   //        count_decrease_satura = count_decrease_satura + 1;
     //   //      }
     //   previous_satu = full_saturation_indi;
     //
     // }





     if(event_to_do == "speciation" ||event_to_do == "contraction" || event_to_do == "expansion" ){
       //
       //      cout << "this full_saturation_indi " << full_saturation_indi << " at cyle: " << cycles << endl;
       accumul_satu_values.push_back(full_saturation_indi);
       //     cout << "lenghth accumul_satu_values: " << accumul_satu_values.size() << endl;
       //      for (int ij = 0; ij < accumul_satu_values.size(); ++ij){
       //        cout << " "<< accumul_satu_values[ij] ;
       //      }
       // cout << endl;
       if(accumul_satu_values.size() > 3000)
       {
         cout << " over limit in vector length" << endl;
         break;
       }
       if(accumul_satu_values.size() == 3000)
       {
         double this_variance;
         this_variance = calculate_variance(accumul_satu_values);
         counting_to_print_variance = counting_to_print_variance + 1;
         if(counting_to_print_variance == 100 )
         {
           cout << "Variance in average range size: " << this_variance << endl;
           counting_to_print_variance = 0;
         }

         if(this_variance == 0){
           cout << "equilibrium found: " << endl;
           equilibrium_achieved = true;
           if(do_change_map_rates == false ){

             break;
           }

           if(do_change_map_rates  && pending_change_in_rates == false)
           {
             cout << "second equilibrium found" << endl;
             break;
           }

         }
         else
         {
           accumul_satu_values.erase(accumul_satu_values.begin());
         }
       }
     }

     if(equilibrium_achieved && do_change_map_rates && pending_change_in_rates)
     {
       change_temperature_map(x_max,y_max,map_temperature_vector2,map1);
       gamma = second_gamma;
       mu =  second_mu;
       if(time_percent_stop_after_first_equilibrium_and_disturbance != 0)
       {

         first_equi_at = t;
         richness_at_equilibrium = id_alive_species.size();
         populations_at_equilibrium = total_num_populations;
         individuals_at_equilibrium = total_indviduals;
         simulated_time = t + (t * (time_percent_stop_after_first_equilibrium_and_disturbance/100)); // to stop the simulation some time (some % of time) after the change in map/rates.
       }

       pending_change_in_rates = false;
       accumul_satu_values.clear();
       //cout << "size_ accumul_satu_values.clear() "<<  accumul_satu_values.size() << endl;
       cout << "___changing local extirpation rates and/or map temperature__" << endl;

     }

     // if(equilibrium_achieved && do_change_map_rates && pending_change_in_rates == false  && full_saturation_indi == stop_at_saturation)
     // {
     //   cout << "---- Equilibrium was reached, then the changes of mu or map temperature were changed, saturation decreased until stop_at_saturation " << endl;
     //   break;
     // }
     //
     // if(do_change_map_rates == false && full_saturation_indi == stop_at_saturation){
     //
     //   cout << "--- saturation reached the required level, no change in rates/maps nor equilibrium was met" << endl;
     //   break;
     // }


     if(id_alive_species.size() >= max_spp)
     {
       lambda = 0;
     }




     // end of it
     //if (t >= simulated_time || id_alive_species.size() >= max_spp)
     if (t >= simulated_time )
     {
       cout << "time: " << t << " cycle: " << cycles << " richness:" << id_alive_species.size() <<  " populations: " << total_num_populations << " indviduals: " << total_indviduals<<  " ind_saturation %: " << full_saturation_indi << endl;
       // cout << "total abundance: " << total_num_populations << "..and computed from elevation info:" << (populations_highlands +populations_intermediate1 +populations_intermediate2 + populations_lowlands) << endl;
       cout << "_________time is up" << endl;
       break;
     }
     //if ((round(t) - round(t_previous_cycle)) > 0)
     // cout << "t: " << t << " "<< t_previous_cycle << endl;
     // cout << (t + (t * 0.0005))  << endl;
     if (   (t_previous_cycle + (t_previous_cycle * 0.00001)) < t )
     {
       cout << "time: " << t << " cycle: " << cycles << " richness:" << id_alive_species.size() <<  " populations: " << total_num_populations << " indviduals: " << total_indviduals<< " ind_saturation %: " << full_saturation_indi << endl;      // cout << "total abundance: " << total_num_populations << "..and computed from elevation info:" << (populations_highlands +populations_intermediate1 +populations_intermediate2 + populations_lowlands) << endl;
     }

     // Here, the give_me_random function will pick a position of the id_alive_species vector
     // that Id, will be the element of all_species vector, and it will match its ID.
     // Imagine gime_me_random gives a 2, and this 2 element points to the id 8 inside
     // id_alive species (which can be 1,2,8,9). Then, that "8" will correspond to the 8th
     // element in all_species, which has id of 8.

     // cout << "			THE id of species doing something: " << species_to_do << endl;

     if (all_species[species_to_do - 1].id != species_to_do)
     {
       cout << "Error: it should be the same id" << endl;
       return 0;
     }

     species_to_do = species_to_do - 1; // because it is an index



     // cout << "calculation_probabilities.gammas_total " << calculation_probabilities.gammas_total  << endl;
     // cout << "calculation_probabilities.geneflow_rate_total " << calculation_probabilities.geneflow_rate_total  << endl;
     // cout << "calculation_probabilities.lambdas_total " << calculation_probabilities.lambdas_total  << endl;
     // cout << "calculation_probabilities.mus_total " << calculation_probabilities.mus_total  << endl;
     // cout << "calculation_probabilities.popchange_rate_total " << calculation_probabilities.popchange_rate_total  << endl;
     // cout << "calculation_probabilities.mutation_rate_total " << calculation_probabilities.mutation_rate_total  << endl;

     double total_rate_traitevol;
     //   I take the total rate to do trait evolution
     total_rate_traitevol = ((calculation_probabilities.popchange_rate_total * q)/100.0);
     // to pick and event
     discrete_distribution <int> events_probabilities_to_pick({calculation_probabilities.gammas_total, calculation_probabilities.mus_total, calculation_probabilities.popchange_rate_total, calculation_probabilities.lambdas_total,calculation_probabilities.geneflow_rate_total,calculation_probabilities.mutation_rate_total,total_rate_traitevol});
     // cout << "gammas " << calculation_probabilities.gammas_total << "mus " << calculation_probabilities.mus_total << "qs "<<  calculation_probabilities.qs_total << "lambdas " << calculation_probabilities.lambdas_total << endl;

     std::string all_events[7];
     all_events[0] = "expansion";
     all_events[1] = "contraction";
     all_events[2] = "pop_change";
     all_events[3] = "speciation";
     all_events[4] = "gene_flow";
     all_events[5] =  "mutation";
     all_events[6] =  "trait_evolution";


     event_to_do = all_events[events_probabilities_to_pick(generator)];

     vector <std::string> list_events_to_do;

     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("contraction");
     list_events_to_do.push_back("contraction");
     list_events_to_do.push_back("contraction");
     list_events_to_do.push_back("contraction");
     list_events_to_do.push_back("speciation");


     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("expansion");
     list_events_to_do.push_back("contraction");
     list_events_to_do.push_back("contraction");
     list_events_to_do.push_back("contraction");
     list_events_to_do.push_back("contraction");
     list_events_to_do.push_back("speciation");
     list_events_to_do.push_back("expansion");


     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("gene_flow");
     // list_events_to_do.push_back("speciation");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("contraction");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("gene_flow");
     // list_events_to_do.push_back("gene_flow");
     // list_events_to_do.push_back("gene_flow");
     // list_events_to_do.push_back("pop_change");
     // list_events_to_do.push_back("pop_change");
     // list_events_to_do.push_back("pop_change");
     // list_events_to_do.push_back("pop_change");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     //
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("gene_flow");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("gene_flow");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("pop_change");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("contraction");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("pop_change");
     // list_events_to_do.push_back("gene_flow");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("pop_change");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("mutation");
     // list_events_to_do.push_back("expansion");
     // list_events_to_do.push_back("gene_flow");
     // list_events_to_do.push_back("contraction");
     // list_events_to_do.push_back("pop_change");
     // list_events_to_do.push_back("speciation");
     // list_events_to_do.push_back("mutation");

     if(cycles < list_events_to_do.size()){
       //event_to_do = list_events_to_do[cycles - 1];   // to DELETE
     }
     // cout << "                       event_to_do: " << event_to_do << endl;


     // all_species[species_to_do].find_patches_distribution(); // to DELETE




     if (event_to_do == "expansion")
     {

       if(colonization_rangesize_unlinked)
       {
         int random_spp_speciation;
         species_to_do = give_me_random_uniform(0, id_alive_species.size() - 1);
         //stop("stop here");
       }
       // cout << "                  i will expand" << endl;

       //  cout << "                   species  BEFORE expansion: " << all_species[species_to_do].presence.size() << endl;
       all_species[species_to_do].happening_expansion(x_max, y_max, use_k, restiction_par, map1, colonization_depen_temperature, alleles_adaptation_coef,t); // restriction par will be either k or trait dissimilarity;
       //  cout << "                    species  AFTER expansion: " << all_species[species_to_do].presence.size() << endl;

       all_species[species_to_do] = all_species[species_to_do]; // this line updates the all_species vector
       total_expansion_events = total_expansion_events + 1;
     }
     if (event_to_do == "speciation")
     {
       //  cout << "                   i will speciate" << endl;

       if(speciation_rangesize_unlinked)
       {
         int random_spp_speciation;
         species_to_do = give_me_random_uniform(0, id_alive_species.size() - 1);
         //stop("stop here");
       }
       int rich_in_here;
       rich_in_here = all_species.size();

       if(vicariant_speciation){
         vector <contiguous_patches> list_patches;
         list_patches = all_species[species_to_do].find_patches_distribution(y_max,x_max);
         if(all_species[species_to_do].range > 2){
           happening_speciation(y_max,x_max, all_species, alleles_adaptation_coef, species_to_do, t, full_saturation_indi,map1, vicariant_speciation);
           all_species[species_to_do] = all_species[species_to_do]; // this line updates the all_species vector
           total_speciation_events = total_speciation_events + 1;


           if(all_species.size() != (rich_in_here + 1)){

             stop("speciation did not work out");
           }
         }


       } else {
         happening_speciation(y_max,x_max, all_species, alleles_adaptation_coef, species_to_do, t, full_saturation_indi,map1, vicariant_speciation);
         all_species[species_to_do] = all_species[species_to_do]; // this line updates the all_species vector
         total_speciation_events = total_speciation_events + 1;
       }




     }

     if (event_to_do == "trait_evolution"){

       all_species[species_to_do].happening_trait_evolution(mean_normal_distribution_traitevol,sd_normal_distribution_traitevol);
       all_species[species_to_do] = all_species[species_to_do]; // this line updates the all_species vector
       total_traitevolution_events = total_traitevolution_events + 1;
     }
     if (event_to_do == "gene_flow")
     {
       // cout << "                   i will gene_flow" << endl;

       all_species[species_to_do].happening_gene_flow(y_max,x_max,percentage_flow,map1);

       all_species[species_to_do] = all_species[species_to_do]; // this line updates the all_species vector
       attempted_geneflow_events = attempted_geneflow_events + 1;
     }

     if (event_to_do == "mutation")
     {
       // cout << "                   i will mutation" << endl;

       all_species[species_to_do].happening_mutation_this_species();
       all_species[species_to_do] = all_species[species_to_do]; // this line updates the all_species vector
       total_mutation_events = total_mutation_events + 1;
     }

     if (event_to_do == "pop_change")
     {
       // cout << "                  i will pop_change" << endl;

       all_species[species_to_do].happening_population_popchange_this_species(sd_normal_distribution_pop_change, map1, alleles_adaptation_coef);
       all_species[species_to_do] = all_species[species_to_do]; // this line updates the all_species vector
       total_popchange_events = total_popchange_events + 1;
     }

     if (event_to_do == "contraction")
     {
       //  cout << "                  I will contract range " << endl;
       all_species[species_to_do].happening_contraction(t, map1, extirpation_depen_temperature);
       total_contraction_events = total_contraction_events + 1;
       //all_species[species_to_do] = all_species[species_to_do]; // this line updates the all_species vector
     }
     vector_events_tookplace.push_back(event_to_do);

     // cout << "super summary from here" << endl;
     // for (int iij = 0; iij < all_species.size(); ++iij){
     //   cout << "id: "<< all_species[iij].id << " range " <<  all_species[iij].range <<" x_coordinate_last_event  " <<  all_species[iij].x_coordinate_last_event << "y_coordinate_last_event " << all_species[iij].y_coordinate_last_event << endl;
     //   cout << "new_species.temperature_optimum.size() " << all_species[iij].temperature_optimum.size() << endl;
     //   cout << "new_species.southernmost "<< all_species[iij].southernmost << "new_species.northernmost " << all_species[iij].northernmost << endl;
     //   cout << "total_rate " <<all_species[iij].total_rate << endl;
     //   cout << "new_species.birth "<< all_species[iij].birth << "new_species.trait_state " << all_species[iij].trait_state << endl;
     //   cout << "new_species.total_pop_size "<< all_species[iij].total_pop_size << endl;
     //   for (int i = 0; i < all_species[iij].presence.size(); ++i)
     //   {
     //     cout << "presence here in  " << all_species[iij].presence[i].x << " " << all_species[iij].presence[i].y  << endl;
     //     cout <<"populations_this_species: " << all_species[iij].populations_this_species[i].pop_size << endl;
     //     cout <<"temperature_optimum: " << all_species[iij].temperature_optimum[i] << endl;
     //     cout <<"computed_rate_based_on_temperature: " << all_species[iij].computed_rate_based_on_temperature[i] << endl;
     //
     //   }
     // }


     // cout << "northermost: "<< all_species[species_to_do].northernmost << "south: " << all_species[species_to_do].southernmost << endl;
     if (total_num_populations == 1 && event_to_do == "contraction")
     {
       cout << "total annihilation of the clade" << endl;
       break;
     }


     bool problem_zero_popsize;
     problem_zero_popsize = false;
     for(int iji = 0; iji < all_species.size(); ++iji)
     {
       species do_this_species;
       do_this_species = all_species[iji];
       if(do_this_species.alive)
       {

         for(int iij = 0; iij < do_this_species.presence.size(); ++iij)
         {
           population_structure do_this_population;
           do_this_population = do_this_species.populations_this_species[iij];
           // cout << "born at: " << do_this_species.birthplace.x << "," << do_this_species.birthplace.y << endl;
           if(do_this_population.pop_size <= 0)
           {
             problem_zero_popsize = true;
             cout << "this is supposed to be zero: " << do_this_population.pop_size << endl;
           }
         }
       }

     }
     //show_all_species_data(all_species);
     if(problem_zero_popsize)
     {
       cout << "+++++++++++++++++++ some population is below zero right after event of  " << event_to_do << endl;
       stop("some issue with population below zero");
       break;
     } // end of checks
   } // End of While loop

   int final_richness;
   final_richness = 0;
   int final_numb_pop;
   final_numb_pop = 0;
   int final_indviduals;
   final_indviduals = 0;
   for(int iji = 0; iji < all_species.size(); ++iji)
   {
     if(all_species[iji].alive)
     {
       total_geneflow_events = total_geneflow_events + all_species[iji].succesful_geneflow_events;
       final_richness = final_richness + 1;
       final_numb_pop = final_numb_pop + all_species[iji].range;
       final_indviduals =  final_indviduals + all_species[iji].total_pop_size;
     }
   }


   cout << "____ Summary____" << endl;
   if(pending_change_in_rates == false)
   {
     cout << "change in rates and/or temperature did take place" << endl;
     if(time_percent_stop_after_first_equilibrium_and_disturbance  != 0){
       cout << "first_equi_at: " << first_equi_at << endl;
       cout << "and it let the model run some more time and stop at: " << simulated_time << " which is "<< time_percent_stop_after_first_equilibrium_and_disturbance <<" % more of the time the simulation had already run for" << endl;
       cout << " richness_at_equilibrium: " << richness_at_equilibrium << " populations_at_equilibrium: " << populations_at_equilibrium << " individuals_at_equilibrium: " << individuals_at_equilibrium <<endl;

     }

   }
   cout << "events took place: " << endl;
   cout << "total_expansion_events " << total_expansion_events << endl;
   cout << "total_contraction_events " << total_contraction_events << endl;
   cout << "total_mutation_events " << total_mutation_events << endl;
   cout << "total_traitevolution_events" << total_traitevolution_events << endl;
   cout << "attempted_geneflow_events " << attempted_geneflow_events << endl;
   cout << "total_geneflow_events " << total_geneflow_events << endl;
   cout << "total_speciation_events " << total_speciation_events << endl;
   cout << "total_popchange_events " << total_popchange_events << endl;


   cout << "time: " << t << " cycle: " << cycles << " richness:" << final_richness <<  " populations: " << final_numb_pop << " indviduals: " << final_indviduals<< " ind_saturation %: " << full_saturation_indi << endl;
   to_show_richness_map(y_max,x_max,show_richness_map,all_species,map1);

   bool no_failure;
   no_failure = final_check(y_max,x_max,all_species,map1);
   if(no_failure)
   {
     cout << "OK" << endl;
     cout << "end of simulation, now we are in Github too" << endl;
   }

   // for RcPP

   List model_output = List::create();
   model_output = get_me_output(y_max,x_max,all_species,t);

   // end for Rcpp
   //
   //       int model_output;
   //    model_output = 3;
   return model_output;
 }
