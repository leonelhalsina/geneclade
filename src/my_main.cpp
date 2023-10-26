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
#include "stocc.h"
using namespace std;
using namespace Rcpp;

//' @export
 // [[Rcpp::export]]
 List try_random_wallenius (int number_alelles, IntegerVector mlist_vector, NumericVector wlist_vector,int items_to_take){

   int32_t xlist[number_alelles];     // output vector




   //int32_t mlist[] = {30,60,50,20};  // number of balls of each color

   int32_t mlist[number_alelles];  // number of balls of each color

   for(int ii = 0; ii < mlist_vector.size(); ++ii){
     mlist[ii] = mlist_vector[ii];
   }
   //
   // mlist[0] = this_is_vector[0];
   // mlist[1] = this_is_vector[1];
   // mlist[2] = this_is_vector[2];
   // mlist[3] = this_is_vector[3];
   //cout<<mlist[3]<<endl;
//
//    vector <double> wlist_vector;
//    wlist_vector = {1,5,6,2};

   double wlist[number_alelles];
     wlist [0] = wlist_vector[0];   // weight of each color
     wlist [1] = wlist_vector[1];   // weight of each color
     wlist [2] = wlist_vector[2];   // weight of each color
     wlist [3] = wlist_vector[3];   // weight of each color

   int32_t n = items_to_take;                       // number of balls to pick

   StochasticLib3 randomwallenius((int32_t)time(0));

   randomwallenius.MultiWalleniusNCHyp(xlist, mlist, wlist, n, number_alelles);
cout << "wlist"<<wlist[3] << endl;
cout << "wlist"<<mlist[3] << endl;
cout << "wlist_vector"<<wlist_vector[3] << endl;
   // update sums
   // int i;
   // double my_probab;
   //
   // for (i=0; i<number_alelles; i++) {
   //   my_probab = my_probab + xlist[i];
   // }
   //

   List list_all_species_distribution = List::create();
   list_all_species_distribution.push_back(xlist[0]);
   list_all_species_distribution.push_back(xlist[1]);
   list_all_species_distribution.push_back(xlist[2]);
   list_all_species_distribution.push_back(xlist[3]);
   List model_output = List::create(Named("Distribution") = list_all_species_distribution);

   // end for Rcpp
   return model_output;

 }
