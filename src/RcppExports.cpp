// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dotRcpp_I
NumericVector dotRcpp_I(NumericVector x1, NumericVector x2);
RcppExport SEXP _geneclade_dotRcpp_I(SEXP x1SEXP, SEXP x2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x2(x2SEXP);
    rcpp_result_gen = Rcpp::wrap(dotRcpp_I(x1, x2));
    return rcpp_result_gen;
END_RCPP
}
// do_simulation
List do_simulation(IntegerVector map_elevation_vector, IntegerVector map_k_vector, IntegerVector map_temperature_vector, IntegerVector map_temperature_vector2, bool extirpation_depen_temperature, bool colonization_depen_temperature, int x_max, int y_max, IntegerVector all_x, IntegerVector all_y, IntegerVector all_IDs, IntegerVector all_parents, NumericVector all_births, NumericVector all_deaths, NumericVector all_traits, IntegerVector all_ranges, IntegerVector all_alleles, IntegerVector all_alleles_neutral, IntegerVector all_popsize, int number_spp, int the_seed, double mutation_rate, double percentage_flow, double geneflow_rate, double popchange_rate, NumericVector the_gammas, NumericVector the_mus, double q, double lambda, bool species_trait_state_gamma, double sd_normal_distribution_traitevol, double mean_normal_distribution_traitevol, double sd_normal_distribution_pop_change, double starting_time, double simulated_time, int max_spp, int maximum_cycles, bool use_k, double restiction_par, std::string show_richness_map, double v, IntegerVector alleles_adaptation_coef2, bool do_change_map_rates, bool vicariant_speciation);
RcppExport SEXP _geneclade_do_simulation(SEXP map_elevation_vectorSEXP, SEXP map_k_vectorSEXP, SEXP map_temperature_vectorSEXP, SEXP map_temperature_vector2SEXP, SEXP extirpation_depen_temperatureSEXP, SEXP colonization_depen_temperatureSEXP, SEXP x_maxSEXP, SEXP y_maxSEXP, SEXP all_xSEXP, SEXP all_ySEXP, SEXP all_IDsSEXP, SEXP all_parentsSEXP, SEXP all_birthsSEXP, SEXP all_deathsSEXP, SEXP all_traitsSEXP, SEXP all_rangesSEXP, SEXP all_allelesSEXP, SEXP all_alleles_neutralSEXP, SEXP all_popsizeSEXP, SEXP number_sppSEXP, SEXP the_seedSEXP, SEXP mutation_rateSEXP, SEXP percentage_flowSEXP, SEXP geneflow_rateSEXP, SEXP popchange_rateSEXP, SEXP the_gammasSEXP, SEXP the_musSEXP, SEXP qSEXP, SEXP lambdaSEXP, SEXP species_trait_state_gammaSEXP, SEXP sd_normal_distribution_traitevolSEXP, SEXP mean_normal_distribution_traitevolSEXP, SEXP sd_normal_distribution_pop_changeSEXP, SEXP starting_timeSEXP, SEXP simulated_timeSEXP, SEXP max_sppSEXP, SEXP maximum_cyclesSEXP, SEXP use_kSEXP, SEXP restiction_parSEXP, SEXP show_richness_mapSEXP, SEXP vSEXP, SEXP alleles_adaptation_coef2SEXP, SEXP do_change_map_ratesSEXP, SEXP vicariant_speciationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type map_elevation_vector(map_elevation_vectorSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type map_k_vector(map_k_vectorSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type map_temperature_vector(map_temperature_vectorSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type map_temperature_vector2(map_temperature_vector2SEXP);
    Rcpp::traits::input_parameter< bool >::type extirpation_depen_temperature(extirpation_depen_temperatureSEXP);
    Rcpp::traits::input_parameter< bool >::type colonization_depen_temperature(colonization_depen_temperatureSEXP);
    Rcpp::traits::input_parameter< int >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< int >::type y_max(y_maxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type all_x(all_xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type all_y(all_ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type all_IDs(all_IDsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type all_parents(all_parentsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type all_births(all_birthsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type all_deaths(all_deathsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type all_traits(all_traitsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type all_ranges(all_rangesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type all_alleles(all_allelesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type all_alleles_neutral(all_alleles_neutralSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type all_popsize(all_popsizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_spp(number_sppSEXP);
    Rcpp::traits::input_parameter< int >::type the_seed(the_seedSEXP);
    Rcpp::traits::input_parameter< double >::type mutation_rate(mutation_rateSEXP);
    Rcpp::traits::input_parameter< double >::type percentage_flow(percentage_flowSEXP);
    Rcpp::traits::input_parameter< double >::type geneflow_rate(geneflow_rateSEXP);
    Rcpp::traits::input_parameter< double >::type popchange_rate(popchange_rateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type the_gammas(the_gammasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type the_mus(the_musSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type species_trait_state_gamma(species_trait_state_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sd_normal_distribution_traitevol(sd_normal_distribution_traitevolSEXP);
    Rcpp::traits::input_parameter< double >::type mean_normal_distribution_traitevol(mean_normal_distribution_traitevolSEXP);
    Rcpp::traits::input_parameter< double >::type sd_normal_distribution_pop_change(sd_normal_distribution_pop_changeSEXP);
    Rcpp::traits::input_parameter< double >::type starting_time(starting_timeSEXP);
    Rcpp::traits::input_parameter< double >::type simulated_time(simulated_timeSEXP);
    Rcpp::traits::input_parameter< int >::type max_spp(max_sppSEXP);
    Rcpp::traits::input_parameter< int >::type maximum_cycles(maximum_cyclesSEXP);
    Rcpp::traits::input_parameter< bool >::type use_k(use_kSEXP);
    Rcpp::traits::input_parameter< double >::type restiction_par(restiction_parSEXP);
    Rcpp::traits::input_parameter< std::string >::type show_richness_map(show_richness_mapSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles_adaptation_coef2(alleles_adaptation_coef2SEXP);
    Rcpp::traits::input_parameter< bool >::type do_change_map_rates(do_change_map_ratesSEXP);
    Rcpp::traits::input_parameter< bool >::type vicariant_speciation(vicariant_speciationSEXP);
    rcpp_result_gen = Rcpp::wrap(do_simulation(map_elevation_vector, map_k_vector, map_temperature_vector, map_temperature_vector2, extirpation_depen_temperature, colonization_depen_temperature, x_max, y_max, all_x, all_y, all_IDs, all_parents, all_births, all_deaths, all_traits, all_ranges, all_alleles, all_alleles_neutral, all_popsize, number_spp, the_seed, mutation_rate, percentage_flow, geneflow_rate, popchange_rate, the_gammas, the_mus, q, lambda, species_trait_state_gamma, sd_normal_distribution_traitevol, mean_normal_distribution_traitevol, sd_normal_distribution_pop_change, starting_time, simulated_time, max_spp, maximum_cycles, use_k, restiction_par, show_richness_map, v, alleles_adaptation_coef2, do_change_map_rates, vicariant_speciation));
    return rcpp_result_gen;
END_RCPP
}
// try_random_wallenius
List try_random_wallenius(int number_alelles, IntegerVector mlist_vector, NumericVector wlist_vector, int items_to_take);
RcppExport SEXP _geneclade_try_random_wallenius(SEXP number_alellesSEXP, SEXP mlist_vectorSEXP, SEXP wlist_vectorSEXP, SEXP items_to_takeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type number_alelles(number_alellesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type mlist_vector(mlist_vectorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wlist_vector(wlist_vectorSEXP);
    Rcpp::traits::input_parameter< int >::type items_to_take(items_to_takeSEXP);
    rcpp_result_gen = Rcpp::wrap(try_random_wallenius(number_alelles, mlist_vector, wlist_vector, items_to_take));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_geneclade_dotRcpp_I", (DL_FUNC) &_geneclade_dotRcpp_I, 2},
    {"_geneclade_do_simulation", (DL_FUNC) &_geneclade_do_simulation, 44},
    {"_geneclade_try_random_wallenius", (DL_FUNC) &_geneclade_try_random_wallenius, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_geneclade(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
