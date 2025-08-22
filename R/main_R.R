#' Prepares an object to be used as starting point for a simulation. It requieres the user to provide spatial and geographic information of each species.
#' @title Initialisation for a simulation when defaults are not adequate.
#' @param all_x Map coordinates (vector) in x for the first populations in the simulation.
#' @param all_y Map coordinates (vector) in y for the first populations in the simulation.
#' @param number_spp number of species to start the simulation with.
#' @param all_traits initial trait state, vector, same length than number_spp
#' @param all_births birthdate of each species, vector, same length than number_spp
#' @param all_IDs ids of each species in the system, vector, same length than number_spp
#' @param all_parents vector of the ids for the parents to the species, same length than number_spp
#' @param all_deaths vector of record if deaths of species, 0 means alive, same length than number_spp
#' @param all_ranges vector of number of cells where each species is present. Currently, the model only handle species starting with range size of one.
#' @param all_alleles vector of allelic frequencies for allele a1, b1, c1...to y1.
#' @param all_alleles_neutral vector of allelic frequencies for allele a2, b2, c2...to y2.
#' @param all_popsize vector of population size for each species, it should be consistent with all_alleles and all_alleles_neutral.
#' @param starting_time starting simulated time, it should be consistent with the times in all_births and all_deaths arguments.
#' @return List of variables needed to start a simulation.
#' @examples
#'
#' number_spp <- 190
#' all_x <- c(rep(2:20,10))#3#3#10 #  initial population column
#'all_y <- c(rep(3,19), rep(5,19),rep(7,19),rep(9,19),rep(11,19),
#'           rep(13,19),rep(15,19),rep(17,19),rep(19,19),rep(21,19))#3 #3#10 # initial population row
#'all_traits <- c(rep(100,190)) # initial trait value
#'all_births <- c(rep(0,190))
#'all_IDs <- c(1:190)
#'all_parents <- c(rep(0,190))
#'all_deaths <- c(rep(0,190))
#'all_ranges <- c(rep(1,190))
#'all_alleles <- rep(6,15) # allelic frequencies for allele a, b, c, d, and e.
#'all_alleles_neutral <- rep(6,15)
#'all_popsize <- c(rep(90,190))#population size for the first population. It should be the same than sum(all_alleles)
#' starting_time <- 20
#'make_advanced_initialization (number_spp,
#'                                         all_x,
#'                                         all_y,
#'                                         all_traits,
#'                                         all_births,
#'                                         all_IDs,
#'                                         all_parents,
#'                                         all_deaths,
#'                                         all_ranges,
#'                                         all_alleles,
#'                                         all_alleles_neutral,
#'                                         all_popsize,
#'                                         starting_time)
#'
#'
#' @export
make_advanced_initialization <- function (number_spp,
                                          all_x,
                                          all_y,
                                          all_traits,
                                          all_births,
                                          all_IDs,
                                          all_parents,
                                          all_deaths,
                                          all_ranges,
                                          all_alleles,
                                          all_alleles_neutral,
                                          all_popsize,
                                          starting_time){


  if(any(all_ranges > 1)){
    stop("Currently, the model only handle species starting with range size of one.")
  }

  advanced_initialization <- list()
  advanced_initialization[[1]] <- number_spp
  advanced_initialization[[2]] <- all_x
  advanced_initialization[[3]] <- all_y
  advanced_initialization[[4]] <- all_traits
  advanced_initialization[[5]] <- all_births
  advanced_initialization[[6]] <- all_IDs
  advanced_initialization[[7]] <- all_parents
  advanced_initialization[[8]] <- all_deaths
  advanced_initialization[[9]] <- all_ranges
  advanced_initialization[[10]] <- all_alleles
  advanced_initialization[[11]] <- all_alleles_neutral
  advanced_initialization[[12]] <- all_popsize
  advanced_initialization[[13]] <- starting_time
  return(advanced_initialization)
}


#' Runs a spatially explicit simulation where alleles, populations and species are modelled. Demographic, dispersal and evolutionary processes are simulated in continuous time.
#' @title Run geneclade model using rates for processes taking place at ecological and evolutionary time scales in a spatial context.
#' @param position_start_x Map coordinate in x for the first population in the simulation.
#' @param position_start_y Map coordinate in y for the first population in the simulation.
#' @param advanced_initialization When initilization needs to be different from default one (one sp with one population of size 150, all alles in the same frequency), user needs to input a list created with make_advanced_initialization() function. Default for this argument is NULL
#' @param max_spp Maximum expected number of species.
#' @param simulated_time The time the simulation will run for. The scale of the simulated time is highly dependent on selected rates.
#' @param condition_to_stop The simulation will stop when the number of species is met ("richness") or when the simulated time is up ("time"). Default is "richness".
#' @param stop_time_after_change Indicates the % of time that the simulation will be let run further, considering the time it needed to reach the first equilibirum. Default is 0 which means that the simulation will run until reaching a second equilibrium after map change (see below).
#' @param time_slices A vector with the simulated time points at which the current state of the simulation will be saved in memory to retrive at the end of the simulation. Notice that you might need to know more or less the temporal scale of the total simulation.
#' @param maximum_cycles Useful parameter to allow a speficic number of events (every cycle has an event e.g., colonization, geneflow) before stopping the simulation.
#' @param rate_speciation Per-population rate of speciation rate (lambda)
#' @param rate_colonisation Per-population rate of colonisation rate (gamma). Dispersal of individuals of one cell to an adjacent which is below carrying capacity K.
#' @param rate_extirpation Per-population rate local extinction rate (mu). For a given species, a population is extirpated from a cell.
#' @param rate_geneflow Per-population rate of geneflow (sigma). One population exchange individuals with an adjacent one.
#' @param rate_demographicchange Vector of length 2. First element is the per-population rate of change in population size (delta). Population can grow or shrink. Second element is the standard deviation for the normal distribution from which the next population size will be drawn. Example rate_demographicchange <- c(0.01,60)
#' @param rate_traitevolution Vector of length 2. First element is the per-population rate of change in trait value. Second element is the standard deviation for the normal distribution from which the next trait value  will be drawn. Example rate_traitevolution <- c(0.01,1). Notice that all individuals of a species have the same trait state.
#' @param rate_mutation Per-population rate of mutation (q) where one individual changes its two loci to different ones.
#' @param alleles_optimum_enviroment Vector with 25 integer values representing the environmental optimum for each allele from a1 to y1. Notice that if all 25 are the same, the model will feature no differences in local adaptation and will therefore be neutral. These values should make sense with the cell values in the map_environment_1, see below.
#' @param percentage_geneflow This arguments regulates how large, in each geneflow event, is the exchange of individuals between two populations. Because it is a percentage, the exchange can be assymetrical. Default is 10.
#' @param vicariant_speciation Boolean. Whether or not species speciate by vicariance. False will result in point-mutation speciation.
#' @param manual_speciation_events_timing In case it is needed, for some reason, to force the simulation to have speciation events at given times, this vector should provide the time for those events. This means that the length of this vector will indicate the total number of speciation events. This will automatically turn rate_speciation = 0. This option is disable by default (i.e., manual_speciation_events_timing)
#' @param growth_only Boolean. Should population fluctiations only lead to population growth? Useful to ensure local saturation of species.
#' @param unlink_range_to If one wishes to make speciation, colonization, or both independent from the number of populations use unlink_range_to <- c("colonization","speciation"). We suggest not using this option. Default is NULL.
#' @param x_max Length of x axis in the map.
#' @param y_max Length of y axis in the map.
#' @param map_k_1 dataframe to be the first map to run simulations on (dimensions x_max and y_max). Inhabitable cells marked with -9. Cell values represent local carrying capacity K.
#' @param map_k_2 dataframe to be the second map to run simulations on (dimensions x_max and y_max). Inhabitable cells marked with -9. Cell values represent local carrying capacity K. This map will replace map_k_2 at THIS TIME!!!!. If there is no interest in changing maps, please do: map_k_2 <- map_k_1
#' @param map_environment_1 dataframe to be the first map to run simulations on (dimensions x_max and y_max). Inhabitable cells marked with -9. Cell values represent local environmental conditions e.g., temperature
#' @param map_environment_2 dataframe to be the second map to run simulations on (dimensions x_max and y_max). Inhabitable cells marked with -9. Cell values represent local environmental conditions e.g., temperature. This map will replace map_environment_1 at THIS TIME!!!!.If there is no interest in changing maps, please do: map_environment_2 <- map_environment_1


#' @return List of species with detailed information on the variables tracked over time, ready for processing with THESE FUNCTIONS.
#' @examples
#'# Example of how to set the arguments for a Maximum Likelihood search.
#'library(geneclade)
#' You can use this output in the plotting function: plot_biogeo_reconst()
#' @export



run_simulation <- function(position_start_x,
                           position_start_y,
                           advanced_initialization = NULL,
                           max_spp,
                           simulated_time,
                           condition_to_stop = "richness",
                           stop_time_after_change = 0,
                           time_slices,
                           maximum_cycles = 2e+09,
                           rate_speciation,
                           rate_colonisation,
                           rate_extirpation,
                           rate_geneflow,
                           rate_demographicchange,
                           rate_traitevolution,
                           rate_mutation,
                           alleles_optimum_enviroment,
                           percentage_geneflow = 10,
                           vicariant_speciation,
                           manual_speciation_events_timing = 0,
                           growth_only = TRUE,
                           unlink_range_to = NULL,
                           x_max,
                           y_max,
                           map_k_1,
                           map_k_2,
                           map_environment_1,
                           map_environment_2){
  # initial

  if(is.null(advanced_initialization)){
    starting_time <- 0
    number_spp <- 1
    all_x <- position_start_x
    all_y <- position_start_y
    all_traits <- 100 # initial trait value
    all_births <- 0
    all_IDs <- 1
    all_parents <- 0
    all_deaths <- 0
    all_ranges <- 1
    all_alleles <- rep(6,25) # allelic frequencies for allele a1, b1, c1...to y1.
    all_alleles_neutral <- rep(6,25) # allelic frequencies for allele a2, b2, c2...to y2.
    all_popsize <- 150 #population size for the first population. It should be the same than sum(all_alleles)


  } else {
    if(class(advanced_initialization) != "list"){
      stop("advanced_initialization needs to be list, please use the function make_advanced_initialization()")
    }
    number_spp <- advanced_initialization[[1]]
    all_x <- advanced_initialization[[2]]
    all_y <- advanced_initialization[[3]]
    all_traits <- advanced_initialization[[4]]
    all_births <- advanced_initialization[[5]]
    all_IDs <- advanced_initialization[[6]]
    all_parents <- advanced_initialization[[7]]
    all_deaths <- advanced_initialization[[8]]
    all_ranges <- advanced_initialization[[9]]
    all_alleles <- advanced_initialization[[10]]
    all_alleles_neutral <- advanced_initialization[[11]]
    all_popsize <- advanced_initialization[[12]]
    starting_time <- advanced_initialization[[13]]
  }

  # simulation control
  max_spp <- max_spp
  simulated_time <- simulated_time
  time_slices <- time_slices
  condition_to_stop <- condition_to_stop
  maximum_cycles <- maximum_cycles
  time_percent_stop_after_first_equilibrium_and_disturbance <- stop_time_after_change



  # rates
  lambda <- rate_speciation
  the_gammas <- rate_colonisation
  the_mus <- rate_extirpation
  geneflow_rate <- rate_geneflow
  popchange_rate <- rate_demographicchange[1]
  q <- rate_traitevolution[1]
  mutation_rate <- rate_mutation
  if(manual_speciation_events_timing != 0){
    lambda <- 0
  }

  # adjustment of processes
  percentage_flow <- percentage_geneflow
  vicariant_speciation <- vicariant_speciation
  alleles_adaptation_coef2 <- alleles_optimum_enviroment

  # map features
  x_max <- x_max
  y_max <- y_max

  text_k_map <- map_k_1
  text_k_map2 <- map_k_2
  map_temperature <- map_environment_1
  map_temperature2 <- map_environment_2


  mean_normal_distribution_traitevol <- 0 # currently not used as the trait state value is used as the mean for the normal distribution to sample from.
  sd_normal_distribution_traitevol  <- rate_traitevolution[2] #I used 2 , trait evolution parameters
  sd_normal_distribution_pop_change <- rate_demographicchange[2] # I used 57


  map_k_vector <- transform_map_vector (text_k_map)
  map_k_vector2 <- transform_map_vector (text_k_map2)
  map_temperature_vector <- transform_map_vector(map_temperature)
  map_temperature_vector2 <- transform_map_vector(map_temperature2)


  restiction_par <- 2  # this parameter might not be important
  do_change_map_rates <- TRUE
  the_seed <- sample(1:9000,1)
  use_k <- TRUE
  show_richness_map <- "none"

  colonization_depen_temperature <- TRUE
  extirpation_depen <- "random" #"temperature" "random" "popsize"
  species_trait_state <- "both"#"gamma" "geneflow" "both"

  v <- 5 # likely to be removed


  speciation_rangesize_unlinked <- FALSE
  colonization_rangesize_unlinked <- FALSE

  if(!is.null(unlink_range_to)){
    for(i in 1:length(unlink_range_to)){
      if(unlink_range_to[i] == "colonization"){
        colonization_rangesize_unlinked <- TRUE
      }
      if(unlink_range_to[i] == "speciation"){
        speciation_rangesize_unlinked <- TRUE
      }
    }
  }





  do_simulation(map_k_vector,
                map_k_vector2,
                map_temperature_vector,
                map_temperature_vector2,
                extirpation_depen,
                colonization_depen_temperature,
                x_max,
                y_max,
                all_x,
                all_y,
                all_IDs,
                all_parents,
                all_births,
                all_deaths,
                all_traits,
                all_ranges,
                all_alleles,
                all_alleles_neutral,
                all_popsize,
                number_spp,
                the_seed,
                mutation_rate,
                percentage_flow,
                geneflow_rate,
                popchange_rate,
                the_gammas,
                the_mus,
                q,
                lambda,
                species_trait_state,
                sd_normal_distribution_traitevol,
                mean_normal_distribution_traitevol,
                sd_normal_distribution_pop_change,
                growth_only,
                starting_time,
                simulated_time,
                max_spp,
                maximum_cycles,
                use_k,
                restiction_par,
                show_richness_map,
                v,
                alleles_adaptation_coef2,
                do_change_map_rates,
                vicariant_speciation,
                speciation_rangesize_unlinked,
                colonization_rangesize_unlinked,
                time_percent_stop_after_first_equilibrium_and_disturbance,
                condition_to_stop,
                time_slices,
                manual_speciation_events_timing)

}

#
# all_x <- 17#3#3#10 #  initial population column
# all_y <- 17#3 #3#10 # initial population row
#
# position_start_x <- 17
#
# position_start_y <- 17
# max_spp <- 20
# simulated_time <- 250
#
# time_slices <- c(10,20)
#
# rate_speciation <- 0.00005
# rate_colonisation <- 5
# rate_extirpation <- 0.0000001
# rate_geneflow <- 0.01
# rate_demographicchange <- 1
# rate_traitevolution <- 0
# rate_mutation <- rate_geneflow/10
# alleles_optimum_enviroment <-  rep(9,25)
# vicariant_speciation <- TRUE
# x_max <- 42
# y_max <- 42
# map_k_1 <- read.table(paste0("k_map_uniform.txt"))
# map_k_2 <- map_k_1
# map_environment_1 <- read.table(paste0("temperature_map_uniform.txt"))
# map_environment_2 <- map_environment_1
#
#
# maximum_cycles <- 200000
# simulation_raw <- run_simulation (position_start_x,
#                          position_start_y,
#                          advanced_initialization = NULL,
#                          max_spp,
#                          simulated_time,
#                          condition_to_stop = "richness",
#                          stop_time_after_change = 0,
#                          time_slices,
#                          maximum_cycles = maximum_cycles,
#                          rate_speciation,
#                          rate_colonisation,
#                          rate_extirpation,
#                          rate_geneflow,
#                          rate_demographicchange,
#                          rate_traitevolution,
#                          rate_mutation,
#                          alleles_optimum_enviroment,
#                          percentage_geneflow = 10,
#                          vicariant_speciation,
#                          manual_speciation_events_timing = 0,
#                          growth_only = TRUE,
#                          unlink_range_to = NULL,
#                          x_max,
#                          y_max,
#                          map_k_1,
#                          map_k_2,
#                          map_environment_1,
#                          map_environment_2)
#
#
# for(ikk in 1:length(simulation_raw)){
#
# }
# simulation_raw_this_timeslice <- simulation_raw[[ikk]]
# simulation_raw_this_timeslice <- make_one_list_fromRaw(simulation_raw_this_timeslice)
#
# time_simulated_from_output <- simulation_raw_this_timeslice$total_time[[1]]
# list_species_from_cpp <- make_list_species_fromcpp_toR(simulation_raw_this_timeslice)
