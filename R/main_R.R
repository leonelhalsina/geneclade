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
#' @param max_spp Total expected number of species.
#' @param simulated_time The time the simulation will run for. The scale of the simulated time is highly dependent on selected rates.
#' @param condition_to_stop The simulation will stop when the number of species is met ("richness") or when the simulated time is up ("time"). Default is "richness".
#' @param time_slices A vector with the simulated time points at which the current state of the simulation will be saved in memory to retrive at the end of the simulation. Notice that you might need to know more or less the temporal scale of the total simulation.
#' @param maximum_cycles Useful parameter to allow a speficic number of events (every cycle has an event e.g., colonization, geneflow) before stopping the simulation.
#' @param rate_speciation Per-population rate of speciation rate (lambda)
#' @param rate_colonisation Per-population rate of colonisation rate (gamma). Dispersal of individuals of one cell to an adjacent which is below carrying capacity K.
#' @param rate_extirpation Per-population rate local extinction rate (mu). For a given species, a population is extirpated from a cell.
#' @param rate_geneflow Per-population rate of geneflow (sigma). One population exchange individuals with an adjacent one.
#' @param rate_demographicchange Per-population rate of change in population size (delta). Population can grow or shrink.
#' @param rate_traitevolution Per-population rate of change in trait value. Their is one value per species.
#' @param rate_mutation Per-population rate of mutation (q) where one individual changes its two loci to different ones.

#' @param percentage_geneflow This arguments regulates how large, in each geneflow event, is the exchange of individuals between two populations. Because it is a percentage, the exchange can be assymetrical. Default is 10.
#' @param vicariant_speciation Boolean. Whether or not species speciate by vicariance. False will result in point-mutation speciation.
#' @param x_max Length of x axis in the map
#' @param y_max Length of y axis in the map
#' @param map_k_1 dataframe to be the first map to run simulations on (dimensions x_max and y_max). Inhabitable cells marked with -9. Cell values represent local carrying capacity K.
#' @param map_k_2 dataframe to be the second map to run simulations on (dimensions x_max and y_max). Inhabitable cells marked with -9. Cell values represent local carrying capacity K. This map will replace map_k_2 at THIS TIME!!!!

#' @param map_environment_1 dataframe to be the first map to run simulations on (dimensions x_max and y_max). Inhabitable cells marked with -9. Cell values represent local environmental conditions e.g., temperature
#' @param map_environment_2 dataframe to be the second map to run simulations on (dimensions x_max and y_max). Inhabitable cells marked with -9. Cell values represent local environmental conditions e.g., temperature. This map will replace map_environment_1 at THIS TIME!!!!


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
                           time_slices,
                           maximum_cycles = 2e+09,
                           rate_speciation,
                           rate_colonisation,
                           rate_extirpation,
                           rate_geneflow,
                           rate_demographicchange,
                           rate_traitevolution,
                           rate_mutation,
                           percentage_geneflow = 10,
                           vicariant_speciation,
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

  # rates
  lambda <- rate_speciation
  the_gammas <- rate_colonisation
  the_mus <- rate_extirpation
  geneflow_rate <- rate_geneflow
  popchange_rate <- rate_demographicchange
  q <- rate_traitevolution
  mutation_rate <- rate_mutation

  # adjustment of processes
  percentage_flow <- percentage_geneflow
  vicariant_speciation <- vicariant_speciation

  # map features
  x_max <- x_max
  y_max <- y_max

  text_k_map <- map_k_1
  text_k_map2 <- map_k_2
  map_temperature <- map_environment_1
  map_temperature2 <- map_environment_2

  ##? potential parameters to go to advanced settings.
  mean_normal_distribution_traitevol <- 0 # currently not used as the trait state value is used as the mean for the normal distribution to sample from.
  sd_normal_distribution_traitevol  <- 2 # trait evolution parameters
  sd_normal_distribution_pop_change <- 57


  map_k_vector <- transform_map_vector (text_k_map)
  map_k_vector2 <- transform_map_vector (text_k_map2)
  map_temperature_vector <- transform_map_vector(map_temperature)
  map_temperature_vector2 <- transform_map_vector(map_temperature2)

  do_simulation(map_k_vector,map_k_vector2,  map_temperature_vector, map_temperature_vector2,
                extirpation_depen,  colonization_depen_temperature, x_max,  y_max,  all_x,  all_y,
                all_IDs,  all_parents,  all_births,  all_deaths,
                all_traits,  all_ranges,  all_alleles, all_alleles_neutral,all_popsize,
                number_spp,  the_seed,
                mutation_rate, percentage_flow, geneflow_rate, popchange_rate, the_gammas,  the_mus,
                q,  lambda, species_trait_state,sd_normal_distribution_traitevol,mean_normal_distribution_traitevol,
                sd_normal_distribution_pop_change,growth_only,  starting_time,  simulated_time,max_spp,
                maximum_cycles,  use_k,  restiction_par, show_richness_map,  v,
                alleles_adaptation_coef2,
                do_change_map_rates,vicariant_speciation,speciation_rangesize_unlinked,
                colonization_rangesize_unlinked,time_percent_stop_after_first_equilibrium_and_disturbance, condition_to_stop,
                time_slices,manual_speciation_events_timing)

}




#i <- 3

#setwd("C:/Users/s06lh9/Research Fellowship/The Science/Lesley group/agentBasedModel/simulations")

do_things <- function(i){
  if(.Platform$OS.type == "unix"){  # when in the cluster
    .libPaths("/uoa/scratch/users/s06lh9/R/x86_64-pc-linux-gnu-library/4.1")
  }


  library("Rcpp")
  source("functions_output.R")
  source("convert.R")
  library(geneclade)
  print(packageVersion("geneclade"))

  max_spp <- 20
  this_lambda <- "low"

  options_what_k <- "uniform" #c("change","unchanged") # "uniform"
  options_what_temp <- "uniform" #c("change","unchanged") "uniform"
  options_vicariance_case <- c("On","Off")
  options_geneflow <- c("low","high")
  options_adaptation <- c(FALSE)
  options_replicates <- 1:30
  all_combinations <- tidyr::expand_grid(options_what_k, options_what_temp, options_vicariance_case,options_adaptation,options_geneflow,options_replicates)

  what_k <- all_combinations$options_what_k [i]
  what_temp <-  all_combinations$options_what_temp [i]
  this_vicariance <- all_combinations$options_vicariance_case[i]
  what_adaptation <- all_combinations$options_adaptation [i]
  this_gene_flow <- all_combinations$options_geneflow [i]
  this_replicate <- all_combinations$options_replicates[i]
  adaptation_on <- what_adaptation




  growth_only <- TRUE
  condition_to_stop <- "richness"# "richness" "time"

  species_trait_state <- "both"#"gamma" "geneflow" "both"


  simulated_time <- 250 #750 #750 #1208#750 #1201 #1250#1420    1250

  starting_one_spp <- TRUE
  #num_sim <- 10

  time_slices <- seq(from=10, to=2000,100)




  the_mus <- c(0.0001,0.0001)#0.1   #0.01. Two values: before and after the change
  the_gammas <- c(5,5) # two values: before and after the change c(0.01,0.01)




  x_max <- 42
  y_max <- 42




  speciation_rangesize_unlinked <- FALSE
  colonization_rangesize_unlinked <- FALSE
  colonization_depen_temperature <- TRUE
  extirpation_depen <- "random" #"temperature" "random" "popsize"


  do_change_map_rates <- TRUE






  percentage_flow <- 10
  popchange_rate <- 1 # rate of population change
  q <- 0 # sort of rate of trait evolution.  it is a % of the total rate of popchange_rate


  #manual_speciation_events_timing <- c(50,300,600,1200,1500,1600,1700,1800,2000,2100) # it should come with lambda <-0.  manual_speciation_events_timing <- 0



  manual_speciation_events_timing <- 0






  # if == 0, it means that the simulation will run until reaching a second equilibrium after disturbance. If
  # it is > 0, it will be taken as the % of time that will be let run further, considering
  # the time it needed to reach the first equilibirum. To be used with options_type_scenario == "dynamic" / do_change_map_rates == TRUE
  time_percent_stop_after_first_equilibrium_and_disturbance <- 0

  #average_equi_richness <- readRDS("richness_equi_times.RDS")





  if(what_k == "uniform"){
    text_k_map <-  read.table(paste0("k_map_uniform.txt")) # "k_map_uniform_low" "k_map_uniform_high"
    text_k_map2 <-  text_k_map
  }

  if(what_temp == "uniform"){
    #adaptation_on <- FALSE
    map_temperature <- read.table(paste0("temperature_map_uniform.txt"))
    map_temperature2 <- map_temperature
  }

  if(what_k == "change"){
    text_k_map <-  read.table(paste0("k_map_gradient_strong.txt")) # "k_map_uniform_low" "k_map_uniform_high"
    text_k_map2 <-   read.table(paste0("k_map_gradient_mild.txt"))
  }
  if(what_k ==  "unchanged"){
    text_k_map <-   read.table(paste0("k_map_gradient_strong.txt")) # "k_map_uniform_low" "k_map_uniform_high"
    text_k_map2 <- text_k_map
  }

  if(what_temp == "change"){
    map_temperature <- read.table(paste0("temperature_map_gradient_strong.txt"))
    map_temperature2 <- read.table(paste0("temperature_map_gradient_mild.txt"))

  }

  if(what_temp == "unchanged"){
    #adaptation_on <- FALSE
    map_temperature <- read.table(paste0("temperature_map_gradient_mild.txt"))
    map_temperature2 <- map_temperature
  }

  the_v <-  5   # 5 or 5000 the smaller the value the steeper. A difference in temperature of 10 with v = 50, brings the fitness to 0.4
  map_elevation <- read.table("elevation_map.txt")

  #origin <- "highlands" # "highlands" ; "intermediate1"; "intermediate2" ; "lowlands"



  if(this_vicariance == "On"){
    vicariant_speciation <- TRUE
  } else {
    vicariant_speciation <- FALSE
  }



  #lambda <-0

  #popchange_rate <- 0
  if(this_lambda == "low"){
    lambda <- c(0.00005) #0.0005 #0.001 #0.00005
    #starting_one_spp <- FALSE

  }
  if(this_lambda == "mid"){

    lambda <- c(0.001) #0.0005 #0.001 #0.00005
  }
  if(this_lambda == "high"){
    lambda <- c(0.0001) #0.0005 #0.001 #0.00005

  }



  if(this_gene_flow == "low"){
    geneflow_rate <- 0.01 #c(0.01, 0.05,0.5,0.01, 0.05,0.5)[i]
  }

  if(this_gene_flow == "mid"){
    geneflow_rate <- 0.1  #c(0.01, 0.05,0.5,0.01, 0.05,0.5)[i]
  }
  if(this_gene_flow == "high"){
    geneflow_rate <- 0.1 #c(0.01, 0.05,0.5,0.01, 0.05,0.5)[i]
  }


  mutation_rate <- geneflow_rate/10


  v <- 5 # likely to be removed


  use_k <- TRUE
  restriction_parameter_k_traitdissimilarity <- 10

  show_richness_map <- "total_abundance" # "both" or "richness" or "total_abundance" or "none"
  show_richness_map <- "both" # "both" or "richness" or "total_abundance" or "none"

  starting_time <- 0
  # average_richness_this_case <- average_equi_richness[which(average_equi_richness$localK_case ==localK &
  #                                                         average_equi_richness$geneflow_case == this_gene_flow &
  #                                                         average_equi_richness$lambda_case == this_lambda &
  #                                                         average_equi_richness$type_map == this_map &
  #                                                         average_equi_richness$type_scenario == "static"),]$equi_richness

  #max_spp <- 10000 #round(as.numeric(average_richness_this_case) * 0.7)



  #time_start_simul <- Sys.time()

  #if(exists("simulation_raw") == FALSE){


  #sourceCpp("dynamic_model_cpp_gradient.cpp")

  map_k_vector <- transform_map_vector (text_k_map)
  map_k_vector2 <- transform_map_vector (text_k_map2)
  map_temperature_vector <- transform_map_vector(map_temperature)
  map_temperature_vector2 <- transform_map_vector(map_temperature2)
  map_elevation_vector <- transform_map_vector(map_elevation)




  #alleles_adaptation_coef2 <- c(13,9,8,3,2) # it should be same length than all_alleles. None should match the cell temperature!

  if(adaptation_on){
    alleles_adaptation_coef2 <- 1:25# it should be same length than all_alleles. None should match the cell temperature!
  } else {
    alleles_adaptation_coef2 <- rep(9,25) # it should be same length than all_alleles. None should match the cell temperature!


  }

  restiction_par <- 2  # this parameter might not be important



  succesful_simulation <- 0
  pres <- NULL
  avg <- NULL
  thesd <- NULL

  the_seed <- sample(1:9000,1)
  cat("                                      this is the seed: ",the_seed,"\n")
  cat("                   This is the simulation number", succesful_simulation + 1, "\n")



  the_seed <- 54
  set.seed(11)


  maximum_cycles <- 2000000000
  #maximum_cycles <- 750000
  #maximum_cycles <- 150000
  #######




  # store_entire_simulation <- list()
  # if(this_scenario == "dynamic" ){
  #   do_change_map_rates <- TRUE
  # }
  #
  # if(this_scenario == "static"){
  #   do_change_map_rates <- FALSE
  # }



  # simultation_name <- paste0("simul_K",what_k,"_T",what_temp,"_vicariance",this_vicariance)
  simulation_name <- paste0("simul_geneflow",this_gene_flow,"_vicariance",this_vicariance,"_",this_replicate)

  cat(simulation_name)
  if(starting_one_spp){
    number_spp <- 1
    all_x <- 17#3#3#10 #  initial population column
    all_y <- 17#3 #3#10 # initial population row
    all_traits <- 100 # initial trait value
    all_births <- 0
    all_IDs <- 1
    all_parents <- 0
    all_deaths <- 0
    all_ranges <- 1
    all_alleles <- rep(6,25) # allelic frequencies for allele a, b, c, d, and e.
    all_alleles_neutral <- rep(6,25)
    all_popsize <- 150 #population size for the first population. It should be the same than sum(all_alleles)

  } else {
    warning("needs adjusting after adding more alleles")
    number_spp <- 190
    all_x <- c(rep(2:20,10))#3#3#10 #  initial population column
    all_y <- c(rep(3,19), rep(5,19),rep(7,19),rep(9,19),rep(11,19),
               rep(13,19),rep(15,19),rep(17,19),rep(19,19),rep(21,19))#3 #3#10 # initial population row
    all_traits <- c(rep(100,190)) # initial trait value
    all_births <- c(rep(0,190))
    all_IDs <- c(1:190)
    all_parents <- c(rep(0,190))
    all_deaths <- c(rep(0,190))
    all_ranges <- c(rep(1,190))
    all_alleles <- rep(6,15) # allelic frequencies for allele a, b, c, d, and e.
    all_alleles_neutral <- rep(6,15)
    all_popsize <- c(rep(90,190))#population size for the first population. It should be the same than sum(all_alleles)
    lambda <- 0

  }
  store_entire_simulation <- list()



  #for(ij in 1:num_sim){

  simulation_raw <- do_simulation(map_elevation_vector,  map_k_vector,map_k_vector2,  map_temperature_vector, map_temperature_vector2,
                                  extirpation_depen,  colonization_depen_temperature, x_max,  y_max,  all_x,  all_y,
                                  all_IDs,  all_parents,  all_births,  all_deaths,
                                  all_traits,  all_ranges,  all_alleles, all_alleles_neutral,all_popsize,
                                  number_spp,  the_seed,
                                  mutation_rate, percentage_flow, geneflow_rate, popchange_rate, the_gammas,  the_mus,
                                  q,  lambda, species_trait_state,sd_normal_distribution_traitevol,mean_normal_distribution_traitevol,
                                  sd_normal_distribution_pop_change,growth_only,  starting_time,  simulated_time,max_spp,
                                  maximum_cycles,  use_k,  restiction_par, show_richness_map,  v,
                                  alleles_adaptation_coef2,
                                  do_change_map_rates,vicariant_speciation,speciation_rangesize_unlinked,
                                  colonization_rangesize_unlinked,time_percent_stop_after_first_equilibrium_and_disturbance, condition_to_stop,
                                  time_slices,manual_speciation_events_timing)


  cat(simulation_name)
  store_entire_simulation[[1]] <- simulation_raw
  #}

  saveRDS(store_entire_simulation,file = paste0(simulation_name,".RDS"))
  process_things(i)

}


# # to plot things
# library(ggpubr)
# library(ggplot2)
#
#
# simulation_raw_this_timeslice <- make_one_list_fromRaw(simulation_raw[[length(simulation_raw)]])
#
# #simulation_raw_this_timeslice <- make_one_list_fromRaw(simulation_raw[[20]])
#
#
# time_simulated_from_output <- simulation_raw_this_timeslice$total_time[[1]]
# list_species_from_cpp <- make_list_species_fromcpp_toR(simulation_raw_this_timeslice)
#
#
# richness_map <- make_richness_map(x_max,y_max,list_species_from_cpp)
# colnames(richness_map)[4] <- "Pop_size"
# class(richness_map)
#
# born_x <- all_x
#
# born_y <- all_y
# draw_map(richness_map,x_max,y_max,born_x,born_y,"richness")
#
# draw_map(richness_map,x_max,y_max,born_x,born_y,"pop_size")
# ks.test(rnorm(1000,0,1),rnorm(1000,4,1))
