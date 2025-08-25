build_tree <- function(time_simulated_from_output,
                       list_species_from_cpp,
                                        forced_crown_age = F){
  species_table <- list_species_from_cpp
    Ltable_to_do <- NULL
    for(i in 1:length(species_table)){
      take_this_sp <- species_table[[i]]
      #cat( take_this_sp$Birth, "\n")
      if(class(take_this_sp$per_species_population)[1] != "matrix"){
        pre_death <- time_simulated_from_output - take_this_sp$Death

      } else {
        pre_death <- take_this_sp$Death
      }
      Ltable_to_do <- data.frame(rbind(Ltable_to_do,cbind(birth = time_simulated_from_output - take_this_sp$Birth,
                                                          parent = take_this_sp$Parent,
                                                          id = take_this_sp$ID,
                                                         death = pre_death )))
    }
    Ltable_to_do[which(Ltable_to_do$death == 0),4] = -1
    phyloTree <- L2phylo(Ltable_to_do,dropextinct = T)
  return(phyloTree = phyloTree)
}









make_spp_table <- function(list_species_from_cpp){
  all_spp <- NULL
for(i in 1:length(list_species_from_cpp)){

  cat("doing species: ",i,"\n")
  this_species <- list_species_from_cpp[[i]]
  #### new part

    id_this_sp <- this_species$ID
    range <- this_species$RangeSize
    trait_value <- this_species$trait_state

    all_spp <- rbind(all_spp,
                     cbind(Sp_ID = rep(id_this_sp,range),
                     trait_value = rep(trait_value,range),
                     list_species_from_cpp[[i]]$per_species_population))
}

  return(all_spp)
}


make_one_list_fromRaw <- function(simulation_raw_this_timeslice){

  length(simulation_raw_this_timeslice)
  first_part <- simulation_raw_this_timeslice[[1]]
  second_part <- simulation_raw_this_timeslice[[2]]
  third_part <- simulation_raw_this_timeslice[[3]]
  fourth_part <- simulation_raw_this_timeslice[[4]]
  new_raw <- list()
  for(i in 1:length(first_part)){
    new_raw <- c(new_raw,first_part[i])
    #first_part[19][]
    #names(new_raw[i]) <- names( first_part[i])
  }
  for(i in 1:length(second_part)){
    new_raw <- c(new_raw,second_part[i])
    #first_part[19][]
    #names(new_raw[i]) <- names( first_part[i])
  }
  for(i in 1:length(third_part)){
    new_raw <- c(new_raw,third_part[i])
    #first_part[19][]
    #names(new_raw[i]) <- names( first_part[i])
  }
  for(i in 1:length(fourth_part)){
    new_raw <- c(new_raw,fourth_part[i])
    #first_part[19][]
    #names(new_raw[i]) <- names( first_part[i])
  }
  return(new_raw)
}
make_list_species_fromcpp_toR <- function(simulation_raw){
  list_species_from_cpp <- list()
  for(i in 1:length(simulation_raw$Species_Info$ID)){

    ID <- simulation_raw$Species_Info$ID[i]
    Parent <- simulation_raw$Species_Info$Parent[i]
    Birth <- simulation_raw$Species_Info$Birth[i]
    Birth_southmost <- simulation_raw$Species_Info$Birth_southmost[i]
    Birth_northmost <- simulation_raw$Species_Info$Birth_northmost[i]
    Death <- simulation_raw$Species_Info$Death[i]
    RangeSize <- simulation_raw$Species_Info$RangeSize[i]
    initial_rangesize <- simulation_raw$Species_Info$percentage_parental_range[i]
    swisscheese <- simulation_raw$Species_Info$swisscheese[i]
    trait_state <- simulation_raw$Species_Info$TraitValue[i]

    Total_pop_size <- simulation_raw$Species_Info$populationSize[i]
    geneflow_events <- simulation_raw$Species_Info$geneflow_events[i]
    saturation_grid_birth <- simulation_raw$Species_Info$saturation_grid_birth[i]

    change_northernmost <- simulation_raw$all_change_northernmost[[i]]
    time_change_northernmost <- simulation_raw$all_time_change_northernmost[[i]]
    change_southernmost <- simulation_raw$all_change_southernmost[[i]]
    time_change_southernmost <- simulation_raw$all_time_change_southernmost[[i]]
    northernnmost <- simulation_raw$Species_Info$northernnmost_locat[i]
    southernmost <- simulation_raw$Species_Info$southernmost_locat[i]
    expansion_failure_becauseK <- simulation_raw$Species_Info$expansion_failure_becauseK[i]
    total_change_southernmost <- simulation_raw$Species_Info$total_change_southernmost[i]
    total_change_northernmost <- simulation_raw$Species_Info$total_change_northernmost[i]
    if(Death == 0){ #still alive
      for(ii in 1:RangeSize){
        per_species_population <- cbind(X=simulation_raw$Distribution[[i]][[1]],   Y= simulation_raw$Distribution[[i]][[2]])
        per_species_population <- cbind(per_species_population,Pop_size=simulation_raw$popsize_perPop[[i]][[1]])

        Allele_A1 <- simulation_raw$Allele_A1[[i]][[1]]
        Allele_B1 <- simulation_raw$Allele_B1[[i]][[1]]
        Allele_C1 <- simulation_raw$Allele_C1[[i]][[1]]
        Allele_D1 <- simulation_raw$Allele_D1[[i]][[1]]
        Allele_E1 <- simulation_raw$Allele_E1[[i]][[1]]
        Allele_F1 <- simulation_raw$Allele_F1[[i]][[1]]
        Allele_G1 <- simulation_raw$Allele_G1[[i]][[1]]
        Allele_H1 <- simulation_raw$Allele_H1[[i]][[1]]
        Allele_I1 <- simulation_raw$Allele_I1[[i]][[1]]
        Allele_J1 <- simulation_raw$Allele_J1[[i]][[1]]
        Allele_K1 <- simulation_raw$Allele_K1[[i]][[1]]
        Allele_L1 <- simulation_raw$Allele_L1[[i]][[1]]
        Allele_M1 <- simulation_raw$Allele_M1[[i]][[1]]
        Allele_N1 <- simulation_raw$Allele_N1[[i]][[1]]
        Allele_O1 <- simulation_raw$Allele_O1[[i]][[1]]
        Allele_P1 <- simulation_raw$Allele_P1[[i]][[1]]
        Allele_Q1 <- simulation_raw$Allele_Q1[[i]][[1]]
        Allele_R1 <- simulation_raw$Allele_R1[[i]][[1]]
        Allele_S1 <- simulation_raw$Allele_S1[[i]][[1]]
        Allele_T1 <- simulation_raw$Allele_T1[[i]][[1]]
        Allele_U1 <- simulation_raw$Allele_U1[[i]][[1]]
        Allele_V1 <- simulation_raw$Allele_V1[[i]][[1]]
        Allele_W1 <- simulation_raw$Allele_W1[[i]][[1]]
        Allele_X1 <- simulation_raw$Allele_X1[[i]][[1]]
        Allele_Y1 <- simulation_raw$Allele_Y1[[i]][[1]]


        Allele_A2 <- simulation_raw$Allele_A2[[i]][[1]]
        Allele_B2 <- simulation_raw$Allele_B2[[i]][[1]]
        Allele_C2 <- simulation_raw$Allele_C2[[i]][[1]]
        Allele_D2 <- simulation_raw$Allele_D2[[i]][[1]]
        Allele_E2 <- simulation_raw$Allele_E2[[i]][[1]]
        Allele_F2 <- simulation_raw$Allele_F2[[i]][[1]]
        Allele_G2 <- simulation_raw$Allele_G2[[i]][[1]]
        Allele_H2 <- simulation_raw$Allele_H2[[i]][[1]]
        Allele_I2 <- simulation_raw$Allele_I2[[i]][[1]]
        Allele_J2 <- simulation_raw$Allele_J2[[i]][[1]]
        Allele_K2 <- simulation_raw$Allele_K2[[i]][[1]]
        Allele_L2 <- simulation_raw$Allele_L2[[i]][[1]]
        Allele_M2 <- simulation_raw$Allele_M2[[i]][[1]]
        Allele_N2 <- simulation_raw$Allele_N2[[i]][[1]]
        Allele_O2 <- simulation_raw$Allele_O2[[i]][[1]]
        Allele_P2 <- simulation_raw$Allele_P2[[i]][[1]]
        Allele_Q2 <- simulation_raw$Allele_Q2[[i]][[1]]
        Allele_R2 <- simulation_raw$Allele_R2[[i]][[1]]
        Allele_S2 <- simulation_raw$Allele_S2[[i]][[1]]
        Allele_T2 <- simulation_raw$Allele_T2[[i]][[1]]
        Allele_U2 <- simulation_raw$Allele_U2[[i]][[1]]
        Allele_V2 <- simulation_raw$Allele_V2[[i]][[1]]
        Allele_W2 <- simulation_raw$Allele_W2[[i]][[1]]
        Allele_X2 <- simulation_raw$Allele_X2[[i]][[1]]
        Allele_Y2 <- simulation_raw$Allele_Y2[[i]][[1]]




        # all_H_values <- NULL
        # all_H_selection <- NULL
        # all_H_neutral <- NULL
        # for (j in 1:length(Allele_A)){ # using H from Pegas package
        #   all_H_values <- c(all_H_values,H(c(Allele_A[j],Allele_B[j],Allele_C[j],Allele_D[j],Allele_E[j],
        #                                      Allele_V[j],Allele_W[j],Allele_X[j],Allele_Y[j],Allele_Z[j])))
        #   all_H_selection <- c(all_H_selection,H(c(Allele_A[j],Allele_B[j],Allele_C[j],Allele_D[j],Allele_E[j])))
        #   all_H_neutral <- c(all_H_neutral,H(c(Allele_V[j],Allele_W[j],Allele_X[j],Allele_Y[j],Allele_Z[j])))
        # }
        # per_species_population <- cbind(per_species_population,Allele_A,Allele_B,Allele_C,Allele_D,Allele_E,
        #                                Allele_V,Allele_W,Allele_X,Allele_Y,Allele_Z, H_value_bothLoci=all_H_values,
        #                                H_value_selection = all_H_selection, H_value_neutral = all_H_neutral )

        per_species_population <- cbind(per_species_population,Allele_A1,Allele_B1,Allele_C1,Allele_D1,Allele_E1,
                                        Allele_F1,Allele_G1,Allele_H1,Allele_I1,Allele_J1,
                                        Allele_K1,Allele_L1,Allele_M1,Allele_N1,Allele_O1,
                                        Allele_P1,
                                        Allele_Q1,
                                        Allele_R1,
                                        Allele_S1,
                                        Allele_T1,
                                        Allele_U1,
                                        Allele_V1,
                                        Allele_W1,
                                        Allele_X1,
                                        Allele_Y1,

                                        Allele_A2,Allele_B2,Allele_C2,Allele_D2,Allele_E2,
                                        Allele_F2,Allele_G2,Allele_H2,Allele_I2,Allele_J2,
                                        Allele_K2,Allele_L2,Allele_M2,Allele_N2,Allele_O2,
                                        Allele_P2,
                                        Allele_Q2,
                                        Allele_R2,
                                        Allele_S2,
                                        Allele_T2,
                                        Allele_U2,
                                        Allele_V2,
                                        Allele_W2,
                                        Allele_X2,
                                        Allele_Y2)

      }
    } else {
      per_species_population <- "extinct"
    }
    focal_species <- list(ID=ID,
                          Parent=Parent,
                          Birth=Birth,
                          Death=Death,
                          RangeSize=RangeSize,
                          initial_rangesize = initial_rangesize,
                          trait_state = trait_state,
                          swisscheese = swisscheese,
                          per_species_population=per_species_population,
                          saturation_grid_birth =saturation_grid_birth,
                          Birth_southmost = Birth_southmost,
                          Birth_northmost = Birth_northmost,
                          change_northernmost = change_northernmost,
                          time_change_northernmost = time_change_northernmost,
                          change_southernmost = change_southernmost,
                          time_change_southernmost = time_change_southernmost,
                          northernnmost = northernnmost,
                          southernmost =southernmost,
                          expansion_failure_becauseK = expansion_failure_becauseK,
                          total_change_southernmost = total_change_southernmost,
                          total_change_northernmost = total_change_northernmost,
                          geneflow_events = geneflow_events)

    list_species_from_cpp[[i]] <- focal_species
  }
  return(list_species_from_cpp)
}

#
# calc_community_matrix <- make_community_matrix(x_max,y_max,list_species_from_cpp)
#
# commun_mat <- calc_community_matrix$community_matrix
# calc_community_matrix$total_occupancy
#
# phylotree
#
#
#
#
# commun_mat_for_mpd <- commun_mat[,-c(1,2)]
# rownames(commun_mat_for_mpd) <- 1:nrow(commun_mat)
# phy_dist <- cophenetic(phylotree)
#
#
#
#
make_two_diversities_table <- function(commun_mat_pres_abs,allele_pres){

  neutral_allele_richness <- NULL

  selection_allele_richness <- NULL
  shannon_selection <- NULL
  shannon_neutral <- NULL
  species_richness <- NULL

  mean_shannon_species <- NULL

  for(i in 1:nrow(allele_pres)){

    this_cell <- allele_pres[i,]
    this_cell2 <- commun_mat_pres_abs[i,]

    mean_shannon_species <- c(mean_shannon_species,
                              vegan::diversity(this_cell2[3:ncol(commun_mat_pres_abs)],index="shannon"))

    #
    # shannon_selection <- c(shannon_selection,
    #                        vegan::diversity(c(as.numeric(this_cell[3:27])),index="shannon"))
    # shannon_neutral <- c(shannon_neutral,
    #                      vegan::diversity(c(as.numeric(this_cell[28:52])),index="shannon"))
    #
    # selection_allele_richness <- c(selection_allele_richness,
    #                                length(which(as.numeric(this_cell[3:27]) > 0)))
    # neutral_allele_richness <- c(neutral_allele_richness,
    #                              length(which(as.numeric(this_cell[28:52]) > 0)))
    species_richness <- c(species_richness,
                          sum(as.numeric(this_cell2[3:ncol(this_cell2)])))

  }

  div_table <- cbind(allele_pres,mean_shannon_species,species_richness )
  colnames(div_table) <- c("X","Y","mean_shannon_alleles","sum_alleles",
                           "mean_shannon_species","species_richness")

  div_table$mean_shannon_alleles[which(is.na(div_table$mean_shannon_alleles))] <- 0
  return(div_table)
}

compute_zscore_sar <- function(time_do,missing_replicates,total_replicates,total_escenarios){
  all_zed_sar <- NULL
  extinct_replicate_repli <- NULL
  extinct_replicate_escena <- NULL
  all_div_mat <- list()
  counter <- 0
  names_all_div_mat <- NULL
  for(i in 1:total_escenarios){
    cat("doing the escenario: ",i,"\n")
    options_what_k <- "uniform" #c("change","unchanged") # "uniform"
    options_what_temp <- "uniform" #c("change","unchanged") "uniform"
    options_vicariance_case <- c("On","Off")
    options_geneflow <- c("low","high")
    options_adaptation <- c(FALSE)

    all_combinations <- tidyr::expand_grid(options_what_k, options_what_temp, options_vicariance_case,options_adaptation,options_geneflow)

    what_k <- all_combinations$options_what_k [i]
    what_temp <-  all_combinations$options_what_temp [i]
    this_vicariance <- all_combinations$options_vicariance_case[i]
    what_adaptation <- all_combinations$options_adaptation [i]
    this_gene_flow <- all_combinations$options_geneflow [i]
    #this_replicate <- all_combinations$options_replicates[i]


    #all_replicates <- readRDS(paste0("simulation_",global_case,"_param_",param,".RDS"))
    #rm(all_replicates)

    all_these_avg_richness <- NULL
    for(ikj in 1:total_replicates){
      simulation_name <- paste0("simul_geneflow",this_gene_flow,"_vicariance",this_vicariance,"_",ikj)
      names_all_div_mat <- c(names_all_div_mat,simulation_name)
      counter <- counter + 1
      if(any(simulation_name == missing_replicates) == FALSE){
        all_replicates <- NULL
        try(all_replicates <- readRDS(paste0(simulation_name,".RDS")),
            silent=T)
        if(is.null(all_replicates) == FALSE){

          all_ready <- list()
          #all_replicates <- store_entire_simulation
          #for(ijk in 1:length(all_replicates)){
          ijk <- 1 # one replicate only
          #for(ijk in 1:1){


          #  for(ijk in 1:2){
          simulation_raw <- all_replicates[1][[1]]

          species_loca_richness <- NULL

          #present time is time_do = 0
          # time_do = 1 will do one time slice before present


          do_this_time <- length(simulation_raw) - time_do

          simulation_raw_this_timeslice <- simulation_raw[[do_this_time]] # l
          simulation_raw_this_timeslice <- make_one_list_fromRaw(simulation_raw_this_timeslice)

          time_simulated_from_output <- simulation_raw_this_timeslice$total_time[[1]]
          list_species_from_cpp <- make_list_species_fromcpp_toR(simulation_raw_this_timeslice)

          alive_vector <- NULL
          for(ijkl in 1:length(list_species_from_cpp)){

            if(list_species_from_cpp[[1]]$Death != 0){
              alive_vector <- c(alive_vector,"dead")

            } else {
              alive_vector <- c(alive_vector,"alive")
            }
          }

          if(all(alive_vector=="dead")){
            extinct_replicate_repli <- c(extinct_replicate_repli,ikj)
            extinct_replicate_escena <- c(extinct_replicate_escena,i)
            cat("whole clade died:",ikj, "\n")

          } else {
            calc_community_matrix <- make_community_matrix(x_max,y_max,list_species_from_cpp)



            commun_mat <- calc_community_matrix$community_matrix
            these_avg_richness <- NULL
            rangesize_from_here <- NULL
            commun_mat_pres_abs <- commun_mat
            for(ij in 1:nrow(commun_mat)){

              for(ijk in 1:ncol(commun_mat)){
                if(commun_mat[ij,ijk] > 0){
                  commun_mat_pres_abs[ij,ijk] <- 1
                } else{
                  commun_mat_pres_abs[ij,ijk] <- 0
                }
              }
            }

            commun_mat_pres_abs$X <- commun_mat$X
            commun_mat_pres_abs$Y <- commun_mat$Y
            sum(commun_mat_pres_abs$t1)
            head(commun_mat_pres_abs)

            allele_pres <- calc_community_matrix$alleles_pooled
            allele_pres <- cbind(commun_mat_pres_abs[,c(1,2)],allele_pres)
            all_div_mat[[counter]] <- make_two_diversities_table(commun_mat_pres_abs,allele_pres)

            table_sar <- make_table_sar (commun_mat_pres_abs,sampling_divided_by = 2)
            if(any(table_sar$richness_this_sampling_unit == 0)){
              table_sar$richness_this_sampling_unit[which(table_sar$richness_this_sampling_unit == 0)] <- 1
            }
            table_sar <- log(table_sar)

            zed_sar <- NA
            if(nrow(table_sar) >= 2){
              #print(cbind(table_gda_noNA_no0$radius,table_gda_noNA_no0$H_value_bothLoci))

              fit_sar <- NULL
              try(fit_sar <- sar_power(data = cbind((table_sar$radius),(table_sar$richness_this_sampling_unit)), grid_start = "partial"),
                  silent=T)
              if(is.null(fit_sar) == FALSE){
                try(zed_sar <- as.numeric(fit_sar$sigConf[2,1]),
                    silent=T)
              }

            } else {
              zed_sar <- NA
            }



          }



        }
      }
      if(is.numeric(zed_sar) == FALSE){
        stop("issue here")
      }
      all_zed_sar <- rbind(all_zed_sar, c(zed_sar,this_gene_flow,this_vicariance,ikj))
    }





  }
  colnames(all_zed_sar) <- c("zed_sar","gene_flow","vicariance","replicate")
  all_zed_sar <- as.data.frame(all_zed_sar)
  warning("make sure that all scenarios were analyze")
  return(list(all_zed_sar = all_zed_sar,extinct_replicate_repli = extinct_replicate_repli,
              extinct_replicate_escena = extinct_replicate_escena,
              all_div_mat = all_div_mat,
              names_all_div_mat = names_all_div_mat))
}


make_table_sar <- function(commun_mat_pres_abs,sampling_divided_by){
#  stop("this function is old and out of touch or it is perhaps intended for a single spp")
  do_this_table <- commun_mat_pres_abs
  #sampling_runs is not longer use as I set up the spatial centroid of the geographic distribution
  table_sar <- NULL
  geo_dista <- as.matrix(dist(cbind(do_this_table[,1], do_this_table[,2])))
  # sampled_table_and_centroid <- prepare_sampled_table_and_find_centroid(do_this_table,num_sampled_cells=1,only_edges = T,sample_it = T, strict_geographic = T)
  centroid_x <- 20
  centroid_y <- 20
  centroid_pop <- which(do_this_table$X == centroid_x &  do_this_table$Y == centroid_y)
  initial_population <- centroid_pop

  dist_to_initial <- geo_dista[initial_population,]
  #dist_to_initial <- dist_to_initial[-which(dist_to_initial == 0)]
  dist_to_initial_sorted <- sort(unique(round(dist_to_initial)))

  table_sar_this_sampling <- NULL
  accumu_to_sampled_populations <- initial_population
  for(i in 2:length(dist_to_initial_sorted)){
    popul_within_radius <- which(dist_to_initial_sorted[i] >= dist_to_initial)
    if(length(popul_within_radius) == 1){
      popul_within_radius <-  which(dist_to_initial_sorted[i + 1] >= dist_to_initial)
    }


    amount_to_sample <- ceiling(length(popul_within_radius)/sampling_divided_by)
    if(amount_to_sample == 1){
      amount_to_sample <- 2
    }
    #to_sampled_populations <- as.numeric(names(sample(popul_within_radius,amount_to_sample,replace=FALSE)))
    if(length(popul_within_radius) == 1){ # to avoid undesiarable behaviour of sample function when pool length == 1. Though because i is never 1, pool is never < 2
      to_sampled_populations <- as.numeric(popul_within_radius)
    } else {
      to_sampled_populations <- as.numeric((sample(popul_within_radius,amount_to_sample,replace=FALSE)))
    }
    to_sampled_populations
    accumu_to_sampled_populations <- c(accumu_to_sampled_populations,to_sampled_populations)
    #cat("accumu_to_sampled_populations: ", accumu_to_sampled_populations,"\n")
    sampled_populations <- as.matrix(do_this_table[unique(accumu_to_sampled_populations),],ncol=53)
    # print(sampled_populations)
    colnames(sampled_populations) <- colnames(do_this_table)
    pooled_populations <- sampled_populations
    pooled_populations2 <- pooled_populations[,-c(1,2)]
    if(nrow(pooled_populations2) > 1){
      consolidated <- colSums(pooled_populations2)
    } else {
      consolidated <- (pooled_populations2)
    }

    #colSums(pooled_populations)


    #consolidated <- colSums(pooled_populations2)
    consolidated[which(consolidated > 0)] <- 1
    richness_this_sampling_unit <- sum(consolidated)
    cat("richness_this_sampling_unit", richness_this_sampling_unit,"\n")

    the_radius <- 3.1416 * (dist_to_initial_sorted[i] ^ 2)

    table_sar_this_sampling <- as.data.frame(rbind(table_sar_this_sampling,
                                                   cbind(radius = the_radius,
                                                         richness_this_sampling_unit = richness_this_sampling_unit
                                                   )))

  }
  table_gda <- rbind(table_sar, table_sar_this_sampling)



  return(table_gda)
}

local_coexistence <- function(time_do,missing_replicates,total_replicates,total_escenarios){
  all_all_these_avg_richness <- NULL
  extinct_replicate_repli <- NULL
  extinct_replicate_escena <- NULL
  for(i in 1:total_escenarios){
    cat("doing the escenario: ",i,"\n")
    options_what_k <- "uniform" #c("change","unchanged") # "uniform"
    options_what_temp <- "uniform" #c("change","unchanged") "uniform"
    options_vicariance_case <- c("On","Off")
    options_geneflow <- c("low","high")
    options_adaptation <- c(FALSE)

    all_combinations <- tidyr::expand_grid(options_what_k, options_what_temp, options_vicariance_case,options_adaptation,options_geneflow)

    what_k <- all_combinations$options_what_k [i]
    what_temp <-  all_combinations$options_what_temp [i]
    this_vicariance <- all_combinations$options_vicariance_case[i]
    what_adaptation <- all_combinations$options_adaptation [i]
    this_gene_flow <- all_combinations$options_geneflow [i]
    #this_replicate <- all_combinations$options_replicates[i]


    #all_replicates <- readRDS(paste0("simulation_",global_case,"_param_",param,".RDS"))
    #rm(all_replicates)

    all_these_avg_richness <- NULL
    for(ikj in 1:total_replicates){
      simulation_name <- paste0("simul_geneflow",this_gene_flow,"_vicariance",this_vicariance,"_",ikj)

      if(any(simulation_name == missing_replicates) == FALSE){
        all_replicates <- NULL
        try(all_replicates <- readRDS(paste0(simulation_name,".RDS")),
            silent=T)
        if(is.null(all_replicates) == FALSE){

          all_ready <- list()
          #all_replicates <- store_entire_simulation
          #for(ijk in 1:length(all_replicates)){
          ijk <- 1 # one replicate only
          #for(ijk in 1:1){


          #  for(ijk in 1:2){
          simulation_raw <- all_replicates[1][[1]]

          species_loca_richness <- NULL

          #present time is time_do = 0
          # time_do = 1 will do one time slice before present


          do_this_time <- length(simulation_raw) - time_do

          simulation_raw_this_timeslice <- simulation_raw[[do_this_time]] # l
          simulation_raw_this_timeslice <- make_one_list_fromRaw(simulation_raw_this_timeslice)

          time_simulated_from_output <- simulation_raw_this_timeslice$total_time[[1]]
          list_species_from_cpp <- make_list_species_fromcpp_toR(simulation_raw_this_timeslice)

          alive_vector <- NULL
          for(ijkl in 1:length(list_species_from_cpp)){

            if(list_species_from_cpp[[1]]$Death != 0){
              alive_vector <- c(alive_vector,"dead")

            } else {
              alive_vector <- c(alive_vector,"alive")
            }
          }

          if(all(alive_vector=="dead")){
            extinct_replicate_repli <- c(extinct_replicate_repli,ikj)
            extinct_replicate_escena <- c(extinct_replicate_escena,i)
            cat("whole clade died:",ikj, "\n")

          } else {
            calc_community_matrix <- make_community_matrix(x_max,y_max,list_species_from_cpp)



            commun_mat <- calc_community_matrix$community_matrix
            these_avg_richness <- NULL
            rangesize_from_here <- NULL
            commun_mat_pres_abs <- commun_mat
            for(ij in 1:nrow(commun_mat)){

              for(ijk in 1:ncol(commun_mat)){
                if(commun_mat[ij,ijk] > 0){
                  commun_mat_pres_abs[ij,ijk] <- 1
                } else{
                  commun_mat_pres_abs[ij,ijk] <- 0
                }
              }
            }

            local_rich <- rowSums(commun_mat_pres_abs[,3:ncol(commun_mat_pres_abs)])
            these_avg_richness <- NULL
            rangesize_from_here <- NULL
            for(ii in 3:ncol(commun_mat)){

              cells_sp_present <- which(commun_mat[,ii]>0)
              coexisting <- as.numeric(local_rich[cells_sp_present])
              rangesize_from_here <- c(rangesize_from_here,length(cells_sp_present))
              these_avg_richness <- c(these_avg_richness,mean(coexisting))

            }
            all_these_avg_richness <- rbind(all_these_avg_richness, cbind(colnames(commun_mat)[3:ncol(commun_mat)],these_avg_richness,rangesize_from_here,repli_fromhere=rep(ikj,length(these_avg_richness))))

          }



        }
      }
    }





    all_all_these_avg_richness <- rbind(all_all_these_avg_richness,cbind(all_these_avg_richness,escena_fromhere=rep(i,nrow(all_these_avg_richness))))
  }
  all_all_these_avg_richness <- as.data.frame(all_all_these_avg_richness)
  warning("output table has no extinct species")
  return(list(all_all_these_avg_richness = all_all_these_avg_richness,extinct_replicate_repli = extinct_replicate_repli,
              extinct_replicate_escena = extinct_replicate_escena))
}




wrapper_make_gda <- function(do_centroid_pop,a_species_table,range_size_threshold,sampling_runs,sampling_divided_by){
  do_this_table <- a_species_table
  if(nrow(do_this_table) > range_size_threshold){
    computed_table <- make_table_gda(do_centroid_pop,do_this_table,sampling_runs,sampling_divided_by)
    table_gda <- computed_table$table_gda
    sampled_table_and_centroid <- computed_table$sampled_table_and_centroid
    table_gda$fst_bothLoci <- table_gda$fst_bothLoci + 1 # because the power function does not work well with negative values, which are bound to exist when doing log of a small number
    table_gda <- log(table_gda)

    table_gda_noNA_no0 <- NULL
    if(any(is.na(table_gda$H_value_bothLoci))){
      table_gda_noNA_no0 <- table_gda[-which(is.na(table_gda$H_value_bothLoci)),]
    } else {
      table_gda_noNA_no0 <- table_gda
    }

    if(any(table_gda_noNA_no0$H_value_bothLoci == 0)){
      table_gda_noNA_no0 <- table_gda_noNA_no0[-which(table_gda_noNA_no0$H_value_bothLoci == 0),]

    }

    #
    # ggplot(table_gda_noNA_no0, aes(x=log(radius), y=log(H_value_bothLoci)  )) +
    #   geom_point(size = 3)+
    #   # xlim(-6, 1)+
    #   # ylim(-7.5, -5)+
    #   geom_smooth(method=lm)

    zed_H <- NA
    if(nrow(table_gda_noNA_no0) >= 2){
      #print(cbind(table_gda_noNA_no0$radius,table_gda_noNA_no0$H_value_bothLoci))

      fit_H <- NULL
      try(fit_H <- sar_power(data = cbind((table_gda_noNA_no0$radius),(table_gda_noNA_no0$H_value_bothLoci)), grid_start = "partial"),
          silent=T)
      if(is.null(fit_H) == FALSE){
        try(zed_H <- as.numeric(fit_H$sigConf[2,1]),
            silent=T)
      }

    } else {
      zed_H <- NA
    }


    table_gda_noNA_no0 <- NULL
    if(any(is.na(table_gda$fst_bothLoci))){
      table_gda_noNA_no0 <- table_gda[-which(is.na(table_gda$fst_bothLoci)),]
    } else {
      table_gda_noNA_no0 <- table_gda
    }

    if(any(table_gda_noNA_no0$fst_bothLoci <= 0)){
      table_gda_noNA_no0 <- table_gda_noNA_no0[-which(table_gda_noNA_no0$fst_bothLoci <= 0),]

    }
    zed_fst <- NA
    if(nrow(table_gda_noNA_no0) >= 2){
      #print(cbind(table_gda_noNA_no0$radius,table_gda_noNA_no0$fst_bothLoci))

      fit_fst <- NULL
      try(fit_fst <- sar_power(data = cbind((table_gda_noNA_no0$radius),(table_gda_noNA_no0$fst_bothLoci)), grid_start = "partial") ,
          silent=T)

      if(is.null(fit_fst) == FALSE){
        try(zed_fst <- as.numeric(fit_fst$sigConf[2,1]),
            silent=T)
      }

    } else {
      zed_fst <- NA
    }
    #summary(fit)

    zed_ind <- NA
    fit_ind <- NULL
    try(fit_ind <- sar_power(data = cbind((table_gda$radius),(table_gda$ind_here)), grid_start = "partial") ,
        silent=T)

    if(is.null(fit_ind) == FALSE){
      try(zed_ind <- as.numeric(fit_ind$sigConf[2,1]),
          silent=T)
    }


    zed_alleles <- NA
    fit_alleles <- NULL
    try(fit_alleles <- sar_power(data = cbind((table_gda$radius),(table_gda$alleles_bothLoci)), grid_start = "partial") ,
        silent=T)
    if(is.null(fit_alleles) == FALSE){
      try(zed_alleles <- as.numeric(fit_alleles$sigConf[2,1]),
          silent=T)

    }


  } else {
    cat ("small-ranged species for z-score analysis \n")
    zed_H <- NA
    zed_fst <- NA
    zed_ind <- NA
    zed_alleles <- NA
    sampled_table_and_centroid <- NA
  }
  cat("zmar_H: ",zed_H,"zmar_fst: ",zed_fst, "\n")
  return(list(zed_fst = zed_fst,
              zed_H = zed_H,
              zed_ind = zed_ind,
              zed_alleles = zed_alleles,
              sampled_table_and_centroid = sampled_table_and_centroid))
}

keep_long_live_sp <- function(use_this_table,min_num_time_points = 1){
  long_live_sp_table <- NULL


  for(i in 1:length(unique(use_this_table$escena))){


    use_this_table2 <- NULL
    use_this_table2 <-    use_this_table[which(use_this_table$escena == i),]

    use_this_table3 <- NULL
    for(ii in 1:length(unique(use_this_table2$repli))){
      use_this_table3 <- use_this_table2[which(use_this_table2$repli == ii),]
      for(ikk in 1:length(unique(use_this_table3$id))){
        sp_tosee <- unique(use_this_table3$id)[ikk]
        sp_subtable <- use_this_table3[which(use_this_table3$id == sp_tosee),]
        if(all(sp_subtable$alive == "alive")){
          finalsurvival <- rep("extant",nrow(sp_subtable))
        } else {
          finalsurvival <- rep("extinct",nrow(sp_subtable))
        }
        cat("species ",sp_tosee ," has ",nrow(sp_subtable)," timepoints \n")
        if(nrow(sp_subtable) >= min_num_time_points){

          long_live_sp_table <- rbind(long_live_sp_table,
                                      cbind(use_this_table3[which(use_this_table3$id == sp_tosee),],finalsurvival))


        }
      }


    }

  }
  return(long_live_sp_table)
}



avg_phylo_dist_focalspp_per_cell <- function(commun_mat,phylotree){
  cat("   computing focal phylo distance \n")
  #commun_mat2 <- commun_mat


  commun_mat_for_mpd <- commun_mat[,-c(1,2)]
  rownames(commun_mat_for_mpd) <- 1:nrow(commun_mat)
  phy_dist <- cophenetic(phylotree)
  average_focal_mat1 <- cbind(x=commun_mat$X,y=commun_mat$Y)
  rownames(average_focal_mat1) <- rownames(commun_mat)
  vector_avg_mpd <- NULL

  per_cell_richness <- NULL

  for(i in 1:nrow(commun_mat_for_mpd)){

    per_cell_richness <- c(per_cell_richness,length(which(commun_mat_for_mpd[i,] > 0)))
  }

  average_focal_mat1 <- cbind(average_focal_mat1,per_cell_richness)

  for(i in 1:ncol(commun_mat_for_mpd)){
    species_focal <- colnames(commun_mat_for_mpd)[i]
    cat("species: ",species_focal, "\n")
    average_focal_mat1 <- cbind(average_focal_mat1,this_sp = NA)
    colnames(average_focal_mat1)[ncol(average_focal_mat1)] <- species_focal
    #commun_mat_for_mpd <- commun_mat[,-c(1,2)]
    commun_mat_for_mpd_prunned_to_focal <- commun_mat_for_mpd[-which(commun_mat_for_mpd[,i] == 0),]
    #### NRIfocal calculation
    #Standardized effect size of mpd

    table_metric_cell <- as.data.frame(cbind(cell_id=rownames(commun_mat_for_mpd_prunned_to_focal),metric=focal_nri(commun_mat_for_mpd_prunned_to_focal, phy_dist,species_focal, null.model = "richness",
                                                                                                                    abundance.weighted = TRUE, runs = 99, iterations = 100)$mpd.obs.z))

    for(ii in 1:nrow(table_metric_cell)){

      average_focal_mat1[which(rownames(average_focal_mat1) == table_metric_cell$cell_id[ii]),ncol(average_focal_mat1)] <- as.numeric(table_metric_cell$metric[ii])
    }


  }
  #vector_avg_mpd[which(is.na(vector_avg_mpd))] <- 0
  return(average_focal_mat1)
}

avg_phylo_dist_per_cell <- function(commun_mat,phylotree){
  cat("   computing focal phylo distance \n")
  commun_mat_for_mpd <- commun_mat[,-c(1,2)]
  rownames(commun_mat_for_mpd) <- 1:nrow(commun_mat)
  phy_dist <- cophenetic(phylotree)

  uio <- ses.mpd(commun_mat_for_mpd, phy_dist,null.model="taxa.labels")


  table_avg_phylo_cell <- cbind(commun_mat$X,commun_mat$Y,uio$mpd.obs.z)


  return(table_avg_phylo_cell)
}

compute_pic_table <- function(simulation_name,
                              variables_to_do,
                              variable_to_plot_ontree = NULL,
                              filtering_rangeSize,
                              log_transform = FALSE){

  all_pic_table <- NULL
  counter <- 0
  no_replicates <- NULL
  raw_table <- NULL



    this_rep <- as.data.frame(readRDS(paste0("processed",simulation_name,".RDS")))

    time_do <- 0
    re <- this_rep[which(this_rep$time_slice == sort(unique(this_rep$time_slice),decreasing=T)[ time_do + 1]),]

    re <- re[which(re$survival == "extant"),]
    all_raw_data <- readRDS(paste0(simulation_name,".RDS"))[1][[1]]



    simulation_raw_this_timeslice <- all_raw_data[[length(all_raw_data)]]
    simulation_raw_this_timeslice <- make_one_list_fromRaw(simulation_raw_this_timeslice)

    time_simulated_from_output <- simulation_raw_this_timeslice$total_time[[1]]
    list_species_from_cpp <- make_list_species_fromcpp_toR(simulation_raw_this_timeslice)


      tree <-  build_tree_and_sister_pairs(time_simulated_from_output,list_species_from_cpp)$phyloTree
      clade_age <- time_simulated_from_output

      re <- cbind(re,id_full = paste0("t",re$id))

      if(length(tree$tip.label) != nrow(re)){
        stop("issue reading trees")
      }
      this_replicate_filtered <- re[which(re$Range > filtering_rangeSize),]

      tips_to_drop <- NULL
      for(jji in 1:length(tree$tip.label)){
        if(all(this_replicate_filtered$id_full != tree$tip.label[jji])){
          tips_to_drop <- c(tips_to_drop,tree$tip.label[jji])
        }
      }

      tree_filtered <- drop.tip(tree,tips_to_drop)
      spp_richness_from_tree <- NULL
      spp_richness_from_tree <- length(tree_filtered$tip.label)
      pic_table <- NULL
      counter <- counter + 1
      observ_to_remove_na <- NULL
      for( iijj in 1:length(variables_to_do)){
        to_this_variable <- variables_to_do[iijj]
        variable_raw <- this_replicate_filtered[,which(colnames(this_replicate_filtered) == to_this_variable)]
        names(variable_raw) <- this_replicate_filtered$id_full
        if(log_transform){

          variable_raw <- log(variable_raw)
        }

        #remove those observations with NAs
        if(any(is.na(variable_raw))){
          observ_to_remove_na <- c(observ_to_remove_na,which(is.na(variable_raw)))
        }

        if(any(is.infinite(variable_raw))){
          observ_to_remove_na <- c(observ_to_remove_na,which(is.infinite(variable_raw)))
        }

      }
      for( iijj in 1:length(variables_to_do)){
        to_this_variable <- variables_to_do[iijj]
        variable_raw <- this_replicate_filtered[,which(colnames(this_replicate_filtered) == to_this_variable)]
        names(variable_raw) <- this_replicate_filtered$id_full

        if(log_transform){

          variable_raw <- log(variable_raw)
        }
        if(is.null(observ_to_remove_na) == FALSE){
          variable_raw_nafree <- variable_raw[-unique(observ_to_remove_na)]
          tree_filtered_nafree <- drop.tip(tree_filtered, unique(names(observ_to_remove_na)))
        } else {
          variable_raw_nafree <- variable_raw
          tree_filtered_nafree <- tree_filtered
        }

        if(length(tree_filtered_nafree$tip.label) > 3){
          variable_pic <- pic(variable_raw_nafree, tree_filtered_nafree)
          variable_pic <- as.numeric(variable_pic)
          if(is.null(variable_to_plot_ontree) == FALSE){
            if(to_this_variable == variable_to_plot_ontree){
              dotTree(tree_filtered_nafree,variable_raw_nafree,length=10,ftype="i")
            }
          }


        } else {
          variable_pic <- NA
          variable_raw_nafree <- NA
        }

        raw_table <- cbind(raw_table,as.numeric(variable_raw_nafree))
        pic_table <- cbind(pic_table,variable_pic)

      }
  return(list(pic_table = pic_table,
              raw_table = raw_table ))
}


do_r2_GLMM <- function(genetic_diver){
  r2_GLMM <- NULL
  for(i in 1:length(genetic_diver)){
    r2_GLMM <- c(r2_GLMM,
                 r.squaredGLMM(genetic_diver[[i]])[1,1])
  }
  r2_GLMM <- as.numeric(r2_GLMM)
  r2_GLMM_table <- data.frame(model=names(genetic_diver),
                              r2_GLMM = r2_GLMM)
  r2_GLMM_table <- r2_GLMM_table[order(r2_GLMM_table$r2_GLMM,decreasing = TRUE),]
  return(r2_GLMM_table)
}


do_r2_beta <- function(genetic_diver){
  r2_beta <- NULL
  for(i in 1:length(genetic_diver)){
    r2_beta <- c(r2_beta,
                 r2beta(genetic_diver[[i]],method='sgv',partial=FALSE)$Rsq)
  }

  r2_beta_table <- data.frame(model=names(genetic_diver),
                              r2_beta = r2_beta)
  r2_beta_table <- r2_beta_table[order(r2_beta_table$r2_beta,decreasing = TRUE),]
  return(r2_beta_table)
}

do_aic_table <- function(genetic_diver){

  Aics <- NULL
  for(i in 1:length(genetic_diver)){
    Aics <- c(Aics,AIC(genetic_diver[[i]]))

  }
  calculated_AICweights <- ICbweights2(Aics)
  aic_table <- data.frame(model=names(genetic_diver),
                          AIC = Aics,
                          AICweights = calculated_AICweights)
  aic_table <- aic_table[order(aic_table$AICweights,decreasing = TRUE),]
  return(aic_table)
}

genet_struct_evol_rate_calc <- function(phylotree,species_table_summary,trait_to_do = NULL){

  sppname_no_t <- NULL
  for(i in 1:length(phylotree$tip.label)){

    #sppname_no_t <- c(sppname_no_t,as.numeric(str_split_1(phylotree$tip.label[i], "t")[2]))
    sppname_no_t <- c(sppname_no_t,as.numeric(str_split(phylotree$tip.label[i], "t")[[1]][2]))
  }

  species_table_summary_extant <- species_table_summary[which(species_table_summary$alive == "alive"),]
  table_tip_and_spp_name <- data.frame(tip_label = 1:length(species_table_summary_extant$id),
                                       full_sppname=phylotree$tip.label,
                                       sppname = sppname_no_t)# in trees the tip number is different from tiplabel


  the_trait <- species_table_summary[,which(colnames(species_table_summary)==trait_to_do)]
  the_trait <- as.data.frame(cbind(id=species_table_summary$id,the_trait))


  traits_order_according_tree <- NULL
  for(i in 1:nrow(table_tip_and_spp_name)){
    traits_order_according_tree <- c(traits_order_according_tree,
                                     the_trait[which(table_tip_and_spp_name$sppname[i] == the_trait$id),]$the_trait)
  }
  traits_order_according_tree <- as.matrix(traits_order_according_tree,ncol=1)
  rownames(traits_order_according_tree) <- phylotree$tip.label

  # using evorates package, bayesian model
  # fit_evorates <- fit.evorates(tree = phylotree, trait.data = traits_order_according_tree, chains = 1)
  # fit_evorates_toframe <- get.R(fit_evorates,select=1:length(phylotree$tip.label),type="means",simplify=T)
  # rates_tips <- as.data.frame(fit_evorates_toframe) # these rates come with an edge index
  # end of it

  tips_to_remove <- NULL
  traits_order_according_tree_NAremoved <- NULL
  for(ii in 1:nrow(traits_order_according_tree)){
    if(is.na(traits_order_according_tree[ii,])){
      tips_to_remove <- c(tips_to_remove,names(traits_order_according_tree[ii,]))
    } else {
      traits_order_according_tree_NAremoved <- c(traits_order_according_tree_NAremoved,
                                                 traits_order_according_tree[ii,])
    }

  }
  phylotree2 <- ape::drop.tip(phylotree,tips_to_remove)

  # fitBM <- multirateBM(phylotree2,traits_order_according_tree_NAremoved,
  #                      optim="L-BFGS-B", lambda=0.01)
  #

  ## print and plot the results
  # print(fitBM)
  # plot(fitBM,ftype="i",fsize=0.8,lwd=6,
  #      outline=TRUE)


  rates_with_added_NAs <- NULL
  for(ii in 1:length(phylotree$tip.label)){
    #
    # if(any(which(phylotree$tip.label[ii] == names(fitBM$sig)))){
    #   rates_with_added_NAs <- c(rates_with_added_NAs,
    #                             as.numeric(fitBM$sig[which(phylotree$tip.label[ii] == names(fitBM$sig))]))
    # } else {
    rates_with_added_NAs <- c(rates_with_added_NAs,NA)
    #}

  }


  table_tip_and_spp_name_rate <- cbind(table_tip_and_spp_name,
                                       evorates = rates_with_added_NAs)

  vector_terminaltip_length <- setNames(phylotree$edge.length[sapply(1:length(phylotree$tip.label),
                                                                     function(x,y) which (y==x),y=phylotree$edge[,2])],phylotree$tip.label)
  if(any(names(vector_terminaltip_length) != table_tip_and_spp_name_rate$full_sppname)){
    stop("problem with tip labels")

  }
  if(any(names(vector_terminaltip_length) != rownames(traits_order_according_tree))){
    stop("problem with tip labels")

  }
  table_tip_and_spp_name_rate <- cbind(table_tip_and_spp_name_rate,
                                       terminaltip_length = vector_terminaltip_length,
                                       traitvalue = traits_order_according_tree)
  export_rates_order_summary_table <- NULL
  export_for_branch_controlled <- NULL
  export_branch_length_controlled_area <- NULL
  export_branch_length_controlled_tiplength <- NULL
  export_branch_length <- NULL
  for(i in 1:nrow(species_table_summary)){

    take_this_sp <- table_tip_and_spp_name_rate[which(species_table_summary$id[i] == table_tip_and_spp_name_rate$sppname),]
    range_size <- species_table_summary[which(species_table_summary$id == take_this_sp$sppname),]$Range
    TotInd_this_sp <- species_table_summary[which(species_table_summary$id == take_this_sp$sppname),]$TotInd
    export_rates_order_summary_table <- c(export_rates_order_summary_table,
                                          take_this_sp$evorates)

    export_for_branch_controlled <- c(export_for_branch_controlled,
                                      take_this_sp$traitvalue/TotInd_this_sp)

    export_branch_length_controlled_area <-  c(export_branch_length_controlled_area,
                                               (take_this_sp$traitvalue/range_size))

    export_branch_length_controlled_tiplength <-  c(export_branch_length_controlled_tiplength,
                                                    (take_this_sp$traitvalue/take_this_sp$terminaltip_length))
    export_branch_length <- c(export_branch_length,
                              take_this_sp$terminaltip_length)
  }
  return(list(modelled_rates = export_rates_order_summary_table,
              branch_length_controlled = export_for_branch_controlled,
              branch_length_controlled_area = export_branch_length_controlled_area,
              branch_length_controlled_tiplength = export_branch_length_controlled_tiplength,
              branch_length = export_branch_length))
}


classify_into_genera <- function(species_table_summary,phylotree){


  sppname_no_t <- NULL
  for(i in 1:length(phylotree$tip.label)){

    #sppname_no_t <- c(sppname_no_t,as.numeric(str_split_1(phylotree$tip.label[i], "t")[2]))
    sppname_no_t <- c(sppname_no_t,as.numeric(str_split(phylotree$tip.label[i], "t")[[1]][2]))
  }
  species_table_summary_extant <- species_table_summary[which(species_table_summary$alive == "alive"),]

  table_tip_and_spp_name <- data.frame(tip_label = 1:length(species_table_summary_extant$id),
                                       full_sppname=phylotree$tip.label,
                                       sppname = sppname_no_t)# in trees the tip number is different from tiplabel

  table_tip_and_spp_name <- table_tip_and_spp_name[order(table_tip_and_spp_name$sppname,decreasing = FALSE),]



  crown_lineage <- phylotree$edge[1,][1]
  first_split <- phylotree$edge[which(phylotree$edge[,1] ==  crown_lineage),2]

  familyA <- Descendants(phylotree,first_split[1],"tips")[[1]]
  familyB <- Descendants(phylotree,first_split[2],"tips")[[1]]
  if(length(familyA) > 1){
    split_familyA <- phylotree$edge[which(phylotree$edge[,1] ==  first_split[1]),2]
    genus_A_a <-  Descendants(phylotree,split_familyA[1],"tips")[[1]]
    genus_A_b <-  Descendants(phylotree,split_familyA[2],"tips")[[1]]
  } else { # it is a terminal branch then
    genus_A_a <- NULL
    familyA <- NULL
    genus_B_b <- NULL

  }
  if(length(familyB) > 1){
    split_familyB <- phylotree$edge[which(phylotree$edge[,1] ==  first_split[2]),2]

    genus_B_a <-  Descendants(phylotree,split_familyB[1],"tips")[[1]]
    genus_B_b <-  Descendants(phylotree,split_familyB[2],"tips")[[1]]
  } else {# it is a terminal branch then
    genus_B_a <- NULL
    familyB <- NULL
    genus_B_b <- NULL

  }

  family <- NULL
  genus <- NULL
  for(i in 1:nrow(table_tip_and_spp_name)){
    it_is_pending_genus <- TRUE
    it_is_pending_family <- TRUE
    if(any(table_tip_and_spp_name$tip_label[i] == genus_A_a)){
      genus <- c(genus,"genus_A_a")
      it_is_pending_genus <- FALSE
    }
    if(any(table_tip_and_spp_name$tip_label[i] == genus_A_b)){
      genus <- c(genus,"genus_A_b")
      it_is_pending_genus <- FALSE
    }
    if(any(table_tip_and_spp_name$tip_label[i] == genus_B_a)){
      genus <- c(genus,"genus_B_a")
      it_is_pending_genus <- FALSE
    }
    if(any(table_tip_and_spp_name$tip_label[i] == genus_B_b)){
      genus <- c(genus,"genus_B_b")
      it_is_pending_genus <- FALSE
    }

    if(any(table_tip_and_spp_name$tip_label[i] == familyA)){
      family <- c(family,"familyA")
      it_is_pending_family <- FALSE
    }
    if(any(table_tip_and_spp_name$tip_label[i] == familyB)){
      family <- c(family,"familyB")
      it_is_pending_family <- FALSE
    }


    if(is.null(genus_A_a) && it_is_pending_genus){
      genus <- c(genus,"genus_A_a")
      cat("jere genus_A_a\n")
      cat(i)
    }

    if(is.null(genus_B_a) && it_is_pending_genus){
      genus <- c(genus,"genus_B_a")
      cat("genus_B_a \n")
      cat(i)
    }

    if(is.null(familyA) && it_is_pending_family){
      family <- c(family,"familyA")
      cat("jere familyA \n")
      cat(i)
    }
    if(is.null(familyB) && it_is_pending_family){
      family <- c(family,"familyB")
      cat("jere familyB \n")
      cat(i)
    }
    #cat(i,"length:",length(family), family[length(family)],"\n")

  }

  return(list(family = family,
              genus = genus))
}

ICbweights2 <- function(IC){
  bestmodelIC <- min(IC)
  weights <- exp(-0.5*(IC-bestmodelIC))
  weights <- weights/sum(weights)
  return(weights)
}
avg_phylo_dist_spp <- function(commun_mat,phylotree){
  cat("   computing focal phylo distance \n")
  commun_mat_for_mpd <- commun_mat[,-c(1,2)]
  rownames(commun_mat_for_mpd) <- 1:nrow(commun_mat)
  phy_dist <- cophenetic(phylotree)

  vector_avg_mpd <- NULL
  for(i in 1:ncol(commun_mat_for_mpd)){
    species_focal <- colnames(commun_mat_for_mpd)[i]


    commun_mat_for_mpd_prunned_to_focal <- commun_mat_for_mpd[-which(commun_mat_for_mpd[,i] == 0),]


    #### NRIfocal calculation
    #Standardized effect size of mpd

    vector_avg_mpd <- c(vector_avg_mpd,
                        mean(focal_nri(commun_mat_for_mpd_prunned_to_focal, phy_dist,species_focal, null.model = "richness",
                                       abundance.weighted = TRUE, runs = 99, iterations = 100)$mpd.obs.z,na.rm=T) )



  }
  #vector_avg_mpd[which(is.na(vector_avg_mpd))] <- 0
  return(vector_avg_mpd)
}


make_community_matrix <- function(x_max,y_max,list_species_from_cpp){

  #x_max2 <- x_max - 2 # because the map's border
  #y_max2 <- y_max - 2 # because the map's border
  x_max2 <- x_max
  y_max2 <- y_max
  x_coordinates <- NULL
  y_coordinates <- NULL

  for (i in 1:x_max2){
    for (ii in 1:y_max2){
      x_coordinates <- c(x_coordinates, i)
      y_coordinates <- c(y_coordinates,ii)

    }
  }
  all_coordinates <- as.data.frame(cbind(X=x_coordinates,Y=y_coordinates))
  vector_colnames <- c("X","Y")
  abunda_total <- NULL
  all_pres <- NULL
  alleles_pooled <- matrix(c(NA,0),ncol=2,nrow=nrow(all_coordinates),byrow=T)
  for(iijj in 1:length(list_species_from_cpp)){
    if(list_species_from_cpp[[iijj]]$Death == 0){
      vector_colnames <- c(vector_colnames,paste0("t",list_species_from_cpp[[iijj]]$ID))
      pres_abs <- rep(0,nrow(all_coordinates))
      do_this_table <- list_species_from_cpp[[iijj]]$per_species_population
      for(i in 1:nrow(do_this_table)){
        found_community <- which(do_this_table[i,1] == all_coordinates$X & do_this_table[i,2] == all_coordinates$Y)
        pres_abs[found_community] <- do_this_table[i,3]

        alleles_pooled[found_community,1] <- mean(c(alleles_pooled[found_community,1],
                                                    vegan::diversity(c(as.numeric(do_this_table[i,29:ncol(do_this_table)])),index="shannon")),
                                                  na.rm = T)


        alleles_pooled[found_community,2] <- sum(alleles_pooled[found_community,2],
                                                 length(which(as.numeric(do_this_table[i,29:ncol(do_this_table)])!=0)))






        abunda_total <- c(abunda_total, as.numeric(do_this_table[i,3]))
      }
      all_pres <- cbind(all_pres,pres_abs)

    }
  }

  total_occupancy <- length(which(apply(all_pres,1,FUN=sum,na.rm=TRUE) == 0))

  total_abundance_size <- as.numeric(apply(all_pres,2,FUN=sum,na.rm=TRUE)) #

  # sum(abunda_total)
  all_coordinates <- cbind(all_coordinates,all_pres)
  colnames(all_coordinates) <- vector_colnames
  return(list(community_matrix = all_coordinates,
              total_occupancy = total_occupancy,
              total_abundance_size = total_abundance_size,
              alleles_pooled = alleles_pooled))
}


make_table_gda <- function(do_centroid_pop,do_this_table,sampling_runs,sampling_divided_by,this_index = "shannon"){
  #sampling_runs is not longer use as I set up the spatial centroid of the geographic distribution
  table_gda <- NULL
  geo_dista <- as.matrix(dist(cbind(do_this_table[,1], do_this_table[,2])))
  sampled_table_and_centroid <- prepare_sampled_table_and_find_centroid(do_this_table,num_sampled_cells=1,only_edges = T,sample_it = T, strict_geographic = T)



  centroid_pop <- which(do_this_table$X == sampled_table_and_centroid$centroid_population_distance_abundance$X &  do_this_table$Y == sampled_table_and_centroid$centroid_population_distance_abundance$Y)
  if(do_centroid_pop){
    initial_population <- centroid_pop
  } else {
    initial_population <- sample(1:nrow(do_this_table),1)
  }

dist_to_initial <- geo_dista[initial_population,]
  #dist_to_initial <- dist_to_initial[-which(dist_to_initial == 0)]
  dist_to_initial_sorted <- sort(unique(round(dist_to_initial)))

  table_gda_this_sampling <- NULL
  accumu_to_sampled_populations <- initial_population
  for(i in 2:length(dist_to_initial_sorted)){
    popul_within_radius <- which(dist_to_initial_sorted[i] >= dist_to_initial)
    if(length(popul_within_radius) == 1){
      popul_within_radius <-  which(dist_to_initial_sorted[i + 1] >= dist_to_initial)
    }


    amount_to_sample <- ceiling(length(popul_within_radius)/sampling_divided_by)
    if(amount_to_sample == 1){
      amount_to_sample <- 2
    }
    #to_sampled_populations <- as.numeric(names(sample(popul_within_radius,amount_to_sample,replace=FALSE)))
    if(length(popul_within_radius) == 1){ # to avoid undesiarable behaviour of sample function when pool length == 1. Though because i is never 1, pool is never < 2
      to_sampled_populations <- as.numeric(popul_within_radius)
    } else {
      to_sampled_populations <- as.numeric((sample(popul_within_radius,amount_to_sample,replace=FALSE)))
    }
    to_sampled_populations
    accumu_to_sampled_populations <- c(accumu_to_sampled_populations,to_sampled_populations)
    #cat("accumu_to_sampled_populations: ", accumu_to_sampled_populations,"\n")
    sampled_populations <- as.matrix(do_this_table[unique(accumu_to_sampled_populations),],ncol=53)
    # print(sampled_populations)
    colnames(sampled_populations) <- colnames(do_this_table)
    pooled_populations <- sampled_populations
    if(nrow(pooled_populations) > 1){
      pooled_populations <- t(data.frame(apply(pooled_populations,2,FUN=sum,na.rm=TRUE)))
    }

    pooled_populations <- matrix(pooled_populations,nrow = 1, ncol = 53)
    colnames(pooled_populations) <- colnames(do_this_table)

    the_radius <- 3.1416 * (dist_to_initial_sorted[i] ^ 2)
    #the_radius <- nrow(sampled_populations)

    #sampled_populations <- rbind(sampled_populations,sampled_populations)

       hud_fst_neutral <- mean(compute_betaPopulationsSpecificFST (sampled_populations,pairwise = FALSE, locus="neutral")$computed_betas$betaiovl,na.rm=T)
    hud_fst_selection <- mean(compute_betaPopulationsSpecificFST (sampled_populations,pairwise = FALSE, locus="selection")$computed_betas$betaiovl,na.rm=T)
    hud_fst_together <- mean(compute_betaPopulationsSpecificFST (sampled_populations,pairwise = FALSE, locus="together")$computed_betas$betaiovl,na.rm=T)


    ind_here <- sum(sampled_populations[,3])

    alleles_selection <-length(which(as.numeric(pooled_populations[,4:28]) != 0))
    alleles_neutral <-length(which(as.numeric(pooled_populations[,29:ncol(pooled_populations)]) != 0))
    alleles_bothLoci <-    length(which(as.numeric(pooled_populations[,4:ncol(pooled_populations)]) != 0))



    the_hs <- as.data.frame(compute_H_values(pooled_populations,"shannon"))
    table_gda_this_sampling <- as.data.frame(rbind(table_gda_this_sampling,
                                                   cbind(radius = the_radius,
                                                         fst_bothLoci = hud_fst_together,
                                                         fst_selection = hud_fst_selection,
                                                         fst_value_neutral = hud_fst_neutral,
                                                         H_value_bothLoci =   as.numeric(the_hs$H_value_bothLoci),
                                                         H_value_selection = as.numeric(the_hs$H_value_selection),
                                                         H_value_neutral = as.numeric(the_hs$H_value_neutral),
                                                         ind_here = ind_here,
                                                         alleles_selection = alleles_selection,
                                                         alleles_neutral = alleles_neutral,
                                                         alleles_bothLoci = alleles_bothLoci)))






  }
  table_gda <- rbind(table_gda, table_gda_this_sampling)



  return(list(table_gda = table_gda,
              sampled_table_and_centroid = sampled_table_and_centroid))
}

# make_table_gda <- function(do_this_table,sampling_runs,sampling_divided_by,this_index = "shannon"){
#   #sampling_runs is not longer use as I set up the spatial centroid of the geographic distribution
#   table_gda <- NULL
#   geo_dista <- as.matrix(dist(cbind(do_this_table[,1], do_this_table[,2])))
#   sampled_table_and_centroid <- prepare_sampled_table_and_find_centroid(do_this_table,num_sampled_cells=1,only_edges = T,sample_it = T, strict_geographic = T)
#
#   for(io in 1:sampling_runs){
#     if(sampling_runs == 1){
#       centroid_pop <- which(do_this_table$X == sampled_table_and_centroid$centroid_population_distance_abundance$X &  do_this_table$Y == sampled_table_and_centroid$centroid_population_distance_abundance$Y)
#
#     } else {
#       centroid_pop <- sample(1:nrow(do_this_table),1)
#     }
#
#
#
#     initial_population <- centroid_pop
#
#     dist_to_initial <- geo_dista[initial_population,]
#     #dist_to_initial <- dist_to_initial[-which(dist_to_initial == 0)]
#     dist_to_initial_sorted <- sort(unique(round(dist_to_initial)))
#
#     table_gda_this_sampling <- NULL
#     accumu_to_sampled_populations <- initial_population
#     for(i in 2:length(dist_to_initial_sorted)){
#       popul_within_radius <- which(dist_to_initial_sorted[i] >= dist_to_initial)
#       if(length(popul_within_radius) == 1){
#         popul_within_radius <-  which(dist_to_initial_sorted[i + 1] >= dist_to_initial)
#       }
#
#
#       amount_to_sample <- ceiling(length(popul_within_radius)/sampling_divided_by)
#       if(amount_to_sample == 1){
#         amount_to_sample <- 2
#       }
#       #to_sampled_populations <- as.numeric(names(sample(popul_within_radius,amount_to_sample,replace=FALSE)))
#       if(length(popul_within_radius) == 1){ # to avoid undesiarable behaviour of sample function when pool length == 1. Though because i is never 1, pool is never < 2
#         to_sampled_populations <- as.numeric(popul_within_radius)
#       } else {
#         to_sampled_populations <- as.numeric((sample(popul_within_radius,amount_to_sample,replace=FALSE)))
#       }
#       to_sampled_populations
#       accumu_to_sampled_populations <- c(accumu_to_sampled_populations,to_sampled_populations)
#       #cat("accumu_to_sampled_populations: ", accumu_to_sampled_populations,"\n")
#       sampled_populations <- as.matrix(do_this_table[unique(accumu_to_sampled_populations),],ncol=53)
#       # print(sampled_populations)
#       colnames(sampled_populations) <- colnames(do_this_table)
#       pooled_populations <- sampled_populations
#       if(nrow(pooled_populations) > 1){
#         pooled_populations <- t(data.frame(apply(pooled_populations,2,FUN=sum,na.rm=TRUE)))
#       }
#
#       pooled_populations <- matrix(pooled_populations,nrow = 1, ncol = 53)
#       colnames(pooled_populations) <- colnames(do_this_table)
#
#       the_radius <- 3.1416 * (dist_to_initial_sorted[i] ^ 2)
#       #the_radius <- nrow(sampled_populations)
#
#       #sampled_populations <- rbind(sampled_populations,sampled_populations)
#       hud_fst_neutral <- mean(compute_betaPopulationsSpecificFST (sampled_populations,pairwise = FALSE, locus="neutral")$betaiovl,na.rm=T)
#       hud_fst_selection <- mean(compute_betaPopulationsSpecificFST (sampled_populations,pairwise = FALSE, locus="selection")$betaiovl,na.rm=T)
#       hud_fst_together <- mean(compute_betaPopulationsSpecificFST (sampled_populations,pairwise = FALSE, locus="together")$betaiovl,na.rm=T)
#       ind_here <- sum(sampled_populations[,3])
#
#       alleles_selection <-length(which(as.numeric(pooled_populations[,4:28]) != 0))
#       alleles_neutral <-length(which(as.numeric(pooled_populations[,29:ncol(pooled_populations)]) != 0))
#       alleles_bothLoci <-    length(which(as.numeric(pooled_populations[,4:ncol(pooled_populations)]) != 0))
#
#
#
#       the_hs <- as.data.frame(compute_H_values(pooled_populations,"shannon"))
#       table_gda_this_sampling <- as.data.frame(rbind(table_gda_this_sampling,
#                                                      cbind(radius = the_radius,
#                                                            fst_bothLoci = hud_fst_together,
#                                                            fst_selection = hud_fst_selection,
#                                                            fst_value_neutral = hud_fst_neutral,
#                                                            H_value_bothLoci =   as.numeric(the_hs$H_value_bothLoci),
#                                                            H_value_selection = as.numeric(the_hs$H_value_selection),
#                                                            H_value_neutral = as.numeric(the_hs$H_value_neutral),
#                                                            ind_here = ind_here,
#                                                            alleles_selection = alleles_selection,
#                                                            alleles_neutral = alleles_neutral,
#                                                            alleles_bothLoci = alleles_bothLoci)))
#
#
#
#
#
#
#     }
#     table_gda <- rbind(table_gda, table_gda_this_sampling)
#   }
#
#
#   return(list(table_gda = table_gda,
#               sampled_table_and_centroid = sampled_table_and_centroid))
# }



# make_table_gda_old <- function(do_this_table,sampling_runs,sampling_divided_by){
#   table_gda <- NULL
#       geo_dista <- as.matrix(dist(cbind(do_this_table[,1], do_this_table[,2])))
#   for(jjj in 1:sampling_runs){
#     initial_population <- sample(1:nrow(do_this_table),1)
#
#     dist_to_initial <- round(geo_dista[initial_population,])
#     dist_to_initial <- dist_to_initial[-which(dist_to_initial == 0)]
#     radiuses <- sort(unique(dist_to_initial))
#     table_gda_this_sampling <- NULL
#     for(i in 1:length(radiuses)){
#       popul_within_radius <- which(radiuses[i] >= dist_to_initial)
#       amount_to_sample <- ceiling(length(popul_within_radius)/sampling_divided_by)
#       to_sampled_populations <- sample(popul_within_radius,amount_to_sample,replace=FALSE)
#
#
#       pooled_populations <- do_this_table[to_sampled_populations,]
#       if(amount_to_sample > 1){
#         pooled_populations <- t(data.frame(apply(pooled_populations,2,FUN=sum,na.rm=TRUE)))
#       }
#
#       pooled_populations <- matrix(pooled_populations,nrow = 1, ncol = 13)
#       colnames(pooled_populations) <- colnames(do_this_table)
#
#       table_gda_this_sampling <- as.data.frame(rbind(table_gda_this_sampling,
#                                                      cbind(radius = radiuses[i],compute_H_values(pooled_populations))))
#
#     }
#     table_gda <- rbind(table_gda, table_gda_this_sampling)
#
#
#   }
#   return(table_gda)
# }


condense_spatial_genDiver <- function(neutral_xy,summarize_by = "mean"){
  if(summarize_by != "mean"){
    stop("not implemente for something else than mean")
  }
  growing_space <- NULL

  run_while <- TRUE
  while(run_while){
    subset <- which(neutral_xy[1,1] ==  as.numeric(neutral_xy[,1]) & neutral_xy[1,2] ==  as.numeric(neutral_xy[,2]))


    if(length(subset) > 2){
      growing_space <- rbind(growing_space,
                             apply(neutral_xy[subset,],2,FUN=mean,na.rm=TRUE))
    } else {
      growing_space <- rbind(growing_space,
                             neutral_xy[subset,])
    }



    neutral_xy <- neutral_xy[-subset,]
    if(is.null(nrow(neutral_xy)) || nrow(neutral_xy) == 0){
      run_while <- FALSE
    }

  }
  return(genetic_diversity_space = growing_space)
}


phylosignal_from_speciestable <- function (table_species_rates_time_NEW){

  metrics_to_do <- c("meanbetaFST_neutral","sd_betaFST_neutral","meanbetaFST_selection","sd_betaFST_selection","Range","H_selection",
                     "rate_expansion_N_early","rate_expansion_N_late","rate_contraction_N_early","rate_contraction_N_late",
                     "rate_expansion_S_early","rate_expansion_S_late","rate_contraction_S_early","rate_contraction_S_late",
                     "gradient_C_S_fst","gradient_C_N_fst","gradient_C_S_H","gradient_C_N_H","total_N_change","total_S_change","penalty_distance","events_expansion_N",
                     "events_contraction_N","events_expansion_S","events_contraction_S")

  names_for_from_phylosignal <- NULL
  for(ij in 1:length(metrics_to_do)){
    variable_to_do <- metrics_to_do[ij]
    names_for_from_phylosignal <- c(names_for_from_phylosignal,
                                    paste0("k_",variable_to_do))
    names_for_from_phylosignal <- c(names_for_from_phylosignal,
                                    paste0("pvalue_",variable_to_do))
  }


  p_ksignalvalue <- NULL
  all_species_table <- NULL

  for(iijjkk in 1:2){

    global_case <- c("lowK","highK")[iijjkk]

    for (iijj in 1:3){
      lambda_case <- c("highLambd","lowLambd","tinylambd")[iijj]

      for(ii in 1:2){
        param <- c(1,3)[ii]

        for(ijk in 1:length(all_replicates)){
          colnames(table_species_rates_time_NEW)

          species_table <- table_species_rates_time_NEW[which(table_species_rates_time_NEW$global_case == global_case), ]
          species_table <- species_table[which(species_table$lambda_case == lambda_case ),]
          species_table <- species_table[which(species_table$param_comb == param ),]
          species_table <- species_table[which(species_table$replicate == ijk ),]



          if(nrow(species_table) > 3){


            try(phylo_tree_and_sisters <- build_tree_and_sister_pairs(time_simulated_from_output,species_table)
                ,silent=TRUE)

            if(is.null(phylo_tree_and_sisters) == FALSE && length(which(species_table$Range > filtering_rangeSize)) > 3){
              phylotree <- phylo_tree_and_sisters$phyloTree


              from_phylosignal <- NULL
              for(ij in 1:length(metrics_to_do)){
                variable_to_do <- metrics_to_do[ij]


                what_with_NAs <- "zero" # "keep" "prune" "zero"
                what_with_NAs <- "keep" # "keep" "prune" "zero"
                what_with_NAs <- "prune" # "keep" "prune" "zero"

                #species_table <- subset(species_table,rate_expansion_N_early > 0)

                vector_sp_remove_tree <- NULL
                vector_sp_withT <- paste0("t",species_table$id)
                for(iuj in 1:length(phylotree$tip.label)){

                  if(all(phylotree$tip.label[iuj] != vector_sp_withT)){
                    vector_sp_remove_tree <- c(vector_sp_remove_tree,phylotree$tip.label[iuj])
                  }

                }
                phylotree1 <- drop.tip(phylotree,vector_sp_remove_tree)

                sorted_variable_and_tree <- extract_variable_sorted(phylotree1,species_table,variable_to_do,what_with_NAs,filtering_rangeSize = filtering_rangeSize)




                phylotree2 <- sorted_variable_and_tree$phylotree2
                if(length(phylotree2$tip.label) > 2){
                  sorted_variable <- sorted_variable_and_tree$sorted_variable
                  variable_to_plot_inTree <- setNames(sorted_variable,
                                                      phylotree2$tip.label)


                  # fit <- fastAnc(phylotree2,variable_to_plot_inTree,vars=TRUE,CI=TRUE)
                  #
                  # obj <- contMap(phylotree2,variable_to_plot_inTree,plot=FALSE)
                  # plot(obj,type="fan",legend=0.7*max(nodeHeights(phylotree2)),
                  #      fsize=c(0.7,0.9))

                  #phenogram(phylotree2,variable_to_plot_inTree,fsize=0.6,spread.costs=c(1,0))

                  #phylosignal_this_variable <- phylosignal(variable_to_plot_inTree,phylotree2, method="lambda", test=TRUE, nsim=999)


                  phylosignal_this_variable <- phylosig(phylotree2,variable_to_plot_inTree, method="lambda", test=TRUE, nsim=999)
                  from_phylosignal <- c(from_phylosignal,
                                        phylosignal_this_variable$lambda,
                                        phylosignal_this_variable$P)
                  spp_prunned <- length(phylotree2$tip.label)
                } else {

                  from_phylosignal <- c(from_phylosignal,
                                        NA,
                                        NA)
                  spp_prunned <- nrow(species_table)
                }


              }
            } else {
              from_phylosignal <- rep(c(NA,NA),length(metrics_to_do))
              spp_prunned <- nrow(species_table)
            }

            names(from_phylosignal) <- names_for_from_phylosignal


            p_ksignalvalue <- rbind(p_ksignalvalue,c(from_phylosignal,
                                                     spp_total = nrow(species_table),
                                                     spp_prunned = spp_prunned,
                                                     global_case = global_case,
                                                     lambda_case = lambda_case,
                                                     replicate = ijk,
                                                     parameter  = param,
                                                     avg_total_N_change = mean(species_table$total_N_change),
                                                     avg_total_S_change = mean(species_table$total_S_change),
                                                     geneflow_events = mean(species_table$geneflow_events,na.rm = T)))

          }


        }
      }
    }
  }
  total_p_ksignalvalue_all_scenarios <- as.data.frame(p_ksignalvalue)
  return(total_p_ksignalvalue_all_scenarios)
}

find_border_changes_within_timewindow <- function(list_species_from_cpp,time_simulated_from_output,
                                                  early_limit_timewindow,late_limit_timewindow,
                                                  meanbetaFST_neutral,sd_betaFST_neutral,
                                                  meanbetaFST_selection,sd_betaFST_selection){

  output_table <- NULL
  for(i in 1:length(list_species_from_cpp)){
    this_species <- list_species_from_cpp[[i]]
    #### new part

    be_born <- time_simulated_from_output - this_species$Birth
    if(this_species$Death != 0){ # extinct spp
      be_dead <- time_simulated_from_output - this_species$Death
    } else {
      be_dead <-  this_species$Death
    }

    mid_age <- (be_born - be_dead) / 2
    ealy_years <- (be_born - be_dead) / 10 # the first 10% percent of the lifetime of the species
    old_years <- (be_born - be_dead) - ( (be_born - be_dead) / 10) # the last 10% percent of the lifetime of the species

    # north border
    vector_change_northernmost <- this_species$change_northernmost
    vector_time_change_northernmost <- this_species$time_change_northernmost
    vector_time_change_northernmost_ageConsidered <- vector_time_change_northernmost  - (time_simulated_from_output - be_born)

    expansion_N_time <- NULL
    contraction_N_time <- NULL
    if(length(vector_change_northernmost) > 0){
      for(iijk in 1:length(vector_time_change_northernmost_ageConsidered)){
        if(vector_change_northernmost[iijk] == 1){

          expansion_N_time <- c(expansion_N_time,vector_time_change_northernmost_ageConsidered[iijk])
        } else {
          contraction_N_time <-  c(contraction_N_time,vector_time_change_northernmost_ageConsidered[iijk])
        }
      }

    }


    expan_N_within_window <-  length(which(early_limit_timewindow >= expansion_N_time & late_limit_timewindow <= expansion_N_time ))
    contrac_N_within_window <-  length(which(early_limit_timewindow >= contraction_N_time & late_limit_timewindow <= contraction_N_time ))


    # south border
    vector_change_southernmost <- this_species$change_southernmost
    vector_time_change_southernmost <- this_species$time_change_southernmost
    vector_time_change_southernmost_ageConsidered <- vector_time_change_southernmost  - (time_simulated_from_output - be_born)



    expansion_S_time <- NULL
    contraction_S_time <- NULL
    if(length(vector_change_southernmost) > 0){
      for(iijk in 1:length(vector_time_change_southernmost_ageConsidered)){
        if(vector_change_southernmost[iijk] == 1){

          expansion_S_time <- c(expansion_S_time,vector_time_change_southernmost_ageConsidered[iijk])
        } else {
          e <-  c(contraction_S_time,vector_time_change_southernmost_ageConsidered[iijk])
        }
      }

    }
    expan_S_within_window <-  length(which(early_limit_timewindow >= expansion_S_time & late_limit_timewindow <= expansion_S_time ))
    contrac_S_within_window <-  length(which(early_limit_timewindow >= contraction_S_time & late_limit_timewindow <= contraction_S_time ))



    output_table <- rbind(output_table,c(this_species$ID,expan_N_within_window,contrac_N_within_window,
                                         expan_S_within_window,contrac_S_within_window))
  }
  colnames(output_table) <- c("species","expan_N_within_window","contrac_N_within_window",
                              "expan_S_within_window","contrac_S_within_window")
  output_table <- cbind(output_table,meanbetaFST_neutral,sd_betaFST_neutral,meanbetaFST_selection,sd_betaFST_selection)
  return(output_table)
}

optimal_route_distance_based <- function(table_for_penality_distance,the_origin_cell){

  origin_cell <- the_origin_cell
  matrix_distance_penalty <- distance(cbind(table_for_penality_distance$X,table_for_penality_distance$Y), method = "euclidean", use.row.names = TRUE)
  colnames(matrix_distance_penalty) <- rownames(table_for_penality_distance)
  rownames(matrix_distance_penalty) <- rownames(table_for_penality_distance)


  origin_cell1 <- origin_cell
  matrix_distance_penalty_to_work <- matrix_distance_penalty
  total_sum <- 0
  origins_all <-NULL
  #shortest path
  for(i in 1:(nrow(matrix_distance_penalty) - 2)){
    origins_all <- c(origins_all,origin_cell)
    from_cell <- matrix_distance_penalty_to_work[which(rownames(matrix_distance_penalty_to_work) == origin_cell),]
    total_sum <- total_sum + as.numeric(sort(from_cell)[2])
    matrix_distance_penalty_to_work <-  matrix_distance_penalty_to_work[-which(rownames(matrix_distance_penalty_to_work) == origin_cell),]
    matrix_distance_penalty_to_work <-  matrix_distance_penalty_to_work[,-which(colnames(matrix_distance_penalty_to_work) == origin_cell)]
    origin_cell <- names(sort(from_cell)[2]) # new origen cell
  }

  total_sum <- total_sum + matrix_distance_penalty_to_work[1,2]
  #longest
  matrix_distance_penalty_to_work <- matrix_distance_penalty
  total_sum2 <- 0
  for(i in 1:(nrow(matrix_distance_penalty) - 2)){

    from_cell <- matrix_distance_penalty_to_work[which(rownames(matrix_distance_penalty_to_work) == origin_cell1),]
    total_sum2 <- total_sum2 + as.numeric(sort(from_cell)[length(from_cell)]) # to pick the last element, the largest
    matrix_distance_penalty_to_work <-  matrix_distance_penalty_to_work[-which(rownames(matrix_distance_penalty_to_work) == origin_cell1),]
    matrix_distance_penalty_to_work <-  matrix_distance_penalty_to_work[,-which(colnames(matrix_distance_penalty_to_work) == origin_cell1)]
    origin_cell1 <- names(sort(from_cell)[length(from_cell)]) # new origin cell
  }

  total_sum2 <- total_sum2 + matrix_distance_penalty_to_work[1,2]

  #latitudinal path
  latitudinal_route <- sum(matrix_distance_penalty[,which(rownames(matrix_distance_penalty) == the_origin_cell)])


  return(list(longest_route = total_sum2,
              optimal = total_sum,
              latitudinal_route = latitudinal_route))
}


order_colonization_pop_route <- function(sampled_table_to_do,nj_res){
  table_for_evol_histo <- sampled_table_to_do

  ancest <- table_for_evol_histo[which(min(table_for_evol_histo$betaFST)==table_for_evol_histo$betaFST),]
  num_equal_ancest <- nrow(ancest)
  if(num_equal_ancest > 1){

    ancest <- table_for_evol_histo[which(min(table_for_evol_histo$betaFST)==table_for_evol_histo$betaFST),][1,]
    table_for_evol_histo[which(min(table_for_evol_histo$betaFST)==table_for_evol_histo$betaFST),][2:num_equal_ancest,]$betaFST <- ancest$betaFST + 0.0000001

  }
  last_pop_colonized <- rownames(ancest) # it is actually the  ancestral population
  table_for_evol_histo <- table_for_evol_histo[-which(rownames(table_for_evol_histo)==last_pop_colonized),]
  cophenetic_dist_matr <- cophenetic(nj_res)
  cophenetic_dist_matr[which(cophenetic_dist_matr < 0)] <- 0


  id_pops <- rownames(sampled_table_to_do)
  colonization_order <- NULL
  colonization_order <- c(colonization_order,last_pop_colonized)
  id_pops <- id_pops[-which(id_pops == last_pop_colonized )]

  for(ii in 1:(nrow(sampled_table_to_do) - 2)){
    # for(ii in 1:3){
    next_based_cophenetic <- cophenetic_dist_matr[which(rownames(cophenetic_dist_matr)==last_pop_colonized),]
    next_based_cophenetic <- names(sort(next_based_cophenetic))
    next_based_cophenetic <- next_based_cophenetic[-which(next_based_cophenetic == last_pop_colonized)]
    next_based_cophenetic <- as.data.frame(cbind(rank1=1:length(next_based_cophenetic),next_based_cophenetic))

    next_based_totalfst <- rownames(table_for_evol_histo[order(table_for_evol_histo$betaFST,decreasing=F),])
    next_based_totalfst <- as.data.frame(cbind(rank2=1:length(next_based_totalfst),next_based_totalfst))




    sort(as.numeric(next_based_cophenetic[,2]))

    sort(as.numeric(next_based_totalfst[,2]))

    total_ranking <- NULL
    for(i in 1:nrow(next_based_cophenetic)){

      take_this_rank <- next_based_totalfst[which(next_based_cophenetic$next_based_cophenetic[i] == next_based_totalfst$next_based_totalfst),]$rank2
      #cat("doing: " ,i,"\n")
      total_ranking <- c(total_ranking,
                         as.numeric(next_based_cophenetic[i,]$rank1) + as.numeric(take_this_rank))
      #cat(as.numeric(next_based_cophenetic[i,]$rank1) + as.numeric(take_this_rank),"\n")
    }
    # cat(length(next_based_cophenetic$next_based_cophenetic), "\n")
    # cat(length(total_ranking),"\n")
    next_total <- as.data.frame(cbind(next_based_cophenetic$next_based_cophenetic,total_ranking=as.numeric(total_ranking)))

    next_total <- next_total[order(as.numeric(next_total$total_ranking),decreasing=F),]

    #print(cophenetic_dist_matr)
    cophenetic_dist_matr <- cophenetic_dist_matr[-which(rownames(cophenetic_dist_matr)==last_pop_colonized),]
    cophenetic_dist_matr <- cophenetic_dist_matr[,-which(colnames(cophenetic_dist_matr)==last_pop_colonized)]
    last_pop_colonized <- next_total[1,1]
    table_for_evol_histo <- table_for_evol_histo[-which(rownames(table_for_evol_histo)==last_pop_colonized),]

    colonization_order <- c(colonization_order,last_pop_colonized)
    id_pops <- id_pops[-which(id_pops == last_pop_colonized )]
    # cat(colonization_order, "\n")

  }
  colonization_order <- c(colonization_order,
                          id_pops)

  table_for_penality_distance <- sampled_table_to_do


  matrix_distance_penalty <- distance(cbind(table_for_penality_distance$X,table_for_penality_distance$Y),
                                      method = "euclidean", use.row.names = TRUE)
  colnames(matrix_distance_penalty) <- rownames(table_for_penality_distance)
  rownames(matrix_distance_penalty) <- rownames(table_for_penality_distance)

  optimal_route_distance_object <- optimal_route_distance_based(table_for_penality_distance,colonization_order[1])
  optimal_route_distance <- optimal_route_distance_object$optimal
  longest_route_distance <- optimal_route_distance_object$longest_route



  total_sum <- 0
  matrix_distance_penalty_to_work <- matrix_distance_penalty
  for(i in 1:(length(colonization_order) - 1)){
    from_cell <- matrix_distance_penalty_to_work[which(rownames(matrix_distance_penalty_to_work) == colonization_order[i]),]
    to_cell <- colonization_order[i + 1]
    total_sum <- total_sum + as.numeric(from_cell[which(names(from_cell) == to_cell )])

  }


  # if(longest_route_distance < optimal_route_distance){
  #   stop("something odd with optimal_route_distance")
  # }
  #
  # if(total_sum < optimal_route_distance || total_sum > longest_route_distance){
  #   stop("something odd with optimal_route_distance in here")
  # }

  #penalty_distance <- abs(optimal_route_distance - total_sum)/(longest_route_distance -  optimal_route_distance)



  penalty_distance <- (abs(total_sum - optimal_route_distance_object$latitudinal_route))/sum(matrix_distance_penalty)
  penalty_distance <- (abs(total_sum - optimal_route_distance))/sum(matrix_distance_penalty)
  #penalty_distance <- (abs(total_sum - optimal_route_distance))
  return(list(the_order = colonization_order,
              penalty_distance = penalty_distance))

}


steeper_gradient <- function(gradient_factor,map_temperature_vector){

  unique_bands <- unique(map_temperature_vector)
  unique_bands <- unique_bands[-which(unique_bands == -9)]

  intermediate_band <- round(median(unique_bands))

  for(ik in 1:length(map_temperature_vector)){
    if(map_temperature_vector[ik] != -9){

      if(map_temperature_vector[ik] > intermediate_band)

        map_temperature_vector[ik] <- map_temperature_vector[ik] + gradient_factor

      if(map_temperature_vector[ik] < intermediate_band)

        map_temperature_vector[ik] <- map_temperature_vector[ik] - gradient_factor
    }

  }
  return(map_temperature_vector)
}



compute_betaPopulationsSpecificFST <- function(table_to_plot,pairwise, locus="together"){ # locus="both" or "neutral" or "selection"
  locus_selection <- NULL
  locus_neutral <- NULL
  population_id <- NULL
  table_to_plot <- as.data.frame(table_to_plot)
  for(i in 1:nrow(table_to_plot)){
    thispopulation <- table_to_plot[i,]
    locus_selection <- c(locus_selection,rep(1,thispopulation$Allele_A1))
    locus_selection <- c(locus_selection,rep(2,thispopulation$Allele_B1))
    locus_selection <- c(locus_selection,rep(3,thispopulation$Allele_C1))
    locus_selection <- c(locus_selection,rep(4,thispopulation$Allele_D1))
    locus_selection <- c(locus_selection,rep(5,thispopulation$Allele_E1))
    locus_selection <- c(locus_selection,rep(6,thispopulation$Allele_F1))
    locus_selection <- c(locus_selection,rep(7,thispopulation$Allele_G1))
    locus_selection <- c(locus_selection,rep(8,thispopulation$Allele_H1))
    locus_selection <- c(locus_selection,rep(9,thispopulation$Allele_I1))
    locus_selection <- c(locus_selection,rep(10,thispopulation$Allele_J1))
    locus_selection <- c(locus_selection,rep(11,thispopulation$Allele_K1))
    locus_selection <- c(locus_selection,rep(12,thispopulation$Allele_L1))
    locus_selection <- c(locus_selection,rep(13,thispopulation$Allele_M1))
    locus_selection <- c(locus_selection,rep(14,thispopulation$Allele_N1))
    locus_selection <- c(locus_selection,rep(15,thispopulation$Allele_O1))

    locus_selection <- c(locus_selection,rep(16,thispopulation$Allele_P1))
    locus_selection <- c(locus_selection,rep(17,thispopulation$Allele_Q1))
    locus_selection <- c(locus_selection,rep(18,thispopulation$Allele_R1))
    locus_selection <- c(locus_selection,rep(19,thispopulation$Allele_S1))
    locus_selection <- c(locus_selection,rep(20,thispopulation$Allele_T1))
    locus_selection <- c(locus_selection,rep(21,thispopulation$Allele_U1))
    locus_selection <- c(locus_selection,rep(22,thispopulation$Allele_V1))
    locus_selection <- c(locus_selection,rep(23,thispopulation$Allele_W1))
    locus_selection <- c(locus_selection,rep(24,thispopulation$Allele_X1))
    locus_selection <- c(locus_selection,rep(25,thispopulation$Allele_Y1))

    locus_neutral <- c(locus_neutral,rep(1,thispopulation$Allele_A2))
    locus_neutral <- c(locus_neutral,rep(2,thispopulation$Allele_B2))
    locus_neutral <- c(locus_neutral,rep(3,thispopulation$Allele_C2))
    locus_neutral <- c(locus_neutral,rep(4,thispopulation$Allele_D2))
    locus_neutral <- c(locus_neutral,rep(5,thispopulation$Allele_E2))
    locus_neutral <- c(locus_neutral,rep(6,thispopulation$Allele_F2))
    locus_neutral <- c(locus_neutral,rep(7,thispopulation$Allele_G2))
    locus_neutral <- c(locus_neutral,rep(8,thispopulation$Allele_H2))
    locus_neutral <- c(locus_neutral,rep(9,thispopulation$Allele_I2))
    locus_neutral <- c(locus_neutral,rep(10,thispopulation$Allele_J2))
    locus_neutral <- c(locus_neutral,rep(11,thispopulation$Allele_K2))
    locus_neutral <- c(locus_neutral,rep(12,thispopulation$Allele_L2))
    locus_neutral <- c(locus_neutral,rep(13,thispopulation$Allele_M2))
    locus_neutral <- c(locus_neutral,rep(14,thispopulation$Allele_N2))
    locus_neutral <- c(locus_neutral,rep(15,thispopulation$Allele_O2))

    locus_neutral <- c(locus_neutral,rep(16,thispopulation$Allele_P2))
    locus_neutral <- c(locus_neutral,rep(17,thispopulation$Allele_Q2))
    locus_neutral <- c(locus_neutral,rep(18,thispopulation$Allele_R2))
    locus_neutral <- c(locus_neutral,rep(19,thispopulation$Allele_S2))
    locus_neutral <- c(locus_neutral,rep(20,thispopulation$Allele_T2))
    locus_neutral <- c(locus_neutral,rep(21,thispopulation$Allele_U2))
    locus_neutral <- c(locus_neutral,rep(22,thispopulation$Allele_V2))
    locus_neutral <- c(locus_neutral,rep(23,thispopulation$Allele_W2))
    locus_neutral <- c(locus_neutral,rep(24,thispopulation$Allele_X2))
    locus_neutral <- c(locus_neutral,rep(25,thispopulation$Allele_Y2))



    population_id <- c(population_id,rep(i,thispopulation$Pop_size))
  }

  if(locus == "together"){
    table_for_bfts <- as.data.frame(cbind(population_id,locus_selection,locus_neutral))
  }

  if(locus == "neutral"){
    table_for_bfts <- as.data.frame(cbind(population_id,locus_neutral,locus_neutral))
  }
  if(locus == "selection"){
    table_for_bfts <- as.data.frame(cbind(population_id,locus_selection,locus_selection))
  }

  if(nrow(table_for_bfts) > 1){
    wc_fst <- NA
    try(wc_fst <- wc(table_for_bfts,diploid=F),
        ,silent=TRUE)
    if(pairwise){
      computed_betas <- pairwise.betas(table_for_bfts,diploid=F)
    } else {
      computed_betas <- betas(table_for_bfts,nboot=100,diploid=F)
    }

  } else {
    computed_betas$betaiovl <- NA
    wc_fst <- NA
  }



  return(list(computed_betas = computed_betas,
              wc_fst = wc_fst))
}


do_averaging_table <- function(average_table_to_plot,table_to_plot){
  averageD_table_to_plot <- average_table_to_plot

  for(ii in 1:nrow(table_to_plot)){

    matched_x <- which(table_to_plot[ii,]$X == averageD_table_to_plot$X)
    matched_y <- which(table_to_plot[ii,]$Y == averageD_table_to_plot$Y)

    if(length(unique(c(matched_x,matched_y))) < length(c(matched_x,matched_y)) ){
      cell_found <- intersect(matched_x,matched_y)

      averaged <- apply(rbind(averageD_table_to_plot[cell_found,],table_to_plot[ii,])
                        ,2,sum)
      averageD_table_to_plot[cell_found,-c(1,2)] <- averaged[-c(1,2)]
    } else{
      cat("this xy:",table_to_plot[ii,]$X,"-",table_to_plot[ii,]$Y, "\n")
      averageD_table_to_plot <- rbind(table_to_plot[ii,],averageD_table_to_plot)
    }
  }

  return(averageD_table_to_plot)
}

process_rangesize_elevation_bands <- function(list_species_from_cpp,band_limits,range_size){

  table_all_species <- NULL
  for( i in 1:length(list_species_from_cpp)){
    workonthis <- list_species_from_cpp[[i]]
    if(workonthis$Death == 0){
      range_in_bands <- 0
      for(ii in 1:length(band_limits)){
        range_band_from_tipdown <- length(which(workonthis$distribution[,1] <=  band_limits[ii]))
        # I now substract the cells the species was present in previous bands
        range_in_bands <- c(range_in_bands,range_band_from_tipdown-sum(range_in_bands))

      }
      range_in_bands <- c(range_in_bands,workonthis$RangeSize - sum(range_in_bands)) # to do the lowest band
      range_in_bands <- range_in_bands[-1]

      if(any(c(workonthis$RangeHighlands,workonthis$RangeIntermediate1,
               workonthis$RangeIntermediate2,workonthis$RangeLowlands) != range_in_bands)){
        stop("elevational range from cpp does not match from distribution info")
      }




      info_one_species <- c(workonthis$ID,workonthis$Birth,range_in_bands)
      table_all_species <- rbind(table_all_species,info_one_species)
    }

  }
  rownames(table_all_species)<-NULL
  if(sum(table_all_species[,-c(1,2)]) != sum(range_size)){
    stop("something wrong with making the elevational bands richnesses")
  }
  return(data.frame(table_all_species))
}



process_traits_ages_idspp_elevation_bands <- function(list_species_from_cpp,map_elevation){
  table_presence_elevation <- NULL

  table_mean_temperature_elevation <- NULL
  table_sd_temperature_elevation <- NULL
  table_mean_temperature_elevation <- NULL
  table_sd_temperature_elevation <- NULL
  species_id <- NULL
  species_age <- NULL
  species_elevational_origin <- NULL
  for (i in 1:length(list_species_from_cpp)){
    this_species <- list_species_from_cpp[[i]]
    this_species_highlands_trait <- NULL
    this_species_intermediate2_trait <- NULL
    this_species_intermediate1_trait <- NULL
    this_species_lowlands_trait <- NULL

    this_species_highlands <- 0
    this_species_intermediate2 <- 0
    this_species_intermediate1 <- 0
    this_species_lowlands<- 0

    if(this_species$Death == 0){

      species_id <- c(species_id,i)
      species_age <- c(species_age,this_species$Birth)
      species_elevational_origin <- c(species_elevational_origin,this_species$elevational_origin)
      for(ii in 1:nrow(this_species$distribution)){
        elevation_from_map <- map_elevation[this_species$distribution[ii,][1],this_species$distribution[ii,][2]]
        if(elevation_from_map == 1){
          this_species_lowlands_trait <- c(this_species_lowlands_trait,this_species$temperature_evolved[ii])
          this_species_lowlands <- this_species_lowlands + 1
        }
        if(elevation_from_map == 2){
          this_species_intermediate2_trait <- c(this_species_intermediate2_trait,this_species$temperature_evolved[ii])
          this_species_intermediate2 <- this_species_intermediate2 + 1
        }
        if(elevation_from_map == 3){
          this_species_intermediate1_trait <- c(this_species_intermediate1_trait,this_species$temperature_evolved[ii])
          this_species_intermediate1 <- this_species_intermediate1 + 1
        }
        if(elevation_from_map == 4){
          this_species_highlands_trait <- c(this_species_highlands_trait,this_species$temperature_evolved[ii])
          this_species_highlands <- this_species_highlands + 1
        }

      }

      if(any(c(this_species$RangeLowlands,this_species$RangeIntermediate2, this_species$RangeIntermediate1,this_species$RangeHighlands)!=
             c(this_species_lowlands,this_species_intermediate2,this_species_intermediate1,this_species_highlands))){
        stop("ranges per elevation from dont match map")
      }



      table_presence_elevation <- rbind(table_presence_elevation,c(this_species_lowlands,this_species_intermediate2,this_species_intermediate1,this_species_highlands))
      table_mean_temperature_elevation <- rbind(table_mean_temperature_elevation,c(
        mean(this_species_lowlands_trait,na.rm = T),
        mean(this_species_intermediate2_trait,na.rm = T),
        mean(this_species_intermediate1_trait,na.rm = T),
        mean(this_species_highlands_trait,na.rm = T))
      )
      table_sd_temperature_elevation <- rbind(table_sd_temperature_elevation,c(
        sd(this_species_lowlands_trait,na.rm = T),
        sd(this_species_intermediate2_trait,na.rm = T),
        sd(this_species_intermediate1_trait,na.rm = T),
        sd(this_species_highlands_trait,na.rm = T))
      )


    }
  }
  column_names <- c("spID","age","eleva_origin","lowlands","intermediate2","intermediate1","highlands")
  table_presence_elevation <- cbind(species_id,species_age,species_elevational_origin,table_presence_elevation)
  table_mean_temperature_elevation <- cbind(species_id,species_age,species_elevational_origin,table_mean_temperature_elevation)
  table_sd_temperature_elevation <- cbind(species_id,species_age,species_elevational_origin,table_sd_temperature_elevation)
  colnames(table_presence_elevation) <- column_names
  colnames(table_mean_temperature_elevation) <- column_names
  colnames(table_sd_temperature_elevation) <- column_names

  return(list(table_presence_elevation=table_presence_elevation,
              table_mean_temperature_elevation=table_mean_temperature_elevation,
              table_sd_temperature_elevation=table_sd_temperature_elevation))
}

extract_variable_sorted <- function(phylotree,species_table,variable_to_do,what_with_NAs = "keep",filtering_rangeSize = NULL){
  sorted_variable <- NULL
  column_variable <- which(variable_to_do == colnames(species_table))

  species_table1 <-  species_table[which(species_table$Range > filtering_rangeSize),]
  tips_to_remove_underRangesize <- species_table[-which(species_table$Range > filtering_rangeSize),]$id
  tips_to_remove_underRangesize <- paste0("t",tips_to_remove_underRangesize)
  phylotree1 <-  drop.tip(phylotree,tips_to_remove_underRangesize)

  for(i in 1:length(phylotree1$tip.label)){

    id_to_find <-  as.numeric(strsplit(phylotree1$tip.label[i],split="t")[[1]][2])


    sorted_variable <- c(sorted_variable,species_table1[which(species_table1$id == id_to_find),column_variable])

  }

  if(what_with_NAs == "keep"){
    sorted_variable <- sorted_variable
    phylotree2 <- phylotree1
  }
  if(what_with_NAs == "zero"){
    sorted_variable[which(is.na(sorted_variable))] <- 0
    phylotree2 <- phylotree1
  }

  if(what_with_NAs == "prune"){

    spp_to_remove <- phylotree1$tip.label[which(is.na(sorted_variable))]
    phylotree2 <-  drop.tip(phylotree1,spp_to_remove)
    if(any(is.na(sorted_variable))){
      sorted_variable <- sorted_variable[-which(is.na(sorted_variable))]
    } else {
      sorted_variable <- sorted_variable
    }




  }


  return(list(sorted_variable = sorted_variable,
              phylotree2 = phylotree2))

}



treeANDLtable<- function(L,set_simulated_time){
  age <- set_simulated_time
  # age <- max(c(LTable[,1],LTable[,4]))
  tableConver <- L
  tableConver[tableConver[,4]==0,4] <- -1 # to go with Rampal's notation
  LTable <- matrix(c(tableConver[,2],tableConver[,3],tableConver[,1],tableConver[,4]),ncol=4)
  #age <- max(c(LTable[,1],LTable[,4]))
  LTable[,1] <- age - c(LTable[,1])
  notmin1 <- which(LTable[,4] != -1)
  LTable[notmin1,4] <- age - c(LTable[notmin1,4])
  LTable[which(LTable[,4] == age + 1),4] <- -1
  phyloTree <- L2phylo(LTable[,1:4],dropextinct = T)
  return(list(LTable=LTable,phyloTree=phyloTree))
}



sort_the_range <- function(phylotree,range_size){

  # there is a problem with function sort() when working with species id of the shape "t33"
  id_species_withT <- phylotree$tip.label
  species_id_numbersOnly <- NULL
  for(i in 1:length(phylotree$tip.label)){

    all_species_id <- strsplit(id_species_withT,split="t")
    species_id_numbersOnly <- c(species_id_numbersOnly,as.numeric(all_species_id[[i]][2]))
  }
  species_id_new <- NULL
  for(i in 1:length(phylotree$tip.label)){
    species_id_new <- c(species_id_new, paste0("t",sort(species_id_numbersOnly)[i]))
  }
  range_size2 <- range_size[-which(range_size==0)]
  table_range_size2_sp <- as.data.frame(cbind(species_id_new,range_size2))
  return(sortingtraits(table_range_size2_sp,phylotree))
}

sort_the_traits <- function (phylotree,traits_to_sort){

  # there is a problem with function sort() when working with species id of the shape "t33"
  id_species_withT <- phylotree$tip.label
  species_id_numbersOnly <- NULL
  for(i in 1:length(phylotree$tip.label)){

    all_species_id <- strsplit(id_species_withT,split="t")
    species_id_numbersOnly <- c(species_id_numbersOnly,as.numeric(all_species_id[[i]][2]))
  }
  species_id_new <- NULL
  for(i in 1:length(phylotree$tip.label)){
    species_id_new <- c(species_id_new, paste0("t",sort(species_id_numbersOnly)[i]))
  }

  table_traits_sp <- as.data.frame(cbind(species_id_new,traits_to_sort))

  sortingtraits(table_traits_sp,phylotree)
  return(sortingtraits(table_traits_sp,phylotree))
}


richness_overTime <- function(Ltable_to_do,simulated_time,num_time_slices){
  alive_species_overtime <- NULL
  time_to_divide <- seq(0,simulated_time,by = simulated_time/num_time_slices )[-1]
  for(i in 1:length(time_to_divide)){
    row_counting <-  which(Ltable_to_do[,2] <= time_to_divide[i])

    LTable <- Ltable_to_do[row_counting,]

    if(class(LTable)=="matrix"){
      # Check whether any death did not happen within the time frame
      if(any(LTable[,4] > time_to_divide[i]) == TRUE){
        LTable[which(LTable[,4] > time_to_divide[i]),4] <- 0 # as was alive withiin time frame
      }
      alive_species_overtime <- c(alive_species_overtime,length(which(LTable[,4] == 0)))

    } else{
      alive_species_overtime <- c(alive_species_overtime, 1)

    }

  }
  return(list(alive_species_overtime = alive_species_overtime,
              timeslices=time_to_divide))
}

#
# compute_gradient2 <- function (matrix_gradient){
#
#   if(ncol(matrix_gradient) != 4){
#     stop("this function is for 4 elevational bands")
#   }
#
#   average_richness_difference_elevation <- NULL
#   total_richness <- rowSums(matrix_gradient)
#   for(i in 1: nrow(matrix_gradient)){
#     average_richness_difference_elevation <- c (average_richness_difference_elevation,mean(c(matrix_gradient[i,4] - matrix_gradient[i,3],
#                                                                                              matrix_gradient[i,3] - matrix_gradient[i,2],
#                                                                                              matrix_gradient[i,2] - matrix_gradient[i,1])))
#
#   }
#
#
#   average_richness_difference_elevation <- average_richness_difference_elevation/total_richness * 100
#   return(list(average_richness_difference_elevation=average_richness_difference_elevation,
#               total_richness = total_richness))
# }

compute_gradient2 <- function (matrix_gradient){
  percentage <- NULL
  peak_at <- NULL

  if(colnames(matrix_gradient)[1]!= "richness_highlands"){
    stop("re-order your matrix")
  }
  if(ncol(matrix_gradient) != 4){
    stop("this function is for 4 elevational bands")
  }

  peak <- max(matrix_gradient[1,])
  elevation_band_peak <- which(peak == matrix_gradient[1,])
  lowest <- min(matrix_gradient[1,])
  total_richness <- sum(matrix_gradient[1,])

  if(any(elevation_band_peak == 4)){ # peak in Lowlands
    peak_at <- "lowlands"
  }

  if(any(elevation_band_peak == 3)){
    peak_at <- "intermediate2"
  }
  if(any(elevation_band_peak == 2)){
    peak_at <- "intermediate1"
  }

  if(any(elevation_band_peak == 1)){
    peak_at <- "highlands"
  }

  average_richness_difference_elevation <- (sum(peak-matrix_gradient[1,])/3)/total_richness


  return(list(average_richness_difference_elevation=average_richness_difference_elevation,
              total_richness = total_richness,
              peak_at = peak_at))
}


draw_map <- function(table_to_plot,x_max,y_max,born_x,born_y,finalmap_this){
  #table_to_plot <- as.data.frame(list_species_from_cpp[[species_to_do]]$per_species_population)
  if(is.null(born_x)){
    born_x <- 0
    born_y <- 0
  }
  if(finalmap_this != "richness" &&  finalmap_this != "pop_size"){

    table_to_plot <- cbind(table_to_plot,compute_H_values(table_to_plot,this_index = "shannon"))
    table_to_plot$Allele_A1 <- table_to_plot$Allele_A1/table_to_plot$Pop_size
    table_to_plot$Allele_B1 <- table_to_plot$Allele_B1/table_to_plot$Pop_size
    table_to_plot$Allele_C1 <- table_to_plot$Allele_C1/table_to_plot$Pop_size
    table_to_plot$Allele_D1 <- table_to_plot$Allele_D1/table_to_plot$Pop_size
    table_to_plot$Allele_E1 <- table_to_plot$Allele_E1/table_to_plot$Pop_size
    table_to_plot$Allele_F1 <- table_to_plot$Allele_F1/table_to_plot$Pop_size
    table_to_plot$Allele_G1 <- table_to_plot$Allele_G1/table_to_plot$Pop_size
    table_to_plot$Allele_H1 <- table_to_plot$Allele_H1/table_to_plot$Pop_size
    table_to_plot$Allele_I1 <- table_to_plot$Allele_I1/table_to_plot$Pop_size
    table_to_plot$Allele_J1 <- table_to_plot$Allele_J1/table_to_plot$Pop_size
    table_to_plot$Allele_K1 <- table_to_plot$Allele_K1/table_to_plot$Pop_size
    table_to_plot$Allele_L1 <- table_to_plot$Allele_L1/table_to_plot$Pop_size
    table_to_plot$Allele_M1 <- table_to_plot$Allele_M1/table_to_plot$Pop_size
    table_to_plot$Allele_N1 <- table_to_plot$Allele_N1/table_to_plot$Pop_size
    table_to_plot$Allele_O1 <- table_to_plot$Allele_O1/table_to_plot$Pop_size

    table_to_plot$Allele_P1 <- table_to_plot$Allele_P1/table_to_plot$Pop_size
    table_to_plot$Allele_Q1 <- table_to_plot$Allele_Q1/table_to_plot$Pop_size
    table_to_plot$Allele_R1 <- table_to_plot$Allele_R1/table_to_plot$Pop_size
    table_to_plot$Allele_S1 <- table_to_plot$Allele_S1/table_to_plot$Pop_size
    table_to_plot$Allele_T1 <- table_to_plot$Allele_T1/table_to_plot$Pop_size
    table_to_plot$Allele_U1 <- table_to_plot$Allele_U1/table_to_plot$Pop_size
    table_to_plot$Allele_V1 <- table_to_plot$Allele_V1/table_to_plot$Pop_size
    table_to_plot$Allele_W1 <- table_to_plot$Allele_W1/table_to_plot$Pop_size
    table_to_plot$Allele_X1 <- table_to_plot$Allele_X1/table_to_plot$Pop_size
    table_to_plot$Allele_Y1 <- table_to_plot$Allele_Y1/table_to_plot$Pop_size


    table_to_plot$Allele_A2 <- table_to_plot$Allele_A2/table_to_plot$Pop_size
    table_to_plot$Allele_B2 <- table_to_plot$Allele_B2/table_to_plot$Pop_size
    table_to_plot$Allele_C2 <- table_to_plot$Allele_C2/table_to_plot$Pop_size
    table_to_plot$Allele_D2 <- table_to_plot$Allele_D2/table_to_plot$Pop_size
    table_to_plot$Allele_E2 <- table_to_plot$Allele_E2/table_to_plot$Pop_size
    table_to_plot$Allele_F2 <- table_to_plot$Allele_F2/table_to_plot$Pop_size
    table_to_plot$Allele_G2 <- table_to_plot$Allele_G2/table_to_plot$Pop_size
    table_to_plot$Allele_H2 <- table_to_plot$Allele_H2/table_to_plot$Pop_size
    table_to_plot$Allele_I2 <- table_to_plot$Allele_I2/table_to_plot$Pop_size
    table_to_plot$Allele_J2 <- table_to_plot$Allele_J2/table_to_plot$Pop_size
    table_to_plot$Allele_K2 <- table_to_plot$Allele_K2/table_to_plot$Pop_size
    table_to_plot$Allele_L2 <- table_to_plot$Allele_L2/table_to_plot$Pop_size
    table_to_plot$Allele_M2 <- table_to_plot$Allele_M2/table_to_plot$Pop_size
    table_to_plot$Allele_N2 <- table_to_plot$Allele_N2/table_to_plot$Pop_size
    table_to_plot$Allele_O2 <- table_to_plot$Allele_O2/table_to_plot$Pop_size

    table_to_plot$Allele_P2 <- table_to_plot$Allele_P2/table_to_plot$Pop_size
    table_to_plot$Allele_Q2 <- table_to_plot$Allele_Q2/table_to_plot$Pop_size
    table_to_plot$Allele_R2 <- table_to_plot$Allele_R2/table_to_plot$Pop_size
    table_to_plot$Allele_S2 <- table_to_plot$Allele_S2/table_to_plot$Pop_size
    table_to_plot$Allele_T2 <- table_to_plot$Allele_T2/table_to_plot$Pop_size
    table_to_plot$Allele_U2 <- table_to_plot$Allele_U2/table_to_plot$Pop_size
    table_to_plot$Allele_V2 <- table_to_plot$Allele_V2/table_to_plot$Pop_size
    table_to_plot$Allele_W2 <- table_to_plot$Allele_W2/table_to_plot$Pop_size
    table_to_plot$Allele_X2 <- table_to_plot$Allele_X2/table_to_plot$Pop_size
    table_to_plot$Allele_Y2 <- table_to_plot$Allele_Y2/table_to_plot$Pop_size




  }


  # to flip the map
  max_y <- max(table_to_plot$Y) + 1
  # max_x <- max(table_to_plot$X) + 1

  max_y <- y_max + 1
  #  max_x <- x_max + 1
  born_y = max_y - born_y
  table_to_plot$Y = max_y - table_to_plot$Y
  #table_to_plot$X = max_x - table_to_plot$X


  if(finalmap_this == "richness"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = richness ))+
      ggtitle("SpeciesRichness")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }

  if(finalmap_this == "beta_FST"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = betaFST ))+
      ggtitle("beta_FST")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }

  if(finalmap_this == "pop_size"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Pop_size ))+
      ggtitle("pop_size")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "H_value_bothLoci"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = H_value_bothLoci  )) +
      ggtitle("H_value_bothLoci")+
      xlim(1,y_max)+
      ylim(1,x_max) +
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "H_value_selection"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = H_value_selection )) +
      ggtitle("H_value_selection")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "H_value_neutral"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = H_value_neutral)) +
      ggtitle("H_value_neutral")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }

  if(finalmap_this == "allele_a1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_A1  )) +
      ggtitle("allele_A1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_b1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_B1  )) +
      ggtitle("allele_B1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_c1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_C1  )) +
      ggtitle("allele_C1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_d1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_D1  )) +
      ggtitle("allele_D1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_e1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_E1  )) +
      ggtitle("allele_E1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_f1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_F1  )) +
      ggtitle("allele_F1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_g1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_G1  )) +
      ggtitle("allele_g1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_h1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_H1  )) +
      ggtitle("allele_H1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_i1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_I1  )) +
      ggtitle("allele_I1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_j1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_J1  )) +
      ggtitle("allele_J1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_k1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_K1  )) +
      ggtitle("allele_k1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_l1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_L1  )) +
      ggtitle("allele_L1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_m1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_M1  )) +
      ggtitle("allele_M1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_n1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_N1  )) +
      ggtitle("allele_N1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }
  if(finalmap_this == "allele_o1"){
    ggp <- ggplot(table_to_plot, aes(X, Y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = Allele_O1  )) +
      ggtitle("allele_O1")+
      xlim(1,y_max)+
      ylim(1,x_max)+
      geom_point(aes(x = born_x, y = born_y), color = "red", size = 5)
  }



  return(ggp)

}


spatial_pattern_H <- function(x_max,y_max,list_species_from_cpp,species_to_do,finalmap_this){
  table_to_plot <- as.data.frame(list_species_from_cpp[[species_to_do]]$per_species_population)
  midpoint_x <- round(x_max/2)
  midpoint_y <- round(y_max/2)

  distance_population <- NULL
  variable_to_take <- NULL

  for(i in 1:nrow(table_to_plot)){


    eucl_distance <- sqrt( ((table_to_plot$X[i] - midpoint_x)^ 2) +
                             ((table_to_plot$Y[i] - midpoint_y)^ 2)
    )

    distance_population <- c(distance_population,eucl_distance)
    if(finalmap_this == "pop_size"){
      variable_to_take <- c(variable_to_take,table_to_plot$Pop_size[i])
    }
    if(finalmap_this == "H_value_bothLoci"){

      variable_to_take <- c(variable_to_take,table_to_plot$H_value_bothLoci[i])
    }
    if(finalmap_this == "H_value_selection"){

      variable_to_take <- c(variable_to_take,table_to_plot$H_value_selection[i])
    }
    if(finalmap_this == "H_value_neutral"){

      variable_to_take <- c(variable_to_take,table_to_plot$H_value_neutral[i])
    }
  }
  spatial_pattern <- data.frame(distance_population,variable_to_take)
  colnames(spatial_pattern) <- c("Eucl_dist",finalmap_this)
  return(spatial_pattern)
}

compute_H_values <- function(table_per_species_population,this_index = "heterozygosity"){
  table_per_species_population <- as.data.frame(table_per_species_population)

  all_H_values <- NULL
  all_H_selection <- NULL
  all_H_neutral <- NULL

  Allele_A1 <- table_per_species_population$Allele_A1
  Allele_B1 <- table_per_species_population$Allele_B1
  Allele_C1 <- table_per_species_population$Allele_C1
  Allele_D1 <- table_per_species_population$Allele_D1
  Allele_E1 <- table_per_species_population$Allele_E1
  Allele_F1 <- table_per_species_population$Allele_F1
  Allele_G1 <- table_per_species_population$Allele_G1
  Allele_H1 <- table_per_species_population$Allele_H1
  Allele_I1 <- table_per_species_population$Allele_I1
  Allele_J1 <- table_per_species_population$Allele_J1
  Allele_K1 <- table_per_species_population$Allele_K1
  Allele_L1 <- table_per_species_population$Allele_L1
  Allele_M1 <- table_per_species_population$Allele_M1
  Allele_N1 <- table_per_species_population$Allele_N1
  Allele_O1 <- table_per_species_population$Allele_O1
  Allele_P1 <- table_per_species_population$Allele_P1
  Allele_Q1 <- table_per_species_population$Allele_Q1
  Allele_R1 <- table_per_species_population$Allele_R1
  Allele_S1 <- table_per_species_population$Allele_S1
  Allele_T1 <- table_per_species_population$Allele_T1
  Allele_U1 <- table_per_species_population$Allele_U1
  Allele_V1 <- table_per_species_population$Allele_V1
  Allele_W1 <- table_per_species_population$Allele_W1
  Allele_X1 <- table_per_species_population$Allele_X1
  Allele_Y1 <- table_per_species_population$Allele_Y1


  Allele_A2 <- table_per_species_population$Allele_A2
  Allele_B2 <- table_per_species_population$Allele_B2
  Allele_C2 <- table_per_species_population$Allele_C2
  Allele_D2 <- table_per_species_population$Allele_D2
  Allele_E2 <- table_per_species_population$Allele_E2
  Allele_F2 <- table_per_species_population$Allele_F2
  Allele_G2 <- table_per_species_population$Allele_G2
  Allele_H2 <- table_per_species_population$Allele_H2
  Allele_I2 <- table_per_species_population$Allele_I2
  Allele_J2 <- table_per_species_population$Allele_J2
  Allele_K2 <- table_per_species_population$Allele_K2
  Allele_L2 <- table_per_species_population$Allele_L2
  Allele_M2 <- table_per_species_population$Allele_M2
  Allele_N2 <- table_per_species_population$Allele_N2
  Allele_O2 <- table_per_species_population$Allele_O2
  Allele_P2 <- table_per_species_population$Allele_P2
  Allele_Q2 <- table_per_species_population$Allele_Q2
  Allele_R2 <- table_per_species_population$Allele_R2
  Allele_S2 <- table_per_species_population$Allele_S2
  Allele_T2 <- table_per_species_population$Allele_T2
  Allele_U2 <- table_per_species_population$Allele_U2
  Allele_V2 <- table_per_species_population$Allele_V2
  Allele_W2 <- table_per_species_population$Allele_W2
  Allele_X2 <- table_per_species_population$Allele_X2
  Allele_Y2 <- table_per_species_population$Allele_Y2




  if(this_index == "heterozygosity"){
    for (j in 1:length(Allele_A1)){ # using H from Pegas package
      all_H_values <- c(all_H_values,pegas::H(c(Allele_A1[j],Allele_B1[j],Allele_C1[j],Allele_D1[j],Allele_E1[j],
                                                Allele_F1[j],Allele_G1[j],Allele_H1[j],Allele_I1[j],Allele_J1[j],
                                                Allele_K1[j],Allele_L1[j],Allele_M1[j],Allele_N1[j],Allele_O1[j],
                                                Allele_P1[j],
                                                Allele_Q1[j],
                                                Allele_R1[j],
                                                Allele_S1[j],
                                                Allele_T1[j],
                                                Allele_U1[j],
                                                Allele_V1[j],
                                                Allele_W1[j],
                                                Allele_X1[j],
                                                Allele_Y1[j],



                                                Allele_A2[j],Allele_B2[j],Allele_C2[j],Allele_D2[j],Allele_E2[j],
                                                Allele_F2[j],Allele_G2[j],Allele_H2[j],Allele_I2[j],Allele_J2[j],
                                                Allele_K2[j],Allele_L2[j],Allele_M2[j],Allele_N2[j],Allele_O2[j],
                                                Allele_P2[j],
                                                Allele_Q2[j],
                                                Allele_R2[j],
                                                Allele_S2[j],
                                                Allele_T2[j],
                                                Allele_U2[j],
                                                Allele_V2[j],
                                                Allele_W2[j],
                                                Allele_X2[j],
                                                Allele_Y2[j])))


      all_H_selection <- c(all_H_selection,pegas::H(c(Allele_A1[j],Allele_B1[j],Allele_C1[j],Allele_D1[j],Allele_E1[j],
                                                      Allele_F1[j],Allele_G1[j],Allele_H1[j],Allele_I1[j],Allele_J1[j],
                                                      Allele_K1[j],Allele_L1[j],Allele_M1[j],Allele_N1[j],Allele_O1[j],
                                                      Allele_P1[j],
                                                      Allele_Q1[j],
                                                      Allele_R1[j],
                                                      Allele_S1[j],
                                                      Allele_T1[j],
                                                      Allele_U1[j],
                                                      Allele_V1[j],
                                                      Allele_W1[j],
                                                      Allele_X1[j],
                                                      Allele_Y1[j]
      )))

      all_H_neutral <- c(all_H_neutral,pegas::H(c(Allele_A2[j],Allele_B2[j],Allele_C2[j],Allele_D2[j],Allele_E2[j],
                                                  Allele_F2[j],Allele_G2[j],Allele_H2[j],Allele_I2[j],Allele_J2[j],
                                                  Allele_K2[j],Allele_L2[j],Allele_M2[j],Allele_N2[j],Allele_O2[j],
                                                  Allele_P2[j],
                                                  Allele_Q2[j],
                                                  Allele_R2[j],
                                                  Allele_S2[j],
                                                  Allele_T2[j],
                                                  Allele_U2[j],
                                                  Allele_V2[j],
                                                  Allele_W2[j],
                                                  Allele_X2[j],
                                                  Allele_Y2[j])))
    }
  }


  if(this_index == "shannon"){
    for (j in 1:length(Allele_A1)){ # using H from Pegas package
      all_H_values <- c(all_H_values,vegan::diversity(c(Allele_A1[j],Allele_B1[j],Allele_C1[j],Allele_D1[j],Allele_E1[j],
                                                        Allele_F1[j],Allele_G1[j],Allele_H1[j],Allele_I1[j],Allele_J1[j],
                                                        Allele_K1[j],Allele_L1[j],Allele_M1[j],Allele_N1[j],Allele_O1[j],
                                                        Allele_P1[j],
                                                        Allele_Q1[j],
                                                        Allele_R1[j],
                                                        Allele_S1[j],
                                                        Allele_T1[j],
                                                        Allele_U1[j],
                                                        Allele_V1[j],
                                                        Allele_W1[j],
                                                        Allele_X1[j],
                                                        Allele_Y1[j],



                                                        Allele_A2[j],Allele_B2[j],Allele_C2[j],Allele_D2[j],Allele_E2[j],
                                                        Allele_F2[j],Allele_G2[j],Allele_H2[j],Allele_I2[j],Allele_J2[j],
                                                        Allele_K2[j],Allele_L2[j],Allele_M2[j],Allele_N2[j],Allele_O2[j],
                                                        Allele_P2[j],
                                                        Allele_Q2[j],
                                                        Allele_R2[j],
                                                        Allele_S2[j],
                                                        Allele_T2[j],
                                                        Allele_U2[j],
                                                        Allele_V2[j],
                                                        Allele_W2[j],
                                                        Allele_X2[j],
                                                        Allele_Y2[j]),index="shannon"))
      all_H_selection <- c(all_H_selection,vegan::diversity(c(Allele_A1[j],Allele_B1[j],Allele_C1[j],Allele_D1[j],Allele_E1[j],
                                                              Allele_F1[j],Allele_G1[j],Allele_H1[j],Allele_I1[j],Allele_J1[j],
                                                              Allele_K1[j],Allele_L1[j],Allele_M1[j],Allele_N1[j],Allele_O1[j],
                                                              Allele_P1[j],
                                                              Allele_Q1[j],
                                                              Allele_R1[j],
                                                              Allele_S1[j],
                                                              Allele_T1[j],
                                                              Allele_U1[j],
                                                              Allele_V1[j],
                                                              Allele_W1[j],
                                                              Allele_X1[j],
                                                              Allele_Y1[j]),index="shannon"))
      all_H_neutral <- c(all_H_neutral,vegan::diversity(c(Allele_A2[j],Allele_B2[j],Allele_C2[j],Allele_D2[j],Allele_E2[j],
                                                          Allele_F2[j],Allele_G2[j],Allele_H2[j],Allele_I2[j],Allele_J2[j],
                                                          Allele_K2[j],Allele_L2[j],Allele_M2[j],Allele_N2[j],Allele_O2[j],
                                                          Allele_P2[j],
                                                          Allele_Q2[j],
                                                          Allele_R2[j],
                                                          Allele_S2[j],
                                                          Allele_T2[j],
                                                          Allele_U2[j],
                                                          Allele_V2[j],
                                                          Allele_W2[j],
                                                          Allele_X2[j],
                                                          Allele_Y2[j]),index="shannon"))
    }
  }


  return(cbind(H_value_bothLoci=all_H_values,H_value_selection=all_H_selection,H_value_neutral=all_H_neutral))
}





per_species_summary <- function (list_species_from_cpp,time_simulated_from_output,
                                 map_temperature,filtering_rangeSize,
                                 num_sampled_cells,start_window = NULL,
                                 end_window = NULL){
  # end_window is closer to the present, so it could be 0
  all_id <- NULL
  all_parent <- NULL
  all_birth <- NULL
  all_Birth_northmost <- NULL
  all_Birth_southmost <- NULL
  all_temperature_born <- NULL
  all_saturation_grid_birth <- NULL
  all_death <- NULL
  all_rangeSize <- NULL
  all_totalArea <- NULL
  all_initial_range <- NULL
  all_swisscheese <- NULL
  all_totalIndividuals <- NULL
  all_H_selection <- NULL
  all_H_selection_sd <- NULL
  all_H_neutral <- NULL
  all_H_neutral_sd <- NULL
  all_alive <- NULL
  all_geneflow_events <- NULL
  all_trait_state <- NULL

  all_H_bothLoci <- NULL
  all_H_bothLoci_sd <- NULL
  all_northernnmost <- NULL
  all_southernmost <- NULL
  all_expansion_failure_becauseK <- NULL
  all_total_change_southernmost <- NULL
  all_total_change_northernmost <- NULL

  meanbetaFST_neutral <- NULL
  sd_betaFST_neutral <- NULL
  meanbetaFST_selection <- NULL
  sd_betaFST_selection <- NULL
  meanbetaFST_bothLoci <- NULL
  gradient_C_S_H <- NULL
  gradient_C_N_H <- NULL


  wc_fst_bothLoci <- NULL
  wc_fst_neutral <- NULL
  wc_fst_selection <- NULL

  gradient_C_S_fst <- NULL
  gradient_C_N_fst <- NULL
  reconst_S_expansion <- NULL
  reconst_N_expansion <- NULL
  penalty_distance <- NULL

  all_expansion_north_events_first_half <- NULL
  all_contraction_north_events_first_half <- NULL
  all_expansion_north_events_second_half <- NULL
  all_contraction_north_events_second_half <- NULL

  all_expansion_south_events_first_half <- NULL
  all_contraction_south_events_first_half <- NULL
  all_expansion_south_events_second_half <- NULL
  all_contraction_south_events_second_half <- NULL


  all_rate_expansion_N_early <- NULL
  all_rate_expansion_N_late <- NULL

  all_rate_contraction_N_early <- NULL
  all_rate_contraction_N_late <- NULL

  all_rate_expansion_S_early <- NULL
  all_rate_expansion_S_late <- NULL

  all_rate_contraction_S_early <- NULL
  all_rate_contraction_S_late <- NULL

  total_N_change <- NULL
  total_S_change <- NULL


  all_events_expansion_N <- NULL
  all_events_contraction_N <- NULL
  all_events_expansion_S <- NULL
  all_events_contraction_S <- NULL
  all_zed_fst <- NULL
  all_zed_H <- NULL
  all_zed_alleles <- NULL
  all_zed_ind <- NULL

  all_zed_NC_fst <- NULL
  all_zed_NC_H <- NULL
  all_zed_NC_ind <- NULL
  all_zed_NC_alleles <- NULL

  all_spp_richness <- 0
  for(i in 1:length(list_species_from_cpp)){
    cat("doing species: ",i,"\n")
    this_species <- list_species_from_cpp[[i]]
    #### new part

    be_born <- time_simulated_from_output - this_species$Birth
    if(this_species$Death != 0){ # extinct spp
      be_dead <- time_simulated_from_output - this_species$Death
    } else {
      be_dead <-  this_species$Death
      all_spp_richness <- all_spp_richness + 1
    }

    mid_age <- (be_born - be_dead) / 2
    ealy_years <- (be_born - be_dead) / 10 # the first 10% percent of the lifetime of the species
    old_years <- (be_born - be_dead) - ( (be_born - be_dead) / 10) # the last 10% percent of the lifetime of the species

    vector_change_northernmost <- this_species$change_northernmost
    vector_time_change_northernmost <- this_species$time_change_northernmost
    vector_time_change_northernmost_ageConsidered <- vector_time_change_northernmost  - (time_simulated_from_output - be_born)

    expansion_N_time <- NULL
    contraction_N_time <- NULL
    if(length(vector_change_northernmost) > 0){
      for(iijk in 1:length(vector_time_change_northernmost_ageConsidered)){
        if(vector_change_northernmost[iijk] == 1){

          expansion_N_time <- c(expansion_N_time,vector_time_change_northernmost_ageConsidered[iijk])
        } else {
          contraction_N_time <-  c(contraction_N_time,vector_time_change_northernmost_ageConsidered[iijk])
        }
      }

    }


    rate_expansion_N_early <- length(which(expansion_N_time <= ealy_years))/ealy_years
    rate_expansion_N_late <- length(which(expansion_N_time >= old_years))/old_years

    rate_contraction_N_early <- length(which(contraction_N_time <= ealy_years))/ealy_years
    rate_contraction_N_late <- length(which(contraction_N_time >= old_years))/old_years

    ###
    vector_change_southernmost <- this_species$change_southernmost
    vector_time_change_southernmost <- this_species$time_change_southernmost
    vector_time_change_southernmost_ageConsidered <- vector_time_change_southernmost  - (time_simulated_from_output - be_born)

    ##
    if(length(this_species$per_species_population) != 1){ # extinct species
      z <- matrix(this_species$per_species_population[,c(1,2)],, ncol = 2)
      if(nrow( z ) >= 3){
        # chull(z)
        # plot(z, cex = 0.5)
        # polygon(z[chull(z), ])
        my_polygon <- z[chull(z), ]
        my_polygon2 <- rbind(my_polygon,my_polygon[1,])
        polygon3 <- st_polygon(list(my_polygon2))
        totalArea <- st_area (polygon3)
      } else {
        totalArea <- NA
      }
    } else {
      totalArea <- NA
    }





    ##
    expansion_S_time <- NULL
    contraction_S_time <- NULL
    if(length(vector_change_southernmost) > 0){
      for(iijk in 1:length(vector_time_change_southernmost_ageConsidered)){
        if(vector_change_southernmost[iijk] == 1){

          expansion_S_time <- c(expansion_S_time,vector_time_change_southernmost_ageConsidered[iijk])
        } else {
          contraction_S_time <-  c(contraction_S_time,vector_time_change_southernmost_ageConsidered[iijk])
        }
      }

    }
    rate_expansion_S_early <- length(which(expansion_S_time <= ealy_years))/ealy_years
    rate_expansion_S_late <-length(which(expansion_S_time >= old_years))/old_years

    rate_contraction_S_early <- length(which(contraction_S_time <= ealy_years))/ealy_years
    rate_contraction_S_late <-length(which(contraction_S_time >= old_years))/old_years




    all_rate_expansion_N_early <- c(all_rate_expansion_N_early,rate_expansion_N_early)
    all_rate_expansion_N_late <- c(all_rate_expansion_N_late,rate_expansion_N_late)

    all_rate_contraction_N_early <- c(all_rate_contraction_N_early,rate_contraction_N_early)
    all_rate_contraction_N_late <- c(all_rate_contraction_N_late,rate_contraction_N_late)

    all_rate_expansion_S_early <- c(all_rate_expansion_S_early,rate_expansion_S_early)
    all_rate_expansion_S_late <- c(all_rate_expansion_S_late,rate_expansion_S_late)

    all_rate_contraction_S_early <- c(all_rate_contraction_S_early,rate_contraction_S_early)
    all_rate_contraction_S_late <- c(all_rate_contraction_S_late,rate_contraction_S_late)



    total_N_change <- c(total_N_change,this_species$Birth_northmost - this_species$northernnmost )

    total_S_change <- c(total_S_change,this_species$southernmost  - this_species$Birth_southmost  )


    if(is.null(expansion_N_time) == FALSE){
      all_events_expansion_N <- c(all_events_expansion_N,
                                  length(which(start_window >= expansion_N_time & expansion_N_time >= end_window)))

      all_events_contraction_N <- c(all_events_contraction_N,
                                    length(which(start_window >= contraction_N_time & contraction_N_time >= end_window)))

      all_events_expansion_S <- c(all_events_expansion_S,
                                  length(which(start_window >= expansion_S_time & expansion_S_time >= end_window)))

      all_events_contraction_S <- c(all_events_contraction_S,
                                    length(which(start_window >= contraction_S_time & contraction_S_time >= end_window)))

    } else {
      all_events_expansion_N <- c(all_events_expansion_N,
                                  NA)

      all_events_contraction_N <- c(all_events_contraction_N,
                                    NA)

      all_events_expansion_S <- c(all_events_expansion_S,
                                  NA)

      all_events_contraction_S <- c(all_events_contraction_S,
                                    NA)

    }


    ##### end of new part

    all_id <- c(all_id,this_species$ID)
    all_parent <- c(all_parent,this_species$Parent)
    all_birth <- c(all_birth,time_simulated_from_output - this_species$Birth)
    all_trait_state <- c(all_trait_state,this_species$trait_state)
    all_initial_range <- c(all_initial_range,this_species$initial_rangesize)
    all_Birth_northmost <- c(all_Birth_northmost,this_species$Birth_northmost)
    all_Birth_southmost <- c(all_Birth_southmost,this_species$Birth_southmost)
    all_expansion_failure_becauseK <- c(all_expansion_failure_becauseK,this_species$expansion_failure_becauseK)
    all_total_change_northernmost <- c(all_total_change_northernmost,this_species$total_change_northernmost)
    all_total_change_southernmost <- c(all_total_change_southernmost,this_species$total_change_southernmost)

    #all_temperature_born <- c(all_temperature_born,map_temperature[this_species$Birth_placeX,this_species$Birth_placeY])
    expansion_contraction_rates <- calculate_rates_expansion_contraction(this_species,time_simulated_from_output)

    all_expansion_north_events_first_half <- c(all_expansion_north_events_first_half,expansion_contraction_rates$expansion_north_events_first_half)
    all_contraction_north_events_first_half <- c(all_contraction_north_events_first_half, expansion_contraction_rates$contraction_north_events_first_half)
    all_expansion_north_events_second_half <- c(all_expansion_north_events_second_half,expansion_contraction_rates$expansion_north_events_second_half)
    all_contraction_north_events_second_half <- c(all_contraction_north_events_second_half,expansion_contraction_rates$contraction_north_events_second_half)

    all_expansion_south_events_first_half <- c(all_expansion_south_events_first_half,expansion_contraction_rates$expansion_south_events_first_half)
    all_contraction_south_events_first_half <- c(all_contraction_south_events_first_half,expansion_contraction_rates$contraction_south_events_first_half)
    all_expansion_south_events_second_half <- c(all_expansion_south_events_second_half,expansion_contraction_rates$expansion_south_events_second_half)
    all_contraction_south_events_second_half <- c(all_contraction_south_events_second_half,expansion_contraction_rates$contraction_south_events_second_half)

    all_saturation_grid_birth <- c(all_saturation_grid_birth,this_species$saturation_grid_birth)
    all_rangeSize <- c(all_rangeSize,this_species$RangeSize)
    all_totalArea <- c(all_totalArea,totalArea)
    all_swisscheese <- c(all_swisscheese,this_species$swisscheese)
    all_geneflow_events <- c(all_geneflow_events, this_species$geneflow_events)
    if(class(this_species$per_species_population)[1] != "matrix"){
      all_alive <- c(all_alive,"extinct")
      all_death <- c(all_death,time_simulated_from_output - this_species$Death)
      all_totalIndividuals <- c(all_totalIndividuals,0)
      all_H_selection <- c(all_H_selection,0)
      all_H_selection_sd <- c(all_H_selection_sd,0)
      all_H_neutral <- c(all_H_neutral,0)
      all_H_neutral_sd <- c(all_H_neutral_sd,0)
      all_H_bothLoci <- c(all_H_bothLoci,0)
      all_H_bothLoci_sd <- c(all_H_bothLoci_sd,0)

    } else {
      this_species$per_species_population <- cbind(this_species$per_species_population,compute_H_values(this_species$per_species_population,this_index="shannon"))

      all_alive <- c(all_alive,"alive")
      all_totalIndividuals <- c(all_totalIndividuals,sum(this_species$per_species_population[,3]))

      all_H_selection <- c(all_H_selection,median(this_species$per_species_population[,55],na.rm=T))
      all_H_selection_sd <- c(all_H_selection_sd, sd(this_species$per_species_population[,55],na.rm=T))
      all_H_neutral <- c(all_H_neutral,median(this_species$per_species_population[,56],na.rm=T))
      all_H_neutral_sd <- c(all_H_neutral_sd, sd(this_species$per_species_population[,56],na.rm=T))

      all_H_bothLoci <- c(all_H_bothLoci,median(this_species$per_species_population[,54],na.rm=T))
      all_H_bothLoci_sd <- c(all_H_bothLoci_sd, sd(this_species$per_species_population[,54],na.rm=T))
      all_death <- c(all_death,this_species$Death)
    }

    table_to_do <- as.data.frame(list_species_from_cpp[[i]]$per_species_population)
    # range_size_threshold <- 15
    # sampling_runs <- 3 # to vary the first population to start the spatial sampling
    # sampling_divided_by <- 1 #a denominator for the sampling : 1 means all populations within a radius are included, 2 is half
    zeds <- wrapper_make_gda(do_centroid_pop = FALSE,table_to_do,range_size_threshold = filtering_rangeSize,sampling_runs = 3,sampling_divided_by=2)
  #non centroid = NC
      NC_zeds <- wrapper_make_gda(do_centroid_pop = TRUE,table_to_do,range_size_threshold = filtering_rangeSize,sampling_runs = 3,sampling_divided_by=2)
    sampled_table_and_centroid <- zeds$sampled_table_and_centroid


    all_zed_NC_fst <- c(all_zed_NC_fst, NC_zeds$zed_fst)
    all_zed_NC_H <- c(all_zed_NC_H, NC_zeds$zed_H)
    all_zed_NC_ind <- c(all_zed_NC_ind,NC_zeds$zed_ind)
    all_zed_NC_alleles <- c(all_zed_NC_alleles,NC_zeds$zed_alleles)


    all_zed_fst <- c(all_zed_fst, zeds$zed_fst)
    all_zed_H <- c(all_zed_H, zeds$zed_H)
    all_zed_ind <- c(all_zed_ind,zeds$zed_ind)
    all_zed_alleles <- c(all_zed_alleles,zeds$zed_alleles)

    if(nrow(table_to_do) > 4){
      computed_betas_neutral  <- compute_betaPopulationsSpecificFST(table_to_do,pairwise=F,locus = "neutral")
      computed_betas_selection  <- compute_betaPopulationsSpecificFST(table_to_do,pairwise=F,locus = "selection")
      computed_betas_bothLoci  <- compute_betaPopulationsSpecificFST(table_to_do,pairwise=F,locus = "together")
      table_to_do2 <- cbind(table_to_do,betaFST = as.numeric(computed_betas_neutral$computed_betas$betaiovl))
      if(any(is.na(table_to_do2$betaFST))){
        table_to_do_NAfree <- table_to_do2[-which(is.na(table_to_do2$betaFST)),]
      } else {
        table_to_do_NAfree <- table_to_do2
      }

      wc_fst_bothLoci <-c(wc_fst_bothLoci,computed_betas_bothLoci$wc_fst$FST)
        wc_fst_neutral <- c(wc_fst_neutral,computed_betas_neutral$wc_fst$FST)
      wc_fst_selection <- c(wc_fst_selection,computed_betas_selection$wc_fst$FST)

      meanbetaFST_bothLoci <- c(meanbetaFST_bothLoci,
                                mean(as.numeric(computed_betas_bothLoci$computed_betas$betaiovl),na.rm=T))
      meanbetaFST_neutral <- c(meanbetaFST_neutral,
                               mean(as.numeric(computed_betas_neutral$computed_betas$betaiovl),na.rm=T))
      sd_betaFST_neutral <- c(sd_betaFST_neutral,
                              sd(as.numeric(computed_betas_neutral$computed_betas$betaiovl),na.rm=T))

      meanbetaFST_selection <- c(meanbetaFST_selection,
                                 mean(as.numeric(computed_betas_selection$computed_betas$betaiovl),na.rm=T))
      sd_betaFST_selection <- c(sd_betaFST_selection,
                                sd(as.numeric(computed_betas_selection$computed_betas$betaiovl),na.rm=T))

    } else {
      table_to_do_NAfree <- data.frame()
      meanbetaFST_neutral <- c(meanbetaFST_neutral,NA)
      meanbetaFST_bothLoci <- c(meanbetaFST_bothLoci,NA)
      sd_betaFST_neutral <- c(sd_betaFST_neutral,NA)
      meanbetaFST_selection <- c(meanbetaFST_selection,NA)
      sd_betaFST_selection <- c(sd_betaFST_selection,NA)

      wc_fst_bothLoci <-c(wc_fst_bothLoci,NA)
      wc_fst_neutral <- c(wc_fst_neutral,NA)
      wc_fst_selection <- c(wc_fst_selection,NA)

    }

    if(list_species_from_cpp[[i]]$Death == 0 && nrow(table_to_do_NAfree) > filtering_rangeSize){# to skip extinct and one-population species


      #sampled_table_and_centroid <- prepare_sampled_table_and_find_centroid(table_to_do,num_sampled_cells,only_edges = T,sample_it = T, strict_geographic = T)

      centroid_population_distance_abundance <- sampled_table_and_centroid$centroid_population_distance_abundance

      longest_vert_line <- names(table(table_to_do_NAfree$X)[which(max(table(table_to_do_NAfree$X))==table(table_to_do_NAfree$X))][1])
      longest_vert_line <- as.numeric(longest_vert_line)
      # sampled_table_to_do <- table_to_do_NAfree[which(table_to_do_NAfree$X == longest_vert_line ),]
      sampled_table_to_do <- table_to_do_NAfree[sample(1:nrow(table_to_do_NAfree),filtering_rangeSize,replace=FALSE),]


      computed_pairwisebetas  <- compute_betaPopulationsSpecificFST(sampled_table_to_do,pairwise=TRUE,locus="neutral")
      colnames(computed_pairwisebetas$computed_betas) <- rownames(sampled_table_to_do)
      rownames(computed_pairwisebetas$computed_betas) <- rownames(sampled_table_to_do)
      nj_res <- NULL
      try(nj_res <- njs(as.dist(computed_pairwisebetas$computed_betas))
          ,silent=TRUE)

      ## IMportant: switch off of this analysis because it might be slow to compute the euclidean
      turn_this_analysis_off <- TRUE

      if(is.null(nj_res) ||  turn_this_analysis_off){
        reconst_S_expansion <- c(reconst_S_expansion,NA)
        reconst_N_expansion <- c(reconst_N_expansion,NA)
        penalty_distance <- c(penalty_distance,NA)
      } else {
        order_colonization_pop_output <- order_colonization_pop_route(sampled_table_to_do,nj_res)

        north_expansion <- 0
        south_expansion <- 0
        for(jij in 1:(length(order_colonization_pop_output$the_order) - 1)){
          pop_orign <- sampled_table_to_do[which(order_colonization_pop_output$the_order[jij] == rownames(sampled_table_to_do)),]
          p_goingto<- sampled_table_to_do[which(order_colonization_pop_output$the_order[jij + 1] == rownames(sampled_table_to_do)),]
          if(pop_orign$Y > p_goingto$Y){
            south_expansion <- south_expansion + 1
          } else{
            north_expansion <- north_expansion + 1
          }
        }

        reconst_S_expansion <- c(reconst_S_expansion,south_expansion/nrow(table_to_do2))
        reconst_N_expansion <-  c(reconst_N_expansion,north_expansion/nrow(table_to_do2))
        penalty_distance <-  c(penalty_distance,order_colonization_pop_output$penalty_distance)

      }


      fst_centroid <- computed_betas_neutral$betaiovl[which(names(computed_betas_neutral$betaiovl) == as.character(rownames(centroid_population_distance_abundance)))]
      fst_centroid <- NULL
      for(ipo in 1:nrow(centroid_population_distance_abundance)){
        fst_centroid <- c(fst_centroid,computed_betas_neutral$betaiovl[which(names(computed_betas_neutral$betaiovl)
                                                                             == as.character(rownames(centroid_population_distance_abundance))[ipo])])

      }
      fst_centroid <- mean(fst_centroid,na.rm=T)
      most_north <- sampled_table_and_centroid$most_north
      most_south <- sampled_table_and_centroid$most_south
      fst_north_edge <- NULL
      for(ipo in 1:nrow(most_north)){
        fst_north_edge <- c(fst_north_edge,computed_betas_neutral$betaiovl[which(names(computed_betas_neutral$betaiovl)
                                                                                 == as.character(rownames(most_north))[ipo])])

      }
      fst_north_edge <- mean(fst_north_edge,na.rm=T)
      fst_south_edge <- NULL
      for(ipo in 1:nrow(most_south)){
        fst_south_edge <- c(fst_south_edge,computed_betas_neutral$betaiovl[which(names(computed_betas_neutral$betaiovl)
                                                                                 == as.character(rownames(most_south))[ipo])])

      }
      fst_south_edge <- mean(fst_south_edge,na.rm=T)


      H_centroid <- mean(this_species$per_species_population[as.numeric(rownames(centroid_population_distance_abundance)),54],na.rm=T )


      H_north_edge <- mean(this_species$per_species_population[as.numeric(rownames(most_north)),54],na.rm=T )

      H_south_edge <- mean(this_species$per_species_population[as.numeric(rownames(most_south)),54],na.rm=T )


      gradient_C_S_fst <- c(gradient_C_S_fst,mean(fst_centroid - fst_south_edge,na.rm=T))
      gradient_C_N_fst <- c(gradient_C_N_fst,mean(fst_centroid - fst_north_edge,na.rm=T))

      gradient_C_S_H <- c(gradient_C_S_H,mean(H_centroid - H_south_edge,na.rm=T))
      gradient_C_N_H <- c(gradient_C_N_H,mean(H_centroid - H_north_edge,na.rm=T))




      fst_one_spp <- cbind(table_to_do,betaFST = as.numeric(computed_betas_neutral$computed_betas$betaiovl))
    } else {

      gradient_C_S_H <- c(gradient_C_S_H,NA)
      gradient_C_N_H <- c(gradient_C_N_H,NA)
      gradient_C_S_fst <- c( gradient_C_S_fst ,NA)
      gradient_C_N_fst <- c(gradient_C_N_fst,NA)
      reconst_S_expansion <- c(reconst_S_expansion,NA)
      reconst_N_expansion <- c(reconst_N_expansion,NA)
      penalty_distance <- c(penalty_distance,NA)
    }



  }




  all_spp_richness <- rep(all_spp_richness,length(all_id))




  species_table <- data.frame( id = all_id,
                               parent = all_parent,
                               birth = all_birth,
                               death = all_death,
                               Range = all_rangeSize,
                               totalArea = all_totalArea,
                               initial_range = all_initial_range,
                               trait_state = all_trait_state,
                               swisscheese = all_swisscheese,
                               TotInd = all_totalIndividuals,
                               H_bothLoci = all_H_bothLoci,
                               H_bothLoci_sd = all_H_bothLoci_sd,
                               H_selection = all_H_selection,
                               H_selection_sd = all_H_selection_sd,
                               H_neutral = all_H_neutral,
                               H_neutral_sd = all_H_neutral_sd,
                               meanbetaFST_neutral = meanbetaFST_neutral,
                               sd_betaFST_neutral = sd_betaFST_neutral,
                               meanbetaFST_selection = meanbetaFST_selection,
                               sd_betaFST_selection = sd_betaFST_selection,
                               meanbetaFST_bothLoci = meanbetaFST_bothLoci,
                               wc_fst_bothLoci = wc_fst_bothLoci,
                               wc_fst_neutral = wc_fst_neutral,
                               wc_fst_selection = wc_fst_selection,
                               gradient_C_S_H = gradient_C_S_H,
                               gradient_C_N_H = gradient_C_N_H,
                               gradient_C_S_fst = gradient_C_S_fst,
                               gradient_C_N_fst = gradient_C_N_fst,
                               reconst_S_expansion = reconst_S_expansion,
                               reconst_N_expansion = reconst_N_expansion,
                               penalty_distance  = penalty_distance,
                               saturation_grid_birth = all_saturation_grid_birth,
                               spp_richness = all_spp_richness,
                               alive = all_alive,
                               Birth_northmost = all_Birth_northmost,
                               Birth_southmost = all_Birth_southmost,
                               expansion_failure_becauseK = all_expansion_failure_becauseK,
                               total_change_southernmost = all_total_change_southernmost,
                               total_change_northernmost = all_total_change_northernmost,
                               #temperature_born = all_temperature_born,
                               geneflow_events = all_geneflow_events,
                               rate_expansion_N_early = all_rate_expansion_N_early,
                               rate_expansion_N_late = all_rate_expansion_N_late,
                               rate_contraction_N_early = all_rate_contraction_N_early,
                               rate_contraction_N_late = all_rate_contraction_N_late,
                               rate_expansion_S_early = all_rate_expansion_S_early,
                               rate_expansion_S_late = all_rate_expansion_S_late,
                               rate_contraction_S_early = all_rate_contraction_S_early,
                               rate_contraction_S_late = all_rate_contraction_S_late,
                               total_N_change = total_N_change,
                               total_S_change = total_S_change,
                               events_expansion_N = all_events_expansion_N,
                               events_contraction_N = all_events_contraction_N,
                               events_expansion_S = all_events_expansion_S,
                               events_contraction_S = all_events_contraction_S,
                               all_zed_fst = all_zed_fst,
                               all_zed_H = all_zed_H,
                               all_zed_alleles = all_zed_alleles,
                               all_zed_ind = all_zed_ind,
                               all_zed_NC_fst = all_zed_NC_fst,
                               all_zed_NC_H = all_zed_NC_H,
                               all_zed_NC_ind = all_zed_NC_ind,
                               all_zed_NC_alleles = all_zed_NC_alleles)
  return (species_table)
}
###############################################################################

calculate_rates_expansion_contraction <- function(this_species,time_simulated_from_output){

  species_birth <- this_species$Birth
  middle_age_species <- species_birth  + ((time_simulated_from_output -  species_birth) / 2)


  expansion_north_events_first_half <- 0
  contraction_north_events_first_half <- 0
  expansion_north_events_second_half <- 0
  contraction_north_events_second_half <- 0
  if(is.na(this_species$change_northernmost[1]) == FALSE){ # no expansion events
    for(i in 1:length( this_species$time_change_northernmost)){

      if(this_species$time_change_northernmost[i] < middle_age_species){  # first half?


        if(this_species$change_northernmost[i] == 1){
          expansion_north_events_first_half <- expansion_north_events_first_half + 1
        }
        if(this_species$change_northernmost[i] == -1){
          contraction_north_events_first_half <- contraction_north_events_first_half + 1
        }



      } else { # second half then

        if(this_species$change_northernmost[i] == 1){
          expansion_north_events_second_half <- expansion_north_events_second_half + 1
        }
        if(this_species$change_northernmost[i] == -1){
          contraction_north_events_second_half <- contraction_north_events_second_half + 1
        }

      }

    }
  }


  ##

  expansion_south_events_first_half <- 0
  contraction_south_events_first_half <- 0
  expansion_south_events_second_half <- 0
  contraction_south_events_second_half <- 0
  if(is.na(this_species$change_southernmost[1]) == FALSE){ # no contraction events
    for(i in 1:length( this_species$time_change_southernmost)){

      if(this_species$time_change_southernmost[i] < middle_age_species){  # first half?


        if(this_species$change_southernmost[i] == 1){
          expansion_south_events_first_half <- expansion_south_events_first_half + 1
        }
        if(this_species$change_southernmost[i] == -1){
          contraction_south_events_first_half <- contraction_south_events_first_half + 1
        }



      } else { # second half then

        if(this_species$change_southernmost[i] == 1){
          expansion_south_events_second_half <- expansion_south_events_second_half + 1
        }
        if(this_species$change_southernmost[i] == -1){
          contraction_south_events_second_half <- contraction_south_events_second_half + 1
        }

      }

    }
  }

  return(list(expansion_north_events_first_half = expansion_north_events_first_half,
              contraction_north_events_first_half = contraction_north_events_first_half,
              expansion_north_events_second_half =  expansion_north_events_second_half,
              contraction_north_events_second_half = contraction_north_events_second_half,
              expansion_south_events_first_half = expansion_south_events_first_half,
              contraction_south_events_first_half = contraction_south_events_first_half,
              expansion_south_events_second_half = expansion_south_events_second_half,
              contraction_south_events_second_half = contraction_south_events_second_half ))



}





build_tree_and_sister_pairs <- function(time_simulated_from_output,
                                        species_table,
                                        forced_crown_age = F){

  if(class(species_table) != "data.frame"){# it is the list of species straight from cpp

    Ltable_to_do <- NULL
    for(i in 1:length(species_table)){
      take_this_sp <- species_table[[i]]
      #cat( take_this_sp$Birth, "\n")
      if(class(take_this_sp$per_species_population)[1] != "matrix"){
        pre_death <- time_simulated_from_output - take_this_sp$Death

      } else {
        pre_death <- take_this_sp$Death
      }
      Ltable_to_do <- data.frame(rbind(Ltable_to_do,cbind(birth = time_simulated_from_output - take_this_sp$Birth,
                                                          parent = take_this_sp$Parent,
                                                          id = take_this_sp$ID,
                                                          death = pre_death )))
    }
    Ltable_to_do[which(Ltable_to_do$death == 0),4] = -1

    phyloTree <- L2phylo(Ltable_to_do,dropextinct = T)
    sister_pairs <- NULL

  } else {

    time_for_Ltable <- time_simulated_from_output
    Ltable_to_do <- as.data.frame(cbind(birth=species_table$birth,parent=species_table$parent,id=species_table$id,death=species_table$death))
    # Ltable_to_do$birth <- time_simulated_from_output - Ltable_to_do$birth
    # Ltable_to_do$death <- time_simulated_from_output - Ltable_to_do$death
    Ltable_to_do[which(Ltable_to_do$death == 0),4] = -1

    phyloTree <- L2phylo(Ltable_to_do,dropextinct = T)

    if(forced_crown_age != FALSE){
      phyloTree <- rescaleTree(phyloTree,forced_crown_age) # Rescalete tree
    }

    ## starting here
    dd<-lapply(1:phyloTree$Nnode+Ntip(phyloTree),function(n,t)
      Descendants(t,n)[[1]],t=phyloTree)
    nodes<-c(1:phyloTree$Nnode+Ntip(phyloTree))[which(sapply(dd,
                                                             length)==2)]
    sisters <- t(sapply(nodes,function(n,t)
      t$tip.label[Descendants(t,n)[[1]]],t=phyloTree))
    rownames(sisters) <- nodes

    range_sister_A <- NULL
    range_sister_B <- NULL
    temp_born_sister_A <- NULL
    temp_born_sister_B <- NULL


    sister_A_expansion_north_events_first_half <- NULL
    sister_A_contraction_north_events_first_half <- NULL
    sister_A_expansion_north_events_second_half <- NULL
    sister_A_contraction_north_events_second_half <- NULL
    sister_A_expansion_south_events_first_half <- NULL
    sister_A_contraction_south_events_first_half <- NULL
    sister_A_expansion_south_events_second_half <- NULL
    sister_A_contraction_south_events_second_half <- NULL


    sister_B_expansion_north_events_first_half <- NULL
    sister_B_contraction_north_events_first_half <- NULL
    sister_B_expansion_north_events_second_half <- NULL
    sister_B_contraction_north_events_second_half <- NULL
    sister_B_expansion_south_events_first_half <- NULL
    sister_B_contraction_south_events_first_half <- NULL
    sister_B_expansion_south_events_second_half <- NULL
    sister_B_contraction_south_events_second_half <- NULL

    for(i in 1:nrow(sisters)){

      sisters[i,1] <- as.numeric (str_split(sisters[i,1],"t")[[1]][2])
      sisters[i,2] <- as.numeric (str_split(sisters[i,2],"t")[[1]][2])


    }


    nodes_age <- tree.age(phyloTree, order = "past", fossil = F, digits = 3)

    divergence_time_sisters <- NULL
    for (i in 1:nrow(sisters)){

      divergence_time_sisters <- c(divergence_time_sisters,
                                   nodes_age[which(rownames(sisters)[i] == nodes_age$elements),1])
    }



    for(i in 1:nrow(sisters)){

      range_sister_A <- c(range_sister_A,species_table[which(species_table$id ==   sisters[i,1]),]$Range)
      range_sister_B <- c(range_sister_B,species_table[which(species_table$id ==   sisters[i,2]),]$Range)


      sister_A_expansion_north_events_first_half <- c(sister_A_expansion_north_events_first_half,
                                                      species_table[which(species_table$id ==  sisters[i,1]),]$expansion_north_events_first_half)


      sister_A_contraction_north_events_first_half <- c(sister_A_contraction_north_events_first_half,
                                                        species_table[which(species_table$id ==  sisters[i,1]),]$contraction_north_events_first_half)
      sister_A_expansion_north_events_second_half <- c(sister_A_expansion_north_events_second_half,
                                                       species_table[which(species_table$id ==  sisters[i,1]),]$expansion_north_events_second_half)
      sister_A_contraction_north_events_second_half <- c(sister_A_contraction_north_events_second_half,
                                                         species_table[which(species_table$id ==  sisters[i,1]),]$contraction_north_events_second_half)
      sister_A_expansion_south_events_first_half <- c(sister_A_expansion_south_events_first_half,
                                                      species_table[which(species_table$id ==  sisters[i,1]),]$expansion_south_events_first_half)
      sister_A_contraction_south_events_first_half <- c(sister_A_contraction_south_events_first_half,
                                                        species_table[which(species_table$id ==  sisters[i,1]),]$contraction_south_events_first_half)
      sister_A_expansion_south_events_second_half <- c(sister_A_expansion_south_events_second_half,
                                                       species_table[which(species_table$id ==  sisters[i,1]),]$expansion_south_events_second_half)
      sister_A_contraction_south_events_second_half <- c(sister_A_contraction_south_events_second_half,
                                                         species_table[which(species_table$id ==  sisters[i,1]),]$contraction_south_events_second_half)


      sister_B_expansion_north_events_first_half <- c(sister_B_expansion_north_events_first_half,
                                                      species_table[which(species_table$id ==  sisters[i,2]),]$expansion_north_events_first_half)
      sister_B_contraction_north_events_first_half <- c(sister_B_contraction_north_events_first_half,
                                                        species_table[which(species_table$id ==  sisters[i,2]),]$contraction_north_events_first_half)
      sister_B_expansion_north_events_second_half <- c(sister_B_expansion_north_events_second_half,
                                                       species_table[which(species_table$id ==  sisters[i,2]),]$expansion_north_events_second_half)

      sister_B_contraction_north_events_second_half <- c(sister_B_contraction_north_events_second_half,
                                                         species_table[which(species_table$id ==  sisters[i,2]),]$contraction_north_events_second_half)
      sister_B_expansion_south_events_first_half <- c(sister_B_expansion_south_events_first_half,
                                                      species_table[which(species_table$id ==  sisters[i,2]),]$expansion_south_events_first_half)
      sister_B_contraction_south_events_first_half <- c(sister_B_contraction_south_events_first_half,
                                                        species_table[which(species_table$id ==  sisters[i,2]),]$contraction_south_events_first_half)
      sister_B_expansion_south_events_second_half <- c(sister_B_expansion_south_events_second_half,
                                                       species_table[which(species_table$id ==  sisters[i,2]),]$expansion_south_events_second_half)
      sister_B_contraction_south_events_second_half <- c(sister_B_contraction_south_events_second_half,
                                                         species_table[which(species_table$id ==  sisters[i,2]),]$contraction_south_events_second_half)


      temp_born_sister_A <- c(temp_born_sister_A,species_table[which(species_table$id ==   sisters[i,1]),]$temperature_born)
      temp_born_sister_B <- c(temp_born_sister_B,species_table[which(species_table$id ==   sisters[i,2]),]$temperature_born)
    }




    sister_pairs <- as.data.frame(cbind(sisters,divergence_time_sisters,range_sister_A,range_sister_B,
                                        temp_born_sister_A,temp_born_sister_B,
                                        sister_A_expansion_north_events_first_half,
                                        sister_A_contraction_north_events_first_half,
                                        sister_A_expansion_north_events_second_half ,
                                        sister_A_contraction_north_events_second_half,
                                        sister_A_expansion_south_events_first_half,
                                        sister_A_contraction_south_events_first_half,
                                        sister_A_expansion_south_events_second_half,
                                        sister_A_contraction_south_events_second_half,


                                        sister_B_expansion_north_events_first_half,
                                        sister_B_contraction_north_events_first_half,
                                        sister_B_expansion_north_events_second_half,
                                        sister_B_contraction_north_events_second_half,
                                        sister_B_expansion_south_events_first_half,
                                        sister_B_contraction_south_events_first_half,
                                        sister_B_expansion_south_events_second_half ,
                                        sister_B_contraction_south_events_second_half))


  }
  return(list(phyloTree = phyloTree,
              sister_pairs = sister_pairs ))
}

make_richness_map <- function(x_max,y_max,list_species_from_cpp){

  x_coordinates <- NULL
  y_coordinates <- NULL
  richness <- NULL
  abundance <- NULL
  for (i in 1:x_max){
    for (ii in 1:y_max){
      x_coordinates <- c(x_coordinates, i)
      y_coordinates <- c(y_coordinates,ii)
      richness <- c(richness,0)
      abundance <- c(abundance,0)
    }
  }

  all_coordinates <- as.data.frame(cbind(X=x_coordinates,Y=y_coordinates,richness=richness,abundance=abundance))

  for(j in 1:length(list_species_from_cpp)){

    this_species <- as.data.frame(list_species_from_cpp[[j]]$per_species_population)

    for(jj in 1:nrow(this_species)){

      x_list <- which(this_species$X[jj] == all_coordinates$X)
      y_list <- which(this_species$Y[jj] == all_coordinates$Y)

      all_coordinates[intersect(y_list,x_list),3] <- all_coordinates[intersect(y_list,x_list),3] + 1
      all_coordinates[intersect(y_list,x_list),4] <- all_coordinates[intersect(y_list,x_list),4]  + this_species$Pop_size[[jj]]
    }
  }

  return(all_coordinates)
}

draw_global_map <- function(to_richness_map,type_global_map){

  # to flip the map
  max_y <- max(to_richness_map$y) + 1

  to_richness_map$y = max_y - to_richness_map$y



  if(type_global_map == "richness"){

    ggp1 <- ggplot(to_richness_map, aes(x, y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = richness ))+
      ggtitle("richness")
  }

  if(type_global_map == "abundances"){
    ggp1 <- ggplot(to_richness_map, aes(x, y)) +                           # Create heatmap with ggplot2
      geom_tile(aes(fill = abundance ))+
      ggtitle("abundance")
  }
  return(ggp1)
}

compute_moran_spatial_correlation <- function(list_species_from_cpp,species_to_do,allele_to_explore = NULL){

  if(all(allele_to_explore != c("A1","B1","C1","D1","E1","F1","G1","H1","I1","J1","K1","L1","M1","N1","O1",
                                "A2","B2","C2","D2","E2","F2","G2","H2","I2","J2","K2","L2","M2","N2","O2"))){
    stop("your allele_to_explore is wrong")
  }

  if(allele_to_explore == "A1"){
    allele_to_explore_column <- 4
  }
  if(allele_to_explore == "B1"){
    allele_to_explore_column <- 5
  }
  if(allele_to_explore == "C1"){
    allele_to_explore_column <- 6
  }
  if(allele_to_explore == "D1"){
    allele_to_explore_column <- 7
  }
  if(allele_to_explore == "E1"){
    allele_to_explore_column <- 8
  }
  if(allele_to_explore == "F1"){
    allele_to_explore_column <- 9
  }
  if(allele_to_explore == "G1"){
    allele_to_explore_column <- 10
  }
  if(allele_to_explore == "H1"){
    allele_to_explore_column <- 11
  }
  if(allele_to_explore == "I1"){
    allele_to_explore_column <- 12
  }
  if(allele_to_explore == "J1"){
    allele_to_explore_column <- 13
  }
  if(allele_to_explore == "K1"){
    allele_to_explore_column <- 14
  }
  if(allele_to_explore == "L1"){
    allele_to_explore_column <- 15
  }
  if(allele_to_explore == "M1"){
    allele_to_explore_column <- 16
  }
  if(allele_to_explore == "N1"){
    allele_to_explore_column <- 17
  }
  if(allele_to_explore == "O1"){
    allele_to_explore_column <- 18
  }
  if (allele_to_explore == "P1") {
    allele_to_explore_column <- 19
  }
  if (allele_to_explore == "Q1") {
    allele_to_explore_column <- 20
  }
  if (allele_to_explore == "R1") {
    allele_to_explore_column <- 21
  }
  if (allele_to_explore == "S1") {
    allele_to_explore_column <- 22
  }
  if (allele_to_explore == "T1") {
    allele_to_explore_column <- 23
  }
  if (allele_to_explore == "U1") {
    allele_to_explore_column <- 24
  }
  if (allele_to_explore == "V1") {
    allele_to_explore_column <- 25
  }
  if (allele_to_explore == "W1") {
    allele_to_explore_column <- 26
  }
  if (allele_to_explore == "X1") {
    allele_to_explore_column <- 27
  }
  if (allele_to_explore == "Y1") {
    allele_to_explore_column <- 28
  }

  ##
  if (allele_to_explore == "B2") {
    allele_to_explore_column <- 30
  }
  if (allele_to_explore == "C2") {
    allele_to_explore_column <- 31
  }
  if (allele_to_explore == "D2") {
    allele_to_explore_column <- 32
  }
  if (allele_to_explore == "E2") {
    allele_to_explore_column <- 33
  }
  if (allele_to_explore == "F2") {
    allele_to_explore_column <- 34
  }
  if (allele_to_explore == "G2") {
    allele_to_explore_column <- 35
  }
  if (allele_to_explore == "H2") {
    allele_to_explore_column <- 36
  }
  if (allele_to_explore == "I2") {
    allele_to_explore_column <- 37
  }
  if (allele_to_explore == "J2") {
    allele_to_explore_column <- 38
  }
  if (allele_to_explore == "K2") {
    allele_to_explore_column <- 39
  }
  if (allele_to_explore == "L2") {
    allele_to_explore_column <- 40
  }
  if (allele_to_explore == "M2") {
    allele_to_explore_column <- 41
  }
  if (allele_to_explore == "N2") {
    allele_to_explore_column <- 42
  }
  if (allele_to_explore == "O2") {
    allele_to_explore_column <- 43
  }
  if (allele_to_explore == "P2") {
    allele_to_explore_column <- 44
  }
  if (allele_to_explore == "Q2") {
    allele_to_explore_column <- 45
  }
  if (allele_to_explore == "R2") {
    allele_to_explore_column <- 46
  }
  if (allele_to_explore == "S2") {
    allele_to_explore_column <- 47
  }
  if (allele_to_explore == "T2") {
    allele_to_explore_column <- 48
  }
  if (allele_to_explore == "U2") {
    allele_to_explore_column <- 49
  }
  if (allele_to_explore == "V2") {
    allele_to_explore_column <- 50
  }
  if (allele_to_explore == "W2") {
    allele_to_explore_column <- 51
  }
  if (allele_to_explore == "X2") {
    allele_to_explore_column <- 52
  }
  if (allele_to_explore == "Y2") {
    allele_to_explore_column <- 53
  }

  table_to_spatial_structure <- list_species_from_cpp[[species_to_do]]$per_species_population[,c(1,2,3,allele_to_explore_column)]
  table_to_spatial_structure <- cbind(table_to_spatial_structure, table_to_spatial_structure[,4]/table_to_spatial_structure[,3])

  table_to_spatial_structure <- table_to_spatial_structure[which(table_to_spatial_structure[,5] != 0),]

  geo_dista <- as.matrix(dist(cbind(table_to_spatial_structure[,1], table_to_spatial_structure[,2])))


  geo_dista_inv <- 1/geo_dista
  diag(geo_dista_inv) <- 0
  return(Moran.I(table_to_spatial_structure[,5],geo_dista_inv,scaled=T))


}


prepare_sampled_table_and_find_centroid <- function(table_to_do,num_sampled_cells,only_edges,sample_it,strict_geographic){

  table_to_plot <- table_to_do
  if(only_edges && sample_it == FALSE){
    stop("it cannot only_edges==TRUE && sample_it == FALSE")
  }
  species_max_x <- max(table_to_plot$X)
  species_max_y <- max(table_to_plot$Y)

  species_min_x <- min(table_to_plot$X)
  species_min_y <- min(table_to_plot$Y)

  centroid_minmax <- c(mean(c(species_max_x,species_min_x)),mean(c(species_max_y,species_min_y)))


  #eucl_dist_pairwise <- distance(cbind(table_to_plot$X,table_to_plot$Y), method = "euclidean", use.row.names = TRUE) # perhaps incorrect


  # this is the distance from the theoretical centroid to each population
  eucl_dist_pairwise <- distance(rbind(as.numeric(centroid_minmax),cbind(table_to_plot$X,table_to_plot$Y)), method = "euclidean", use.row.names = TRUE)

  eucl_dist_pairwise <- eucl_dist_pairwise[,1][-1] # removing the 0 which was the distance of theoretical centroid to itself





  centroid_euclid <- table_to_plot[which(min(eucl_dist_pairwise) == eucl_dist_pairwise),1:2]




  eucl_popinfo <- cbind(table_to_plot,eucl_dist_pairwise)
  eucl_popinfo <- cbind(eucl_popinfo[order(eucl_popinfo$Pop_size,decreasing=T),],rank_abundance=1:nrow(eucl_popinfo))
  eucl_popinfo <- cbind(eucl_popinfo[order(eucl_popinfo$eucl_dist_pairwise,decreasing=F),],rank_distance=1:nrow(eucl_popinfo))

  eucl_popinfo <- cbind(eucl_popinfo,total_ranking = eucl_popinfo$rank_abundance + eucl_popinfo$rank_distance)


  centroid_population_distance_abundance <- eucl_popinfo[order(eucl_popinfo$total_ranking,decreasing=F),]

  # centroid_population_distance_abundance <- eucl_popinfo[which(min(eucl_popinfo$total_ranking)==eucl_popinfo$total_ranking),]
  # centroid_population_distance_abundance <- centroid_population_distance_abundance[1,]

  if(strict_geographic){

    centroid_population_distance_abundance <- eucl_popinfo[order(eucl_popinfo$rank_distance,decreasing=F),]
    centroid_population_distance_abundance <- centroid_population_distance_abundance[1:num_sampled_cells,]
  } else {
    centroid_population_distance_abundance <- eucl_popinfo[order(eucl_popinfo$total_ranking,decreasing=F),]
    centroid_population_distance_abundance <- centroid_population_distance_abundance[1:num_sampled_cells,]
  }


  id_row_cetroid <- rownames(centroid_population_distance_abundance)
  table_to_plot[id_row_cetroid,]


  sample_table_to_plot <- NULL
  most_north <- table_to_plot[order(table_to_plot$Y,decreasing=T),]
  most_north <- most_north[sample(1:(num_sampled_cells + 1),num_sampled_cells,replace=F),]

  most_south <- table_to_plot[order(table_to_plot$Y,decreasing=F),]
  most_south <- most_south[sample(1:(num_sampled_cells + 1),num_sampled_cells,replace=F),]


  if(only_edges){
    sample_table_to_plot <- rbind(most_north,
                                  most_south)
  } else {
    middle <- table_to_plot[order(table_to_plot$Y,decreasing=T),]
    vector_tomatch_id <- c(which(centroid_population_distance_abundance$X == middle$X),which(centroid_population_distance_abundance$Y == middle$Y ))
    id_row_cetroid_in_middle_table <- vector_tomatch_id[which(duplicated(vector_tomatch_id))]

    middle_to_north <- middle[1:id_row_cetroid_in_middle_table,]
    mid_way_this_list <- round(nrow(middle_to_north))/2
    cells_around_middle_this_list <- unique(c( (mid_way_this_list - 3):mid_way_this_list,mid_way_this_list:(mid_way_this_list + 3)))
    middle_to_north <- middle_to_north[sample(cells_around_middle_this_list,3,replace=F),]

    middle_to_south <- middle[id_row_cetroid_in_middle_table:nrow(middle),]
    mid_way_this_list <- round(nrow(middle_to_south))/2
    cells_around_middle_this_list <- unique(c( (mid_way_this_list - 3):mid_way_this_list,mid_way_this_list:(mid_way_this_list + 3)))
    middle_to_south <- middle_to_south[sample(cells_around_middle_this_list,3,replace=F),]
    sample_table_to_plot <- rbind(most_north,middle_to_north,
                                  middle_to_south,most_south)
  }

  # if(all(rownames(sample_table_to_plot) !=  id_row_cetroid)){
  #   sample_table_to_plot <- rbind(sample_table_to_plot,
  #                                 centroid_population_distance_abundance[,1:ncol(sample_table_to_plot)])
  # }
  #
  if(any(rownames(most_south) ==  id_row_cetroid)){

    most_south <- most_south[-which(rownames(most_south) ==  id_row_cetroid),]

  }

  if(any(rownames(most_north) ==  id_row_cetroid)){

    most_north <- most_north[-which(rownames(most_north) ==  id_row_cetroid),]

  }
  # id_row_centroid_in_sampled <- which(rownames(sample_table_to_plot) == id_row_cetroid)
  #
  #
  # eucl_dist_pairwise_sampled <- distance(cbind(sample_table_to_plot$X,sample_table_to_plot$Y), method = "euclidean", use.row.names = TRUE)
  # distance_from_centroid <- eucl_dist_pairwise_sampled[id_row_centroid_in_sampled,]
  #
  # sample_table_to_plot <- cbind(sample_table_to_plot,distance_from_centroid=distance_from_centroid)
  #



  if(sample_it){
    most_north <- most_north
    most_south <- most_south
  } else {
    stop("this last part needs revision ")
    distance_from_centroid <- eucl_dist_pairwise[as.numeric(id_row_cetroid),]
    export_this_table <- cbind(table_to_plot,distance_from_centroid = distance_from_centroid)
  }

  return(list(centroid_population_distance_abundance = centroid_population_distance_abundance,
              most_north = most_north,
              most_south = most_south))
}


## function to create "edge.widthMap" object
edge.widthMap<-function(tree,x,res=100,...){
  if(!inherits(tree,"phylo"))
    stop("tree should be an object of class \"phylo\".")
  tree<-as.phylo(tree)
  #h<-max(nodeHeights(tree))
  #LL<-setNames(seq(0,h,length.out=res),1:res)
  #tree<-map.to.singleton(make.era.map(tree,LL))
  a<-fastAnc(tree,x)
  node.values<-c(x[tree$tip.label],a)
  edge.values<-apply(tree$edge,1,function(e,nv)
    mean(nv[e]),nv=node.values)
  edge.widths<-edge.values
  object<-list(tree=tree,edge.widths=edge.widths,
               node.values=node.values)
  class(object)<-"edge.widthMap"
  object
}

## plot method

adjust_species_with_landchange <- function(list_species_from_cpp,new_map,time_map_change){
  extinct_species_duelandchange <- 0
  population_reduction_duelandchange <- NULL
  list_species_after_landchange <- NULL
  for(i in 1:length(list_species_from_cpp)){
    focal <-  list_species_from_cpp[[i]]
    if(focal$RangeSize != 0){ # alive species

      surviving_cell <- NULL
      for(ii in 1:focal$RangeSize){
        cell_to_check <- c(focal$distribution[ii,2],focal$distribution[ii,1])
        if(new_map[cell_to_check[1],cell_to_check[2]] == 0){ # it is land, safe
          surviving_cell <- c(surviving_cell,ii)
        }

      }
      if(length(surviving_cell) < focal$RangeSize){
        cat("species ",focal$ID, "had a reduction of ",focal$RangeSize-length(surviving_cell), "\n")
        population_reduction_duelandchange <- c(population_reduction_duelandchange, focal$RangeSize-length(surviving_cell))
      }
      if(length(surviving_cell) == 0){
        cat("species ",focal$ID, "died with land change \n")
        extinct_species_duelandchange <- c(extinct_species_duelandchange + 1)
        focal$RangeSize <- 0
        focal$Death <- time_map_change
        focal$distribution <- "extinct"

      } else {
        focal$distribution <- focal$distribution[surviving_cell,]
        focal$RangeSize <- length(surviving_cell)
      }

      list_species_after_landchange [[i]] <- focal

    } else {
      list_species_after_landchange [[i]] <- focal
    }


  }
  return(list(list_species_after_landchange = list_species_after_landchange,
              population_reduction_duelandchange = population_reduction_duelandchange,
              extinct_species_duelandchange = extinct_species_duelandchange))
}




make_list_species_fromR_tocpp <- function(list_species_after_landchange){
  list_ready_forcpp <- list()

  ID <- NULL
  Parent <- NULL
  Birth <- NULL
  Death <- NULL
  TraitValue <- NULL
  RangeSize <- NULL

  distribution_x <- NULL
  distribution_y <- NULL

  for(i in 1:length(list_species_after_landchange)){
    focal <- list_species_after_landchange[[i]]
    #cat(i,"\n")
    ID <- c(ID,focal$ID)
    Parent <- c(Parent,focal$Parent)
    Birth <- c(Birth,focal$Birth)
    Death <- c(Death,focal$Death)
    TraitValue <- c(TraitValue,focal$TraitValue)
    RangeSize <- c(RangeSize,focal$RangeSize)

    if(focal$RangeSize > 0){
      if(focal$RangeSize == 1){
        distribution_x <- c(distribution_x, focal$distribution[1])
        distribution_y <- c(distribution_y, focal$distribution[2])
      } else {
        # for(ii in 1:focal$RangeSize){
        distribution_x <- c(distribution_x, focal$distribution[,1])
        distribution_y <- c(distribution_y, focal$distribution[,2])
        #  }
      }
    }
  }
  list_ready_forcpp <- list(ID=ID,
                            Parent=Parent,
                            Birth=Birth,
                            Death=Death,
                            TraitValue=TraitValue,
                            RangeSize=RangeSize,
                            distribution_x=distribution_x,
                            distribution_y=distribution_y)

  return(list_ready_forcpp)

}


transform_map_vector <- function(input_map){
  map_vector <- NULL
  for(i in 1:nrow(input_map)){
    for(ii in 1:ncol(input_map)){
      map_vector <- c(map_vector,input_map[i,ii])
    }
  }
  return(map_vector)
}





do_some_tests <- function(old_map,new_map,simulation_raw,list_species_after_landchange,list_ready_forcpp){

  change_map_thing <- adjust_species_with_landchange (list_species_from_cpp,new_map,time_map_change)

  alive_species_fromcpp <- length(which(simulation_raw$evolution_trait$Death == 0))
  alive_species_after_landchange <- 0
  for(i in 1:length(list_species_after_landchange)){

    if(list_species_after_landchange[[i]]$RangeSize > 0){
      alive_species_after_landchange <- alive_species_after_landchange + 1
    }
  }

  if((change_map_thing$extinct_species_duelandchange + alive_species_after_landchange) != alive_species_fromcpp){
    stop("problem with extinction")

  }

  all_populations_fromcpp <- sum(simulation_raw$evolution_trait$RangeSize)
  all_populations_fromready_forcpp <- sum(list_ready_forcpp$RangeSize)

  if((all_populations_fromcpp-sum(change_map_thing$population_reduction_duelandchange)) != all_populations_fromready_forcpp){
    stop("problem with range reduction")
  }

  land_reduction <- length(which(old_map==0)) - length(which(new_map==0))

  if(land_reduction > 0){
    if(any(change_map_thing$population_reduction_duelandchange) > land_reduction){
      stop("problem with range reduction 2")
    }
  }



  if(length(list_ready_forcpp$distribution_x) != sum(list_ready_forcpp$RangeSize)){
    stop("range size problem")
  }
  if (length(which(list_ready_forcpp$Death!=0)) != length(which(list_ready_forcpp$RangeSize==0))){
    stop("extinct problem")
  }
}
