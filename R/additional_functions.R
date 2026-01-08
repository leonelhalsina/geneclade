#' Computes population genetic metrics using output from run_geneclade_simulation()
#' @title Population genetic metrics
#' @param this_species_pop_info Table that includes, for a number of species, the xy coordinates of each population as well as their allelic frequencies.
#' @param pairwise Boolean, should the analysis be done in a pairwise fashion?
#' @param locus Select the locus to analysis "selection", "neutral", or "both"
#' @return A list to retrieve Weir and Cockrham Fst and Hudson's Fst.
#' @examples
#'library(geneclade)
#'library(DDD)
#'library(hierfstat)
#'# load a table which is one of the three elements retrived by run_geneclade_simulation().
#'tablespecies_overtime <- get("tablespecies_overtime")
#'this_time_slice <- tablespecies_overtime$`time: 9.5`
#'# Take species with ID = 1
#'this_species <- this_time_slice[which(this_time_slice$Sp_ID == 1),]
#'# remove columns with uneccesary info (xy coordinates for instance)
#'this_species_pop_info <- this_species[,-c(1,2,3,4)]

#'computed_genenetic_metrics <- compute_PopulationsSpecificFST(this_species_pop_info,
#'                                                             pairwise = FALSE,
#'                                                             locus="both")

#'# Hudson's Fst
#'mean(as.numeric(computed_genenetic_metrics$computed_betas$betaiovl),na.rm=TRUE)

#'# Weir and Cockrham Fst
#'computed_genenetic_metrics$wc_fst$FST

#' @export

compute_PopulationsSpecificFST <- function(this_species_pop_info,pairwise = FALSE, locus="both"){ # locus="both" or "neutral" or "selection"
  locus_selection <- NULL
  locus_neutral <- NULL
  population_id <- NULL
  this_species_pop_info <- as.data.frame(this_species_pop_info)
  for(i in 1:nrow(this_species_pop_info)){
    thispopulation <- this_species_pop_info[i,]
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

  if(locus == "both"){
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
    try(wc_fst <-  hierfstat::wc(table_for_bfts,diploid=F),
        ,silent=TRUE)
    if(pairwise){
      computed_betas <- hierfstat::pairwise.betas(table_for_bfts,diploid=F)
    } else {
      computed_betas <- hierfstat::betas(table_for_bfts,nboot=100,diploid=F)
    }

  } else {
    computed_betas$betaiovl <- NA
    wc_fst <- NA
  }



  return(list(computed_betas = computed_betas,
              wc_fst = wc_fst))
}

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
    phyloTree <- DDD::L2phylo(Ltable_to_do,dropextinct = T)
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



transform_map_vector <- function(input_map){
  map_vector <- NULL
  for(i in 1:nrow(input_map)){
    for(ii in 1:ncol(input_map)){
      map_vector <- c(map_vector,input_map[i,ii])
    }
  }
  return(map_vector)
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


