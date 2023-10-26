make_list_species_fromcpp_toR <- function(simulation_raw){
  list_species_from_cpp <- list()
  for(i in 1:length(simulation_raw$evolution_trait$Death)){
    
    
    
    ID <- simulation_raw$evolution_trait$ID[i]
    Parent <- simulation_raw$evolution_trait$Parent[i]
    Birth <- simulation_raw$evolution_trait$Birth[i]
    Death <- simulation_raw$evolution_trait$Death[i]
    elevational_origin <- simulation_raw$evolution_trait$elevational_origin[i]
    TraitValue <- simulation_raw$evolution_trait$TraitValue[i]
    RangeSize <- simulation_raw$evolution_trait$RangeSize[i]
    RangeHighlands <- simulation_raw$evolution_trait$RangeHighlands[i]
    RangeIntermediate1 <-simulation_raw$evolution_trait$RangeIntermediate1[i]
    RangeIntermediate2 <- simulation_raw$evolution_trait$RangeIntermediate2[i]
    RangeLowlands <- simulation_raw$evolution_trait$RangeLowlands[i]
    if(Death == 0){ #still alive
      for(ii in 1:RangeSize){
        distribution <- cbind(simulation_raw$Distribution[[i]][[1]],    simulation_raw$Distribution[[i]][[2]])
      temperature_evolved <- simulation_raw$temperature_evolved[[i]][[1]]
        }
      
    } else {
      distribution <- "extinct"
      temperature_evolved <- "extinct"
    }
    focal_species <- list(ID=ID,
                          Parent=Parent,
                          Birth=Birth,
                          Death=Death,
                          elevational_origin = elevational_origin,
                          TraitValue=TraitValue,
                          RangeSize=RangeSize,
                          RangeHighlands =RangeHighlands,
                          RangeIntermediate1= RangeIntermediate1,
                          RangeIntermediate2= RangeIntermediate2,
                          RangeLowlands = RangeLowlands,
                          distribution=distribution,
                          temperature_evolved = temperature_evolved)
    
    list_species_from_cpp[[i]] <- focal_species
    
  }
  return(list_species_from_cpp)
}



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
