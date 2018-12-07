
##input_data 
path_java_class <- "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/java/MERS-CoV_Korea_ABM/MERS_Korea_IBM"
path_commons_math <- "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/lib/commons-math3-3.6.1.jar"
path_commons_lang <- "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/lib/commons-lang3-3.8.jar"

path_symptom_onset_data <- "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/ABM_Analysis/data/mers_symptom_onset_data.R"
path_Korea_shapefile <-"C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/Data_Analysis/data/KOR_adm"

run_java_ibm <- function( numruns = 1, 
                          stepsize = 0.2, 
                          stoptime = 60, 
                          beta = 0.35, 
                          frac_highrisk = 22/185, 
                          factor_highrisk = 7.9/0.1, 
                          shape_gamma_offspring = 0.2, 
                          hospital_search_radius = 30,
                          delay_move = 2,
                          outbreak_scenario = "2015",
                          under_vaccination_scenario = FALSE,
                          vaccination_scenario = "Distance",
                          vacc_coverage = 0.95,
                          vaccination_target_radius = 60,
                          vacc_eff = 0.7,
                          rel_vacc_eff_post_exposure = 0.5,
                          delay_vacc = 14,
                          threshold_case_vacc = 5,
                          threshold_day_vacc = 10,
                          random_seed = 0 ){
   

   Sys.setenv( JAVA_HOME = "C:/Program Files/Java/jdk-11.0.1/" )
   library( rJava )
   .jinit()
   .jaddClassPath( path_java_class )
   .jaddClassPath( path_commons_math )
   .jaddClassPath( path_commons_lang )
   # output containers
   daily_symptom_onset = matrix( NA, nrow=numruns, ncol=stoptime ) # store sim res of the ibm
   cumul_symptom_onset = matrix( NA, nrow=numruns, ncol=stoptime) # store sim res of the ibm
   var_mean_ratio = rep( NA, numruns ) # store sim res of the ibm
   num_hospital_infected = rep( NA, numruns ) # store sim res of the ibm
   cumul_vacc_dose = rep( NA, numruns )
   vacc_day_adjusted = rep( NA, numruns )
   
   daily_symptom_onset[ , 1] = 1
   # (R0 = 4.259 * beta *( (1-frac_highrisk) + factor_highrisk*frac_highrisk ))# average day before isolation during the initial stage
   # vacc_target_region_id = c( 0L, 1L, 6L, 8L, 11L, 13L, 15L, 16L) #as.integer(0:15)
   # vacc_target_region_id = as.integer(0:16)
   offspring_list = list()
   longitude_affected_hospital_list = list()
   latitude_affected_hospital_list = list()
   
   for( i in 1:numruns ){
      pars <- .jnew( "Parameters" )
      if( numruns > 1 ){
         random_seed = i
      }
      if( !is.null(random_seed) ){
         .jcall( pars, "V", "setRandomSeed", as.integer( random_seed ) )
         # cat( "random seed is being set to", random_seed, "\n" )
      }
      .jcall( pars, "V", "setStepSize", stepsize )
      .jcall( pars, "V", "setStopTime", stoptime )
      .jcall( pars, "V", "setRateTransmit", beta )
      .jcall( pars, "V", "setPropSeekingCareFromOtherHospitals", frac_highrisk )
      .jcall( pars, "V", "setFactorHighRiskTransmissibility", factor_highrisk )
      .jcall( pars, "V", "setShapeGammaOffspring", shape_gamma_offspring )
      .jcall( pars, "V", "setUnderVaccinationScenario", under_vaccination_scenario )
      .jcall( pars, "V", "setOutbreakScenario", outbreak_scenario )
      .jcall( pars, "V", "setVaccinationScenario", vaccination_scenario )
      .jcall( pars, "V", "setRadiusHospitalSearch", hospital_search_radius )
      .jcall( pars, "V", "setDayDelayBeforeMovingToAnotherHospital", delay_move )
      .jcall( pars, "V", "setVaccinationTargetRadius", vaccination_target_radius )
      .jcall( pars, "V", "setVaccCoverage", vacc_coverage )
      .jcall( pars, "V", "setVaccEfficacy", vacc_eff )
      .jcall( pars, "V", "setRelativeVaccEfficacyPostExposure", rel_vacc_eff_post_exposure )
      .jcall( pars, "V", "setMeanDelayVaccineInducedImmunity", delay_vacc )
      .jcall( pars, "V", "setThresholdNumberCaseForVaccinationInitiation", as.integer( threshold_case_vacc ) )
      .jcall( pars, "V", "setThresholdDayVaccinationInitiation", as.integer( threshold_day_vacc ) )
      
      model <- .jnew( "Model" ) 
      arr <- .jcall( model, "[[D", "runModel", pars )
      res <- sapply( arr, .jevalArray )
      
      cumul_symptom_onset[ i,  ] <- res[ 7, ]
      daily_symptom_onset[ i, 2:stoptime ] <- cumul_symptom_onset[ i, 2:stoptime ] - cumul_symptom_onset[ i, 1:(stoptime-1) ]
      var_mean_ratio[ i ] <- res[ 8, stoptime ]
      num_hospital_infected[ i ] <- res[ 9, stoptime ]
      offspring_list[[i]] <- .jcall( model, "[I", "getNumOffspringArray" )
      longitude_affected_hospital_list[[i]] <- .jcall( model, "[D", "getAffectedHospitalLongitudeArray" )
      latitude_affected_hospital_list[[i]] <- .jcall( model, "[D", "getAffectedHospitalLatitudeArray" )
      
      cumul_vacc_dose[ i ] <- .jcall( pars, "I", "getCumulVaccDose" )
      vacc_day_adjusted[ i ] <- .jcall( model, "I", "getDayVaccinationStartAdjusted", pars )
      
   }
   
   list <- list()
   list$cumul_symptom_onset <- cumul_symptom_onset
   list$daily_symptom_onset <- daily_symptom_onset
   list$var_mean_ratio <- var_mean_ratio
   list$num_hospital_infected <- num_hospital_infected
   list$offspring_list <- offspring_list
   list$longitude_affected_hospital_list <- longitude_affected_hospital_list
   list$latitude_affected_hospital_list <- latitude_affected_hospital_list
   list$cumul_vacc_dose <- cumul_vacc_dose
   list$vacc_day_adjusted <- vacc_day_adjusted
   
   return( list )
   
}


#########################################################################################################
## returns the only the final values of outbreak size, vaccine doses, ect
## over iteration
run_ibm_simple <- function( iter=1, ... ){
   # # stopifnot( length(res_var) > 0, length(res_var) > 0 )
   source( "data/mers_symptom_onset_data.R" )
   source( "util/mers_ibm_funcs.R" )  # run_java_ibm function is needed
   # cat( "beta=", beta,"\n" )
   list <- list() 
   for( i in 1:iter ){
      res <- run_java_ibm( ..., random_seed = i )
      list$cumul_symptom_onset[i] <- tail( res$cumul_symptom_onset[1,], 1 )
      list$daily_symptom_onset[i] <- tail( res$daily_symptom_onset[1,], 1 )
      list$var_mean_ratio[i] <- res$var_mean_ratio
      list$num_hospital_infected[i] <- res$num_hospital_infected
      list$offspring_list[i] <- res$offspring_list
      list$longitude_affected_hospital_list[i] <- res$longitude_affected_hospital_list
      list$latitude_affected_hospital_list[i] <- res$latitude_affected_hospital_list
      list$cumul_vacc_dose[i] <- res$cumul_vacc_dose
      list$vacc_day_adjusted[i] <- res$vacc_day_adjusted
   }
   
   filename <- paste0( "out/run_ibm_simple_",format(Sys.time(), "%Y%m%dT%H%M%S"), ".rds")
   saveRDS( list, filename )
   
   return( list )
}


