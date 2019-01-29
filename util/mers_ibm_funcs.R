java_class_path_setup <- function(){
  Sys.setenv( JAVA_HOME = "C:/Program Files/Java/jdk-11.0.1/" )
  library( rJava )
  .jinit()
  .jaddClassPath( "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/java/MERS-CoV_Korea_ABM/MERS_Korea_IBM" )
  .jaddClassPath( "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/lib/commons-math3-3.6.1.jar" )
  .jaddClassPath( "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/lib/commons-lang3-3.8.jar" )
}

##input_data
path_java_home <- "C:/Program Files/Java/jdk-11.0.1/"
path_java_class <- "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/java/MERS-CoV_Korea_ABM/MERS_Korea_IBM"
path_commons_math <- "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/lib/commons-math3-3.6.1.jar"
path_commons_lang <- "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/lib/commons-lang3-3.8.jar"
path_symptom_onset_data <- "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/ABM_Analysis/data/mers_symptom_onset_data.R"
path_Korea_shapefile <-"C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/Data_Analysis/data/KOR_adm"
path_util_func <- "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/ABM_Analysis/util/mers_ibm_funcs.R"


par_dataframe <- function( numruns = 1, 
          stepsize = 0.2, 
          stoptime = 60, 
          beta = 0.35,  
          delay_move = 2,
          frac_highrisk = 22/185, 
          factor_highrisk = 7.9/0.1, 
          shape_gamma_offspring = 0.2, 
          hospital_search_radius = 30,
          outbreak_scenario = "2015",
          under_vaccination_scenario = FALSE,
          vaccination_scenario = "Distance",
          vacc_coverage = 0.9,
          vaccination_target_radius = 30,
          vacc_eff = 0.7,
          rel_vacc_eff_post_exposure = 0.5,
          delay_vacc = 14,
          threshold_case_vacc = 5,
          threshold_day_vacc = 10,
          random_seed = 0 ){
  
df <- data.frame(
    par= c(  "numruns", 
  "stepsize", 
  "stoptime" , 
  "beta"  , 
  "frac_highrisk", 
  "factor_highrisk", 
  "shape_gamma_offspring", 
  "hospital_search_radius",
  "delay_move",
  "outbreak_scenario",
  "under_vaccination_scenario",
  "vaccination_scenario",
  "vacc_coverage",
  "vaccination_target_radius",
  "vacc_eff",
  "rel_vacc_eff_post_exposure",
  "delay_vacc",
  "threshold_case_vacc",
  "threshold_day_vacc",
  "random_seed" ), 
  val = c( numruns, 
           stepsize,
           stoptime, 
           beta, 
           frac_highrisk, 
           factor_highrisk, 
           shape_gamma_offspring, 
           hospital_search_radius,
           delay_move,
           outbreak_scenario,
           under_vaccination_scenario,
           vaccination_scenario,
           vacc_coverage,
           vaccination_target_radius,
           vacc_eff,
           rel_vacc_eff_post_exposure,
           delay_vacc,
           threshold_case_vacc,
           threshold_day_vacc,
           random_seed  ), stringsAsFactors = FALSE )
  
  return( df )
}


run_java_ibm <- function( stepsize = 0.2, 
                          stoptime = 60, 
                          beta = 0.35, 
                          delay_move = 2,
                          frac_highrisk = 22/185, 
                          factor_highrisk = 7.9/0.1, 
                          shape_gamma_offspring = 0.2, 
                          hospital_search_radius = 30,
                          meantime_isolation = c(4.259, 2.400, 0.500),
                          outbreak_scenario = "2015",
                          under_vaccination_scenario = FALSE,
                          vaccination_scenario = "Distance",
                          vacc_coverage = 0.9,
                          vaccination_target_radius = 30,
                          vacc_efficacy = 0.7,
                          rel_vacc_eff_post_exposure = 0.5,
                          delay_vacc = 14,
                          threshold_case_vacc = 5,
                          threshold_day_vacc = 10,
                          dur_vaccination = 10,
                          random_seed = 0 ){
   

   Sys.setenv( JAVA_HOME = "C:/Program Files/Java/jdk-11.0.1/" )
   library( rJava )
   .jinit()
   .jaddClassPath( "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/java/MERS-CoV_Korea_ABM/MERS_Korea_IBM" )
   .jaddClassPath( "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/lib/commons-math3-3.6.1.jar" )
   .jaddClassPath( "C:/Users/jonghoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/lib/commons-lang3-3.8.jar" )
   # output containers
   daily_symptom_onset = rep( NA, stoptime ) # store sim res of the ibm
   cumul_symptom_onset = rep( NA, stoptime) # store sim res of the ibm
   var_mean_ratio = NA # store sim res of the ibm
   num_hospital_infected = NA # store sim res of the ibm
   cumul_vacc_dose = NA
   vacc_day_adjusted = NA
   
   daily_symptom_onset[ 1 ] = 1
   # (R0 = 4.259 * beta *( (1-frac_highrisk) + factor_highrisk*frac_highrisk ))# average day before isolation during the initial stage
   # vacc_target_region_id = c( 0L, 1L, 6L, 8L, 11L, 13L, 15L, 16L) #as.integer(0:15)
   # vacc_target_region_id = as.integer(0:16)
   cumulative_case_list = list()
   offspring_list = list()
   longitude_affected_hospital_list = list()
   latitude_affected_hospital_list = list()
   
    pars <- .jnew( "Parameters" )
    .jcall( pars, "V", "setRandomSeed", as.integer( random_seed ) )
    .jcall( pars, "V", "setStepSize", stepsize )
    .jcall( pars, "V", "setStopTime", stoptime )
    .jcall( pars, "V", "setRateTransmit", beta )
    .jcall( pars, "V", "setMeanTimeToIsolation", meantime_isolation )
    .jcall( pars, "V", "setPropSeekingCareFromOtherHospitals", frac_highrisk )
    .jcall( pars, "V", "setFactorHighRiskTransmissibility", factor_highrisk )
    .jcall( pars, "V", "setShapeGammaOffspring", shape_gamma_offspring )
    .jcall( pars, "V", "setUnderVaccinationScenario", under_vaccination_scenario )
    .jcall( pars, "V", "setOutbreakScenario", as.character( outbreak_scenario ) )
    .jcall( pars, "V", "setVaccinationScenario", vaccination_scenario )
    .jcall( pars, "V", "setRadiusHospitalSearch", hospital_search_radius )
    .jcall( pars, "V", "setDayDelayBeforeMovingToAnotherHospital", delay_move )
    .jcall( pars, "V", "setVaccinationTargetRadius", vaccination_target_radius )
    .jcall( pars, "V", "setVaccCoverage", vacc_coverage )
    .jcall( pars, "V", "setVaccEfficacy", vacc_efficacy )
    .jcall( pars, "V", "setRelativeVaccEfficacyPostExposure", rel_vacc_eff_post_exposure )
    .jcall( pars, "V", "setMeanDelayVaccineInducedImmunity", delay_vacc )
    .jcall( pars, "V", "setDayNeededForVaccination", dur_vaccination )
    .jcall( pars, "V", "setThresholdNumberCaseForVaccinationInitiation", as.integer( threshold_case_vacc ) )
    .jcall( pars, "V", "setThresholdDayVaccinationInitiation", as.integer( threshold_day_vacc ) )
      
    model <- .jnew( "Model" ) 
    arr <- .jcall( model, "[[D", "runModel", pars )
    res <- sapply( arr, .jevalArray )
    
    cumul_symptom_onset <- res[ 7, ]
    daily_symptom_onset[ 2:stoptime ] <- cumul_symptom_onset[ 2:stoptime ] - cumul_symptom_onset[ 1:(stoptime-1) ]
    var_mean_ratio <- res[ 8, stoptime ]
    num_hospital_infected <- res[ 9, stoptime ]
    offspring_list[[1]] <- .jcall( model, "[I", "getNumOffspringArray" )
    longitude_affected_hospital_list[[1]] <- .jcall( model, "[D", "getAffectedHospitalLongitudeArray" )
    latitude_affected_hospital_list[[1]] <- .jcall( model, "[D", "getAffectedHospitalLatitudeArray" )
    
    cumul_vacc_dose <- .jcall( pars, "I", "getCumulVaccDose" )
    vacc_day_adjusted <- .jcall( model, "I", "getDayVaccinationStartAdjusted", pars )
    dur_outbreak <- .jcall( model, "D", "calcOutbreakDuration" )
      
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
    list$dur_outbreak <- dur_outbreak

    return( list )
}

#########################################################################################################
## returns the only the final values of outbreak size, vaccine doses, ect
## over iteration
run_ibm_simple <- function( iter=1, rng_seed_init=0, return_list=FALSE, ... ){
  list <- list()
  cumul_symptom_onset_list <- list()
  for( i in 1:iter ){
    rng_seed <- rng_seed_init + i
    res <- run_java_ibm( ..., random_seed = rng_seed )
    list$cumul_symptom_onset[ i ] <- tail( res$cumul_symptom_onset, 1 )
    # list$daily_symptom_onset[ i ] <- tail( res$daily_symptom_onset, 1 )
    list$var_mean_ratio[ i ] <- res$var_mean_ratio
    list$num_hospital_infected[ i ] <- res$num_hospital_infected
    list$cumul_vacc_dose[ i ] <- res$cumul_vacc_dose
    list$vacc_day_adjusted[ i ] <- res$vacc_day_adjusted
    list$dur_outbreak[ i ] <- res$dur_outbreak
    if( return_list ){
      cumul_symptom_onset_list[[ 1 ]] <- res$cumul_symptom_onset
      list$cumul_symptom_onset_list[ i ] <- cumul_symptom_onset_list
      list$offspring_list[ i ] <- res$offspring_list
      list$longitude_affected_hospital_list[ i ] <- res$longitude_affected_hospital_list
      list$latitude_affected_hospital_list[ i ] <- res$latitude_affected_hospital_list
    }
  }
  return( list )
}

#########################################################################################################
### run functions for approximate Bayesian computtion 

run_ibm_abc <- function( theta ){
  # # stopifnot( length(res_var) > 0, length(res_var) > 0 )
  source( "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/ABM_Analysis/data/mers_symptom_onset_data.R" )
  source( "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/ABM_Analysis/util/mers_ibm_funcs.R" )
  dat <- dat_symptom_onset_daily
  num_hospital_infected_observed = 14
  offspring_var_mean_ratio_observed = 51.769
  total_case_observed = 179
  
  stepsize = 0.2
  stoptime = 60
  beta = 0.35
  frac_highrisk = 22/185
  factor_highrisk = 7.9/0.1
  delay_move = 3
  shape_gamma_offspring = 0.2
  if( !is.na( theta[ 2 ] ) ) beta = theta[ 1 ]
  if( !is.na( theta[ 3 ] ) ) delay_move = theta[ 2 ]
  if( !is.na( theta[ 4 ] ) ) frac_highrisk = theta[ 3 ]
  if( !is.na( theta[ 5 ] ) ) factor_highrisk = theta[ 4 ] # inverse is taken to estimate the param between 0 and 1
  cat( "beta=", beta,"\n" )
  
  res <- run_java_ibm( stepsize = stepsize, 
                       beta = beta, 
                       frac_highrisk = frac_highrisk, 
                       factor_highrisk = factor_highrisk, 
                       shape_gamma_offspring = shape_gamma_offspring, 
                       delay_move = delay_move )
  
  sum_stat <- numeric()
  res_var <- c( "inc", "cum", "var_mean_ratio", "num_hospital" )
  # res_var <- c( "inc", "cum", "var_mean_ratio" )
  if( "inc" %in% res_var ){
    ssq_inc <- sqrt( sum( ( res$daily_symptom_onset - dat )^2 ) )
    sum_stat <- c( sum_stat, ssq_inc )
  }
  if( "cum" %in% res_var ){
    ssq_cumul <- sqrt( ( res$cumul_symptom_onset[stoptime] - total_case_observed )^2 )
    sum_stat <- c( sum_stat, ssq_cumul )
  }
  if( "var_mean_ratio" %in% res_var ){
    ssq_var_mean_ratio <- sqrt( ( res$var_mean_ratio - offspring_var_mean_ratio_observed )^2 )
    sum_stat <- c( sum_stat, ssq_var_mean_ratio )
  }
  if( "num_hospital" %in% res_var ){
    ssq_num_hospital_infected <- sqrt( ( res$num_hospital_infected - num_hospital_infected_observed )^2 )
    sum_stat <- c( sum_stat, ssq_num_hospital_infected )
  }
  return( sum_stat )
}
### run via a cluster of CPU's
run_ibm_abc_cls <- function( theta ){
  # # stopifnot( length(res_var) > 0, length(res_var) > 0 )
  source( "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/ABM_Analysis/data/mers_symptom_onset_data.R" )
  source( "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/ABM_Analysis/util/mers_ibm_funcs.R" )
  dat <- dat_symptom_onset_daily
  num_hospital_infected_observed = 14
  offspring_var_mean_ratio_observed = 51.769
  total_case_observed = 186
  
  stepsize = 0.2
  stoptime = 60
  beta = 0.35
  frac_highrisk = 22/185
  factor_highrisk = 7.9/0.1
  delay_move = 3
  shape_gamma_offspring = 0.2
  
  random_seed <- theta[ 1 ]
  
  if( !is.na( theta[ 2 ] ) ) beta <- theta[ 2 ]
  if( !is.na( theta[ 3 ] ) ) delay_move <- theta[ 3 ]
  if( !is.na( theta[ 4 ] ) ) frac_highrisk <- theta[ 4 ]
  if( !is.na( theta[ 5 ] ) ) factor_highrisk <- theta[ 5 ] 
  
  res <- run_java_ibm( stepsize = stepsize, 
                       beta = beta, 
                       frac_highrisk = frac_highrisk, 
                       factor_highrisk = factor_highrisk, 
                       delay_move = delay_move,
                       shape_gamma_offspring = shape_gamma_offspring,
                       random_seed = random_seed )
  
  # ssq_inc <- sqrt( sum( ( res$daily_symptom_onset - dat )^2 )  / sum( res$daily_symptom_onset^2 ) ) 
  # ssq_cumul <- sqrt( ( res$cumul_symptom_onset[stoptime] - total_case_observed )^2 / res$cumul_symptom_onset[stoptime]^2 )
  # ssq_var_mean_ratio <- sqrt( ( res$var_mean_ratio - offspring_var_mean_ratio_observed )^2 /  res$var_mean_ratio^2 )
  # ssq_num_hospital_infected <- sqrt( ( res$num_hospital_infected - num_hospital_infected_observed )^2 / res$num_hospital_infected^2 )
  
  ssq_cumul <- ( (res$cumul_symptom_onset[stoptime] - total_case_observed) / res$cumul_symptom_onset[stoptime] )^2 
  ssq_var_mean_ratio <- ( (res$var_mean_ratio - offspring_var_mean_ratio_observed)/res$var_mean_ratio )^2
  ssq_num_hospital_infected <- ( (res$num_hospital_infected - num_hospital_infected_observed )/res$num_hospital_infected )^2 
  
  # sum_stat <- c( ssq_inc, ssq_cumul, ssq_var_mean_ratio,  ssq_num_hospital_infected )
  sum_stat <- c( ssq_cumul, ssq_var_mean_ratio,  ssq_num_hospital_infected )
  # sum_stat <- c( ssq_cumul )
  return( sum_stat )
}

