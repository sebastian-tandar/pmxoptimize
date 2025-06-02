#USER-EXPOSED-------
#' Run Simulation and Evaluate
#'
#' Runs simulation on the given treatment parameter and perform objective attainment evaluation.
#'
#' @import dplyr
#' @importFrom reshape2 dcast melt
#' @param pars A (named) vector representing the treatment design parameters to evaluate.  
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}.
#' Must contain the following elements:
#' \itemize{
#'   \item{\code{searchSpace}}{A data frame or matrix of variables to simulate.}
#'   \item{\code{simulationSettings}}{Named list with simulation configuration, including \code{n_subjects}.}
#'   \item{\code{populationSamples}}{Data frame or list of sampled individual ETAs and demographics.}
#'   \item{\code{fun_doseScaling}}{Function used to scale doses based on individual characteristics.}
#'   \item{\code{pkpdModel}}{A model object with a \code{$solve()} method (e.g., a compiled \code{rxode2} model).}
#'   \item{\code{modelParameters}}{Named numeric vector of fixed effect model parameters.}
#'   \item{\code{fun_postProcessing}}{(Optional) A function to transform or post-process simulation results.}
#'   \item{\code{treatmentObjectives}}{Specification of therapeutic targets or criteria for evaluation.}
#' }
#' @param expand_event_matrix Logical. If \code{TRUE}, expands the event matrix to include full observation timepoints. 
#' Default is \code{FALSE}.
#' @param quantile_values Numeric vector. Quantiles to summarize the simulation output (e.g., median, 95% CI). 
#' Default is \code{c(med = 0.5, lb = 0.025, ub = 0.975)}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{penalty_value}}{Named numeric vector of population binary failure rates computed against treatment objectives.}
#'   \item{\code{penalty_value_nonbinary}}{Named numeric vector of penalty contribution per-objective.}
#'   \item{\code{sim_res_summary}}{Simulation result.}
#' }
#' @export
runEvaluation <- function(pars, opt_object, expand_event_matrix=F, quantile_values=c(med=0.5, lb=0.025, ub=0.975)){
  # generate event matrix
  event_matrix <- opt_object[['fun_generateEventTable_general']](pars, rownames(opt_object[['searchSpace']]), opt_object[['simulationSettings']]) 
  
  # expand observation timepoints if prompted
  if(expand_event_matrix){
    event_matrix <- eventMatrix_appendFullObsTimes(event_matrix, opt_object[['simulationSettings']])
  }
  
  # append ETA & demographics sample
  event_matrix <- eventMatrix_appendEtaDemographicsSample(event_matrix,
                                                          opt_object[['populationSample']], 
                                                          n_subjects = opt_object[['simulationSettings']]['n_subjects'],
                                                          dose_scaling_function = opt_object[['fun_doseScaling']])
  
  # run simulation
  sim_res <- opt_object[['pkpdModel']]$solve(opt_object[['modelParameters']], event_matrix) %>% data.frame()
  
  # apply post-processing function if given
  if(!is.null(opt_object[['fun_postProcessing']])){
    sim_res <- opt_object[['fun_postProcessing']](sim_res)
  }
  
  # safekeeping for plotting
  sim_res_plotting <- melt(sim_res, id.vars=c("id", "time")) %>%
    dcast(time+variable~id, value.var='value')
  sim_res_summary <- apply(sim_res_plotting[,3:ncol(sim_res_plotting)], 1, function(x) quantile(x, quantile_values, na.rm=T)) %>%
    t() %>% data.frame() %>% mutate(time = sim_res_plotting$time, variable=sim_res_plotting$variable) 
  
  # adjust column names
  if(!is.null(names(quantile_values))){
    colnames(sim_res_summary)[1:length(quantile_values)] <- names(quantile_values)
  }
  
  # run evaluation
  penalty_value <- evaluation_allObjectives(sim_res, opt_object[['treatmentObjectives']],
                                            event_matrix_eval = event_matrix,
                                            evaluation_format='binary') %>% unlist()
  penalty_value_nonbinary <- evaluation_allObjectives(sim_res, opt_object[['treatmentObjectives']],
                                                      event_matrix_eval = event_matrix) %>% unlist() 
  
  return(list(failure_rate=penalty_value, penalty=penalty_value_nonbinary, simulation=sim_res_summary))
}

#NOT USER-EXPOSED------------
#' @import dplyr
eventMatrix_appendFullObsTimes <- function(event_matrix, sim_settings){
  # collect one representation of observation row
  obs_row <- event_matrix[which(is.na(event_matrix$cmt))[1],]
  
  # collect all observation timepoints that are not yet in the event matrix
  obs_times <- seq(0, sim_settings['time_max'], sim_settings['dt_time'])
  obs_times <- obs_times[!(obs_times %in% event_matrix$time)]
  
  # replicate observation rows through all observation timepoints and bind to the original event matrix
  obs_row <- obs_row[rep(1, length(obs_times)), ] %>% mutate(time = obs_times)
  event_matrix <- rbind.data.frame(event_matrix, obs_row) %>% arrange(time)
  
  return(event_matrix)
}