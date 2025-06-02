#OBJECT OPERATION----------------
#' Initiate an Optimization Task
#'
#' This function creates and initializes a new optimization project object containing the PKPD model,
#' parameters, treatment objectives, and other components required for model-based optimization.
#' It sets up the working directory, stores all relevant inputs, and saves the object for later use.
#'
#' @param work_directory Path to the project folder. The directory will be created if it does not exist.
#' @param pkpd_model An \code{rxode2} model object representing the population PKPD model.
#' @param fixed_parameters A named numeric vector of fixed-effect parameter values. Names must match how the parameters are defined in the \code{pkpd_model}.
#' @param treatment_parameters A data frame with treatment parameters. Must contain three columns: \code{start}, \code{lb} (lower bound), and \code{ub} (upper bound). Row names should represent the parameter names.
#' @param evaluation_table A data frame defining treatment objectives. Must include columns such as \code{evaluation_column}, \code{time_evaluation}, \code{threshold_min}, \code{threshold_max}, \code{type}, \code{source}, and \code{weight}.
#' @param fun_generateGeneralEventTable A function that converts treatment parameters into an event matrix for a single individual. It must accept three arguments: treatment_parameters (numeric vector), treatment_parameters_names (character vector), and sim_settings (named vector).
#' @param simulation_settings Optional named vector specifying simulation parameters such as \code{simulation_tmax}, \code{simulation_dt}, and \code{n_subjects}. Can be \code{NULL} initially and updated later.
#' @param population_sample Optional data frame containing the population sample. Can also be generated using the \code{generatePopulationSample} function.
#' @param fun_dosescaling Optional function used to scale or adjust dosing based on individual characteristics from the population sample. Applied before simulation.
#' @param fun_postprocessing Optional function applied to the simulation results before calculating penalties or evaluation metrics.
#'
#' @return A list object containing all project components needed for optimization. The object is also saved to disk in the specified working directory.
#' @export
initiateOptimizationTask <- function(work_directory, pkpd_model, fixed_parameters,
                                     treatment_parameters, evaluation_table,
                                     fun_generateGeneralEventTable,
                                     simulation_settings = NULL, 
                                     population_sample = NULL, 
                                     fun_dosescaling = NULL,
                                     fun_postprocessing = NULL){
  # ensure path ends with backslash
  if(!grepl("[/\\\\]$", work_directory)){
    work_directory <- paste0(work_directory, "/")
  }
  
  # create work directory path
  dir.create(work_directory, showWarnings = FALSE)
  
  # create object
  opt_object <- list(WorkDirectory = work_directory,
                     pkpdModel = pkpd_model,
                     modelParameters = fixed_parameters,
                     searchSpace = treatment_parameters,
                     treatmentObjectives = evaluation_table,
                     populationSample = population_sample, # can be NULL
                     simulationSettings = simulation_settings, # can be NULL
                     
                     # functions
                     fun_generateEventTable_general = fun_generateGeneralEventTable, 
                     fun_doseScaling = fun_dosescaling, # can be NULL
                     fun_postProcessing = fun_postprocessing, # can be NULL
                     
                     # grid search
                     grid = NULL,
                     paretoGrid = NULL,
                     
                     # markers
                     latest_optimization_run = NULL,
                     grid_run_status = NA,
                     pareto_run_status = NA,
                     populationSample_file = NULL)
  
  # return
  saveProject(opt_object)
  return(opt_object)
}

# object modification
#' Modify Simulation Settings
#'
#' Updates the simulation settings in the optimization task object, including simulation duration,
#' time step, and number of simulated subjects. Optionally saves the updated object to file.
#'
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}.
#' @param simulation_tmax Numeric. The total simulation time.
#' @param simulation_dt Numeric. The time interval between simulation points.
#' @param n_subjects Integer. The number of individuals to simulate.
#' @param save_object Logical. If TRUE, the updated object will be saved in the working directory, rewriting the previous version. Default is TRUE.
#'
#' @return The updated \code{opt_object}.
#' @export
modifySimulationSettings <- function(opt_object, simulation_tmax, 
                                     simulation_dt, n_subjects, save_object=T){
  opt_object[['simulationSettings']] <- c(time_max = simulation_tmax, 
                                          dt_time = simulation_dt, 
                                          n_subjects = n_subjects) 
  
  # save changes
  if(save_object){
    saveProject(opt_object)
  }
  
  return(opt_object)
}

# save and load
#' Save Optimization Task Object
#'
#' Saves the optimization task object to disk as an \code{.Rdata} file, excluding large components
#' such as the population sample and simulation grids.
#'
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}.
#'
#' @return \code{NULL}. The object is saved to the \code{WorkDirectory} defined in the task object. Previously saved object will be rewritten.
#' @export
saveProject <- function(opt_object){
  # remove bulk information
  opt_object[['populationSample']] <- NULL
  opt_object[['grid']] <- NULL
  opt_object[['paretoGrid']] <- NULL
  
  # save object
  save(opt_object, file=paste0(opt_object[['WorkDirectory']], 'OptimizationTask.Rdata'))
  
  return(NULL)
}

#' Load Optimization Task Object
#'
#' Loads an optimization task object from a specified working directory. Optionally reloads
#' large components such as population samples and simulation results.
#'
#' @import rxode2
#' 
#' @param work_directory Path to the project folder containing the saved \code{OptimizationTask.Rdata} file.
#' @param load_all_contents Logical. If TRUE, large contents such as grid search table, pareto grid, and population samples are reloaded. Default is TRUE.
#'
#' @return The loaded optimization task object.
#' @export
loadProject <- function(work_directory, load_all_contents=T){
  # ensure path ends with backslash
  if(!grepl("[/\\\\]$", work_directory)){
    work_directory <- paste0(work_directory, "/")
  }
  
  # load object
  load(paste0(work_directory, 'OptimizationTask.Rdata'))
  
  # load heavier contents
  if(load_all_contents){
    opt_object <- reloadContents(opt_object)
  }
  
  # recompile model object
  opt_object[['pkpdModel']] <- rxode2(opt_object[['pkpdModel']]$model)
  
  return(opt_object)
}

#' Replicate Optimization Task
#'
#' Creates a copy of an existing optimization task in a new working directory, including
#' population samples and simulation grids if available. Useful for testing alternative scenarios.
#'
#' @param opt_object An optimization task object to replicate.
#' @param new_work_directory Path to the new directory where the replicated project will be stored.
#'
#' @return The replicated \code{opt_object} with the new working directory.
#' @export
replicateProject <- function(opt_object, new_work_directory){
  # reload all past content
  opt_object <- reloadContents(opt_object)
  
  # ensure path ends with backslash
  if(!grepl("[/\\\\]$", new_work_directory)){
    new_work_directory <- paste0(new_work_directory, "/")
  }
  
  # create work directory path
  dir.create(new_work_directory, showWarnings = FALSE)
  
  # append new work directory to object
  opt_object[['WorkDirectory']] <- new_work_directory
  opt_object[['populationSample_file']] <- paste0(new_work_directory, "PopulationSample.csv")
  
  # rewrite grid and population sample in the new directory
  if(!is.na(opt_object[['grid_run_status']])){
    write.csv(opt_object[['grid']], file=paste0(opt_object[['WorkDirectory']], "GridSearch.csv"), row.names=F)
  }
  
  if(!is.null(opt_object[['populationSample_file']])){
    write.csv(opt_object[['populationSample']], file=opt_object[['populationSample_file']], row.names=F)
  }
  
  # save object
  saveProject(opt_object)
  return(opt_object)
}

#NON-USER-----------
reloadContents <- function(opt_object){
  # load optimization grid if exist
  if(!is.na(opt_object[['grid_run_status']])){
    opt_object[['grid']] <- read.csv(paste0(opt_object[['WorkDirectory']], "GridSearch.csv"), header=T)
  }
  
  # load pareto grid if exist
  if(!is.na(opt_object[['pareto_run_status']])){
    opt_object[['paretoGrid']] <- read.csv(paste0(opt_object[['WorkDirectory']], "ParetoGrid.csv"), header=T)
  }
  
  # load population sample if exist
  if(!is.null(opt_object[['populationSample_file']])){
    opt_object[['populationSample']] <- read.csv(opt_object[['populationSample_file']], header=T)
  }
  
  return(opt_object)
}
