#GRID SEARCH---------------
#' Draw Grid Search Simulation Agenda
#'
#' Generates a grid of treatment parameter combinations for evaluation using with grid search.
#' Parameter values are created based on the lower and upper bounds defined in the optimization object's
#' \code{searchSpace}. The generated grid is saved to file and stored in the object.
#'
#' @import dplyr
#' @importFrom rlist list.rbind
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}. Must include a defined \code{searchSpace}.
#' @param round_inputs A logical or integer vector indicating which treatment parameters should be rounded to whole numbers. If logical, \code{TRUE} values indicate rounding.
#' @param intersect_points Optional named vector. If provided, generates a sub-grid where one or more parameters are fixed at specified values (e.g. the identified optima). Parameters not mentioned in this vector are allowed to vary.
#' @param grid_resolution Integer. Number of points to generate per parameter range. Default is 50.
#'
#' @return The updated \code{opt_object} with a new \code{grid} element and status. The grid is also saved to \code{GridSearch.csv} in the project directory.
#' @export
drawSimulationGrid <- function(opt_object, round_inputs=c(), 
                               intersect_points=NULL, grid_resolution=50){
  
  # setup from opt object
  dose_inputs <- opt_object[['searchSpace']]
  
  # Create a list where each element is a sequence based on lb and ub
  seq_list <- lapply(seq_len(nrow(dose_inputs)), function(i) {
    seq(dose_inputs$lb[i], dose_inputs$ub[i], length.out = grid_resolution)
  })
  names(seq_list) <- rownames(dose_inputs)
  
  # expand grid to create simulation agenda
  par_values <- expand.grid(seq_list)
  
  # rounding based on input
  if(is.logical(round_inputs)) round_inputs <- which(round_inputs)
  for(i in round_inputs){
    par_values[,i] <- round(par_values[,i], 0)
  }
  par_values <- distinct(par_values)
  
  # select only intersects of the optima
  if(!is.null(intersect_points)){
    sim_agenda <- lapply(c(1:length(intersect_points)), function(xi){
      applied_par_values <- par_values
      applied_par_values[,names(intersect_points)[xi]] <- intersect_points[xi]
      applied_par_values <- distinct(applied_par_values)
      return(applied_par_values)
    }) %>% list.rbind() %>% distinct()
  }else{
    sim_agenda <- par_values
  }
  rownames(sim_agenda) <- NULL
  
  # append to object
  opt_object[['grid']] <- sim_agenda
  opt_object[['grid_run_status']] <- "Not run"
  
  # save grid
  write.csv(sim_agenda, paste0(opt_object[['WorkDirectory']], "GridSearch.csv"), row.names=F)
  
  # save object
  saveProject(opt_object)
  
  return(opt_object)
}

#' Run Grid Search
#'
#' Runs PKPD simulations over a grid of treatment parameter combinations stored in the optimization object.
#' Function \code{drawSimulationGrid} need to first be run before running this function to generate the grid table.
#' The results are evaluated using the objective function, and penalties are written to a result file.
#' Previously evaluated points are skipped unless \code{redo_run = TRUE}.
#' 
#' @import dplyr
#' @import rlist
#' @import parallel
#' @param opt_object An optimization task object that includes a valid simulation grid. The grid should be created using \code{drawSimulationGrid()}.
#' @param output_filename Optional character string specifying the name of the CSV file to save results. Default is \code{"GridSearchResult.csv"}.
#' @param n_cores Integer. Number of CPU cores to use for parallel simulation. If \code{n_cores > 1}, uses \code{mclapply()} for parallel execution. Default is 1 (sequential).
#' @param redo_run Logical. If TRUE, re-evaluates all grid points even if results already exist. Default is FALSE.
#'
#' @return The updated \code{opt_object} with updated \code{grid_run_status}. The results are appended to a CSV file in the working directory.
#' @export
gridSearchEvaluation <- function(opt_object, output_filename=NULL,
                                 n_cores=1, redo_run=F){
  if(is.null(opt_object[['grid']])){
    print("Grid missing! Please first draw the simulation grid using 'drawSimulationGrid'")
  }else{
    # setup
    sim_inputs <- opt_object[['simulationSettings']]
    n_subjects <- sim_inputs[['n_subjects']]
    grid_table <- opt_object[['grid']]
    
    # check previously generated values
    if(is.null(output_filename)){output_filename <- "GridSearchResult.csv"}
    output_file <- paste0(opt_object[['WorkDirectory']], output_filename)
    if(file.exists(output_file) & !redo_run){
      # read previously ran grid
      previous_grid <- read.csv(output_file, header=T) 
      
      if(penalty %in% colnames(previous_grid)){
        previous_grid <- dplyr::select(previous_grid, -penalty)
      }
      
      # summarize condition in grid
      previous_grid$condition <- apply(previous_grid, 1, function(x) paste(x, collapse='_'))
      
      # remove executed runs from the simulation agenda
      grid_table$condition <- apply(grid_table, 1, function(x) paste(x, collapse='_'))
      grid_table <- filter(grid_table, !condition %in% previous_grid$condition) %>% 
        dplyr::select(-condition)
    }else{
      # initiate output file
      empty_df <- data.frame(matrix(ncol = (ncol(grid_table)+
                                              nrow(opt_object[['treatmentObjectives']])), nrow = 0))
      obj_names <- apply(opt_object[['treatmentObjectives']], 1, function(x){
        paste("penalty", x['evaluation_column'], x['type'], sep="_")
      })
      colnames(empty_df) <- c(colnames(grid_table), obj_names)
      write.csv(empty_df, output_file, row.names = FALSE)
    }
    
    # status update
    opt_object[['grid_run_status']] <- "Run started . . . "
    saveProject(opt_object)
    
    # run through the simulation agenda
    if(n_cores<=1){
      nullCallBack <- lapply(c(1:nrow(grid_table)), function(xi){
        penalty_value <- objFn(unlist(grid_table[xi,]), colnames(grid_table),
                               sim_inputs,
                               model=opt_object[['pkpdModel']], 
                               evaluation_matrix=opt_object[['treatmentObjectives']],
                               pkpd_pars=opt_object[['modelParameters']], 
                               fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']],
                               population_samples=opt_object[['populationSample']], n_sub=n_subjects,
                               fun_postprocessing=opt_object[['fun_postProcessing']],
                               verbose_level=0, write_file=NULL,
                               dose_scaling = opt_object[['fun_doseScaling']])
        
        penalty_value <- c(unlist(grid_table[xi,]), penalty_value)
        
        write.table(data.frame(t(penalty_value)), output_file,
                    sep=",", row.names=F, col.names=F, append=T)
        return(NULL)
      })
    }else{
      nullCallBack <- mclapply(c(1:nrow(grid_table)), function(xi){
        penalty_value <- objFn(unlist(grid_table[xi,]), colnames(grid_table),
                               sim_inputs,
                               model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                               pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']],
                               population_samples=opt_object[['populationSample']], n_sub=n_subjects,
                               fun_postprocessing=opt_object[['fun_postProcessing']],
                               verbose_level=0,
                               write_file=NULL,
                               dose_scaling = opt_object[['fun_doseScaling']])
        penalty_value <- c(unlist(grid_table[xi,]), penalty_value)
        
        write.table(data.frame(t(penalty_value)), output_file,
                    sep=",", row.names=F, col.names=F, append=T)
        return(NULL)
      }, mc.cores=n_cores)
    }
  }
  
  # reporting
  opt_object[['grid_run_status']] <- "Completed"
  saveProject(opt_object)
  
  return(opt_object)
}

#PARETO---------------
#' Draw Pareto Optimization Grid
#'
#' Creates a grid of relative objective weights for multi-objective optimization. The grid is built
#' using historical optimization results to identify which objectives vary meaningfully, and constructs
#' combinations of normalized weights across those objectives. The resulting grid is saved to file and
#' stored in the optimization object.
#'
#' Note: \code{\link{runOptimization}} or \code{\link{runOptimization_Multiple}} must be run beforehand to generate the optimization history file required for this function.
#'
#' @import dplyr
#' @import rlist
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}. Must include previously run optimization history.
#' @param grid_values Numeric vector. Defines the set of relative weights to assign to each objective. Default is \code{c(0.1, 0.25, 0.5, 0.75, 1)}.
#' @param opt_history_name Character. The method name used in the previous optimization run (e.g., \code{"PSO-LBFGSB"}), which also serves as the prefix for the CSV file containing optimization history (e.g., \code{"PSO-LBFGSB_OptimizationHistory.csv"}).
#' @param nonzero_fraction_threshold Numeric. Threshold for determining whether an objective varies meaningfully across history. Objectives that are zero in more than \code{(1 - threshold)} of runs are fixed to the lowest priority. Default is \code{0.05}.
#'
#' @return The updated \code{opt_object} with a new \code{paretoGrid} element and status. The grid is also saved to \code{ParetoGrid.csv} in the project directory.
#' @export
drawParetoGrid <- function(opt_object, grid_values=c(0.1, 0.25, 0.5, 0.75, 1), 
                           opt_history_name="PSO-LBFGSB", nonzero_fraction_threshold=0.05){
  # collect objective list
  treatment_objectives <- opt_object[['treatmentObjectives']]
  
  # collect optimization history
  opt_history <- paste0(opt_object[['WorkDirectory']], opt_history_name, "_OptimizationHistory.csv")
  opt_history <- read.csv(opt_history, header=T)
  colnames(opt_history) <- apply(treatment_objectives, 1, function(x) paste(x['evaluation_column'], x['type'], sep="_"))
  
  # select varying objectives that is fully attained in more than 5% of the tested data points
  varying_parameters <- apply(opt_history, 2, function(x){1-mean(x==0)})
  varying_parameters <- varying_parameters>nonzero_fraction_threshold
  
  # create grid
  param_list <- setNames(rep(list(grid_values), length(varying_parameters)), names(varying_parameters))
  grid <- expand_grid(!!!param_list)
  
  # fix priority of non-varying parameter to lowest
  for(xi in c(1:length(varying_parameters))){
    if(!varying_parameters[xi]){
      grid[,xi] <- min(grid_values)
    }
  }
  grid <- distinct(grid)
  
  # normalize to maxima (relative priority)
  grid_normalized <- grid %>%
    as.matrix() %>% t() %>%
    apply(2, function(row) row / max(row)) %>%
    t() %>% as.data.frame()
  
  # Restore original column names
  colnames(grid_normalized) <- colnames(grid)
  
  # Remove duplicate rows
  grid <- distinct(grid_normalized)
  
  # append to optimization object
  opt_object[['paretoGrid']] <- grid
  opt_object[['pareto_run_status']] <- 'Not run'
  
  # write grid
  write.csv(grid, paste0(opt_object[['WorkDirectory']], "ParetoGrid.csv"), row.names=F)
  saveProject(opt_object)
  
  return(opt_object)
}

#' Perform Pareto Grid Optimization
#'
#' Runs multi-objective optimization using a grid of objective weight combinations defined in the \code{paretoGrid}.
#' For each weight setting, a separate optimization is performed, and the resulting parameters and penalty values
#' are saved to file. This approach enables Pareto front exploration across competing objectives.
#'
#' Note: \code{\link{drawParetoGrid}} must be run before executing this function.
#'
#' @import dplyr
#' @import rlist
#' @import parallel
#' 
#' @param opt_object An optimization task object that includes a valid \code{paretoGrid}. The grid should be created using \code{\link{drawParetoGrid}}.
#' @param opt_method Character. Name of the optimization algorithm to use for each Pareto point (e.g., \code{"PSO-LBFGSB"}, \code{"GA"}, \code{"BFGS"}). For list of available method, see \code{\link{runOptimization}}.
#' @param opt_control Named list of control parameters specific to the selected optimization method. See \code{\link{runOptimization}} for supported options.
#' @param output_filename Optional character string specifying the name of the CSV file to store the results. Default is \code{"ParetoResult.csv"}.
#' @param n_cores Integer. Number of CPU cores to use for parallel execution. If \code{n_cores > 1}, simulations are run using \code{mclapply()}. Default is 1 (sequential).
#' @param kurtosis_penalty Numeric. Penalty tuning parameter passed to the objective function. Represents the kurtosis parameter in a logistic function. Default is 50.
#' @param redo_run Logical. If TRUE, all Pareto grid points are re-evaluated even if results already exist. Default is FALSE.
#'
#' @return The updated \code{opt_object} with updated \code{pareto_run_status}. Results are appended to a CSV file in the working directory.
#' @export
paretoGridSearch <- function(opt_object, opt_method,
                             opt_control=list(),
                             output_filename=NULL, n_cores=1,
                             kurtosis_penalty=50, redo_run=F){
  # setup
  sim_inputs <- opt_object[['simulationSettings']]
  n_subjects <- sim_inputs[['n_subjects']]
  pareto_grid <- opt_object[['paretoGrid']]
  
  # check previously generated values
  if(is.null(output_filename)){output_filename <- "ParetoResult.csv"}
  output_file <- paste0(opt_object[['WorkDirectory']], output_filename)
  if(file.exists(output_file) & !redo_run){
    # read previously ran grid
    previous_grid <- read.csv(output_file, header=T)[,c(1:ncol(pareto_grid))]
    
    # summarize condition in grid
    previous_grid$condition <- apply(previous_grid, 1, function(x) paste(x, collapse='_'))
    
    # remove executed runs from the simulation agenda
    pareto_grid$condition <- apply(pareto_grid, 1, function(x) paste(x, collapse='_'))
    pareto_grid <- filter(pareto_grid, !condition %in% previous_grid$condition) %>% 
      dplyr::select(-condition)
  }else{
    # initiate output file
    empty_df <- data.frame(matrix(ncol = (ncol(pareto_grid)+
                                            nrow(opt_object[['searchSpace']])+
                                            nrow(opt_object[['treatmentObjectives']])), nrow = 0))
    colnames(empty_df) <- c(paste0("weights_", colnames(pareto_grid)), 
                            rownames(opt_object[['searchSpace']]), 
                            apply(opt_object[['treatmentObjectives']], 1, function(x) paste(x['evaluation_column'], x['type'], sep=".")))
    write.csv(empty_df, output_file, row.names = FALSE)
  }
  
  # status update
  opt_object[['pareto_run_status']] <- "Run started . . . "
  saveProject(opt_object)
  
  # run optimization per-grid object
  if(n_cores==1){
    # perform run
    nullCallBack <- lapply(c(1:nrow(pareto_grid)), function(xi){
      # adjust objective weights
      current_opt_task <- opt_object
      current_opt_task[['treatmentObjectives']]$weight <- unlist(pareto_grid[xi,])
      
      # run optimization
      opt_result <- runOptimization(current_opt_task, method_name=opt_method, 
                                    kurtosis = kurtosis_penalty, 
                                    verbose=2,
                                    opt_control = opt_control, save_object=F,
                                    return_optimization_result=T)
      
      # get parameter values
      opt_par <- opt_result[c(1:nrow(current_opt_task[['searchSpace']]))] %>% unlist()
      names(opt_par) <- rownames(current_opt_task[['searchSpace']])
      
      # get penalty values
      penalty_values <- objFn(opt_par, 
                              rownames(current_opt_task[['searchSpace']]),
                              sim_inputs,
                              model=current_opt_task[['pkpdModel']], 
                              evaluation_matrix=opt_object[['treatmentObjectives']], # USE UNWEIGHTED MATRIX
                              pkpd_pars=current_opt_task[['modelParameters']], 
                              fun_generalEventMatrix=current_opt_task[['fun_generateEventTable_general']],
                              population_samples=current_opt_task[['populationSample']], n_sub=n_subjects,
                              fun_postprocessing=current_opt_task[['fun_postProcessing']],
                              verbose_level=0, write_file=NULL,
                              dose_scaling = current_opt_task[['fun_doseScaling']])
      
      # append result
      row_write <- c(unlist(pareto_grid[xi,]), opt_par, penalty_values)
      
      # output
      write.table(data.frame(t(row_write)), output_file,
                  sep=",", row.names=F, col.names=F, append=T)
      return(NULL)
    })
  }else{
    # perform run
    library(parallel)
    nullCallBack <- mclapply(c(1:nrow(pareto_grid)), function(xi){
      # adjust objective weights
      current_opt_task <- opt_object
      current_opt_task[['treatmentObjectives']]$weight <- unlist(pareto_grid[xi,])
      
      # run optimization
      opt_result <- runOptimization(current_opt_task, method_name=opt_method, 
                                    kurtosis = kurtosis_penalty, 
                                    verbose=2,
                                    opt_control = opt_control, save_object=F,
                                    return_optimization_result=T)
      
      # get parameter values
      opt_par <- opt_result[c(1:nrow(current_opt_task[['searchSpace']]))] %>% unlist()
      names(opt_par) <- rownames(current_opt_task[['searchSpace']])
      
      # get penalty values
      penalty_values <- objFn(opt_par, 
                              rownames(current_opt_task[['searchSpace']]),
                              sim_inputs,
                              model=current_opt_task[['pkpdModel']], 
                              evaluation_matrix=opt_object[['treatmentObjectives']], # USE UNWEIGHTED MATRIX
                              pkpd_pars=current_opt_task[['modelParameters']], 
                              fun_generalEventMatrix=current_opt_task[['fun_generateEventTable_general']],
                              population_samples=current_opt_task[['populationSample']], n_sub=n_subjects,
                              fun_postprocessing=current_opt_task[['fun_postProcessing']],
                              verbose_level=0, write_file=NULL,
                              dose_scaling = current_opt_task[['fun_doseScaling']])
      
      # append result
      row_write <- c(unlist(pareto_grid[xi,]), opt_par, penalty_values)
      
      # output
      write.table(data.frame(t(row_write)), output_file,
                  sep=",", row.names=F, col.names=F, append=T)
      return(NULL)
    }, mc.cores=n_cores)
  }
  
  # reporting
  opt_object[['pareto_run_status']] <- "Completed"
  saveProject(opt_object)
  
  return(opt_object)
}