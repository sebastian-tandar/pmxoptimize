#COMPLEX - USER-EXPOSED ------------
# optimization
#' Run Treatment Optimization
#'
#' Executes model-based optimization using the specified algorithm to identify optimal treatment parameters. 
#' Supports various global and local optimization methods including PSO (particle swarm optimization), 
#' L-BFGS-B (limited memory BFGS with boundary), GA (genetic algorithm), SA (simulated annealing), DE (differential evolution), ABC (artificial bee colony), and others included in the 
#' from the \code{optimx} package. PSO-LBFGSB can be used to initiate the search with PSO and use L-BFGS-B as final optimizer. 
#' The function uses the information stored in the optimization task object to evaluate candidate solutions via PKPD simulations.
#'
#' @import dplyr
#' @import reshape2
#' @importFrom pso psoptim
#' @importFrom optimx optimx
#' @importFrom GA ga
#' @importFrom optimization optim_sa
#' @importFrom DEoptim DEoptim
#' @importFrom ABCoptim abc_optim
#' 
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}. Must include a valid population sample.
#' @param method_name Character. Name of the optimization algorithm to use. Supported methods include: \code{"PSO-LBFGSB"}, \code{"PSO"}, \code{"LBFGSB"}, \code{"GA"}, \code{"SA"}, \code{"ABC"}, \code{"DE"}, \code{"BFGS"}, \code{"Nelder-Mead"}, \code{"bobyqa"}, and \code{"nlminb"}.
#' @param kurtosis Numeric. A penalty tuning parameter passed to the objective function. Represents the kurtosis parameter in a logistic function. Default is 50.
#' @param verbose Integer (0–3). Level of verbosity during optimization. If set to 3, history of objective values will be written to file.
#' @param opt_control Named list of control parameters specific to the selected optimization method. The available options depend on the selected \code{method_name}:
#' \itemize{
#'   \item \strong{PSO-LBFGSB:}
#'     \describe{
#'       \item{\code{use_PSO} (logical)}{Whether to run PSO pre-fit before L-BFGS-B.}
#'       \item{\code{use_previous_PSO} (logical)}{Whether to use a previously saved PSO result.}
#'       \item{\code{pso_maxit} (integer)}{Maximum number of iterations for PSO.}
#'     }
#'   \item \strong{PSO:}
#'     \describe{
#'       \item{\code{pso_maxit} (integer)}{Maximum number of iterations for PSO.}
#'     }
#'   \item \strong{GA (Genetic Algorithm):}
#'     \describe{
#'       \item{\code{GA_maxit} (integer)}{Maximum number of generations for GA.}
#'       \item{\code{GA_parallel} (logical)}{Whether to run GA in parallel.}
#'     }
#'   \item \strong{ABC (Artificial Bee Colony):}
#'     \describe{
#'       \item{\code{ABC_maxit} (integer)}{Maximum number of iterations.}
#'       \item{\code{ABC_nFoodSource} (integer)}{Number of food sources.}
#'       \item{\code{ABC_foodLimit} (integer)}{Limit for abandoning food sources.}
#'     }
#'   \item \strong{Other methods (e.g., \code{LBFGSB}, \code{BFGS}, \code{DE}, \code{SA}, \code{Nelder-Mead}, \code{bobyqa}, \code{nlminb}):}
#'     \describe{
#'       \item{}{No additional control parameters currently implemented.}
#'     }
#' }
#' @param save_object Logical. If TRUE, the updated object will be saved in the working directory, rewriting the previous version. Default is TRUE.
#' @param return_optimization_result Logical. If TRUE, the optimization result object is returned directly. If FALSE, the updated \code{opt_object} is returned. Default is FALSE.
#'
#' @return Either the updated \code{opt_object} or the optimization result, depending on the value of \code{return_optimization_result}.
#' @export
runOptimization <- function(opt_object,
                            method_name = 'PSO-LBFGSB', 
                            kurtosis = 50,
                            verbose = 2, opt_control=list(), save_object=T,
                            return_optimization_result=F){
  if(is.null(opt_object[['populationSample']])){
    print("Population sample required!")
  }else{
    # setup loading
    sim_settings <- opt_object[['simulationSettings']]
    
    # load previous fit if prompted and if file existed
    if(method_name=='PSO-LBFGSB'){
      if(opt_control[['use_previous_PSO']]){
        if(file.exists(paste0(opt_object[['WorkDirectory']], "PSOprefit.Rdata"))){
          # if file exist, load previous prefit and use as init for BFGS
          load(paste0(opt_object[['WorkDirectory']], "PSOprefit.Rdata"))
        }else{
          # if file does not exist; override 'use_previous_prefit' request
          opt_control[['use_previous_PSO']]=F
        }
      }
    }
    
    # initiate output history track if verbose=3
    if(verbose==3){
      empty_df <- data.frame(matrix(ncol = nrow(opt_object[['treatmentObjectives']]), nrow = 0))
      colnames(empty_df) <- opt_object[['treatmentObjectives']]$evaluation_column
      write_file_name <- paste0(opt_object[['WorkDirectory']], method_name, "_OptimizationHistory.csv")
      write.csv(empty_df, write_file_name, row.names = FALSE)
    }else{
      write_file_name <- NULL
    }
    
    # method selection
    if(method_name=='PSO-LBFGSB'){
      # run PSO pre-fit
      # library(pso)
      print('### PSO ###')
      if(!opt_control[['use_previous_PSO']] & opt_control[['use_PSO']]){
        pso_prefit <- psoptim(opt_object[['searchSpace']]$start, fn=objFn,
                              par_names = rownames(opt_object[['searchSpace']]),
                              simulation_settings=sim_settings,
                              model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                              pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                              population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                              gamma = kurtosis,
                              verbose_level = verbose,
                              write_file = write_file_name,
                              dose_scaling = opt_object[['fun_doseScaling']],
                              fun_postprocessing=opt_object[['fun_postProcessing']],
                              lower=opt_object[['searchSpace']]$lb,
                              upper=opt_object[['searchSpace']]$ub,
                              control=list(maxit=opt_control[['pso_maxit']]))
        pso_prefitPars <- pso_prefit$par
        save(pso_prefit, file=paste0(opt_object[['WorkDirectory']], "PSOprefit.Rdata"))
      }else{
        if(!opt_control[['use_PSO']]){
          # get default starting value if use_PSO is not prompted
          pso_prefitPars <- opt_object[['searchSpace']]$start
        }else{
          # get previously fitted PSO parameter
          pso_prefitPars <- pso_prefit$par
        }
      }
      
      # run L-BFGS-B
      # library(optimx)
      print('### L-BFGS-B ###')
      optRes <- optimx(pso_prefitPars, fn=objFn,
                       par_names = rownames(opt_object[['searchSpace']]),
                       simulation_settings=sim_settings,
                       model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                       pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                       population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                       gamma = kurtosis,
                       verbose_level = verbose,
                       write_file = write_file_name,
                       dose_scaling = opt_object[['fun_doseScaling']],
                       fun_postprocessing=opt_object[['fun_postProcessing']],
                       lower=opt_object[['searchSpace']]$lb,
                       upper=opt_object[['searchSpace']]$ub,
                       method = 'L-BFGS-B')
      save(optRes, file=paste0(opt_object[['WorkDirectory']], "PSO_LBFGSB_mainfit.Rdata"))
    }else if(method_name=='PSO'){
      # library(pso)
      print('### PSO ###')
      optRes <- psoptim(opt_object[['searchSpace']]$start, fn=objFn,
                        par_names = rownames(opt_object[['searchSpace']]),
                        simulation_settings=sim_settings,
                        model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                        pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                        population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                        gamma = kurtosis,
                        verbose_level = verbose,
                        write_file = write_file_name,
                        dose_scaling = opt_object[['fun_doseScaling']],
                        fun_postprocessing=opt_object[['fun_postProcessing']],
                        lower=opt_object[['searchSpace']]$lb,
                        upper=opt_object[['searchSpace']]$ub,
                        control=list(maxit=opt_control[['pso_maxit']]))
      save(optRes, file=paste0(opt_object[['WorkDirectory']], method_name, "_mainfit.Rdata"))
    }else if(method_name=='GA'){
      # library(GA)
      print('### GA ###')
      optRes <- ga(type='real-valued', 
                   fitness = objFn,
                   par_names = rownames(opt_object[['searchSpace']]),
                   simulation_settings=sim_settings,
                   model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                   pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                   population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                   gamma = kurtosis,
                   verbose_level = verbose,
                   write_file = write_file_name,
                   dose_scaling = opt_object[['fun_doseScaling']],
                   fun_postprocessing=opt_object[['fun_postProcessing']],
                   maxiter=opt_control[['GA_maxit']],
                   lower=opt_object[['searchSpace']]$lb,
                   upper=opt_object[['searchSpace']]$ub,
                   parallel=opt_control[['GA_parallel']])
      save(optRes, file=paste0(opt_object[['WorkDirectory']], method_name, "_mainfit.Rdata"))
    }else if(method_name %in% c('SA', 'DE')){
      # methods requiring wrapper function
      # create wrapper function
      wrapperFunction <- function(pars){
        objFn(pars, 
              par_names = rownames(opt_object[['searchSpace']]),
              simulation_settings=sim_settings,
              model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
              pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
              population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
              gamma = kurtosis,
              verbose_level = verbose,
              write_file = write_file_name,
              dose_scaling = opt_object[['fun_doseScaling']],
              fun_postprocessing=opt_object[['fun_postProcessing']])
      }
      
      # run SA
      if(method_name=='SA'){
        # library(optimization)
        print('### SA ###')
        optRes <- optim_sa(fun=wrapperFunction,
                           start=opt_object[['searchSpace']]$start,
                           maximization = F, trace = F,
                           lower=opt_object[['searchSpace']]$lb,
                           upper=opt_object[['searchSpace']]$ub)
      }else{
        # library(DEoptim)
        print("### DE ###")
        optRes <- DEoptim(wrapperFunction,
                          lower=opt_object[['searchSpace']]$lb,
                          upper=opt_object[['searchSpace']]$ub) # currently no control function (add later!)
      }
    }else if(method_name=='ABC'){
      # library(ABCoptim)
      print("### ABC ###")
      optRes <- abc_optim(opt_object[['searchSpace']]$start, fn=objFn,
                          par_names = rownames(opt_object[['searchSpace']]),
                          simulation_settings=sim_settings,
                          model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                          pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                          population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                          gamma = kurtosis,
                          verbose_level = verbose,
                          write_file = write_file_name,
                          dose_scaling = opt_object[['fun_doseScaling']],
                          fun_postprocessing=opt_object[['fun_postProcessing']],
                          lb=opt_object[['searchSpace']]$lb,
                          ub=opt_object[['searchSpace']]$ub,
                          FoodNumber = opt_control[['ABC_nFoodSource']],
                          maxCycle = opt_control[['ABC_maxit']],
                          limit = opt_control[['ABC_foodLimit']])
    }else{
      # library(optimx)
      print(paste('###', method_name, '###', sep=" "))
      
      # rename LBFGSB
      use_method_name <- ifelse(method_name=='LBFGSB', "L-BFGS-B", method_name)
      
      # other methods by optimx - remove boundary when method is BFGS
      if(method_name != 'BFGS'){
        lower_bound <- opt_object[['searchSpace']]$lb
        upper_bound <- opt_object[['searchSpace']]$ub
      }else{
        lower_bound <- NULL
        upper_bound <- NULL
      }
      
      # run optimization
      optRes <- optimx(opt_object[['searchSpace']]$start, fn=objFn,
                       par_names = rownames(opt_object[['searchSpace']]),
                       simulation_settings=sim_settings,
                       model=opt_object[['pkpdModel']], evaluation_matrix=opt_object[['treatmentObjectives']],
                       pkpd_pars=opt_object[['modelParameters']], fun_generalEventMatrix=opt_object[['fun_generateEventTable_general']], 
                       population_samples=opt_object[['populationSample']], n_sub=sim_settings['n_subjects'],
                       gamma = kurtosis,
                       verbose_level = verbose,
                       write_file = write_file_name,
                       dose_scaling = opt_object[['fun_doseScaling']],
                       fun_postprocessing=opt_object[['fun_postProcessing']],
                       lower=lower_bound, upper=upper_bound, 
                       method = use_method_name)
    }
    
    # save output
    save(optRes, file=paste0(opt_object[['WorkDirectory']], method_name, "_mainfit.Rdata"))
    
    # add marker to object
    opt_object[['latest_optimization_run']] <- paste0(opt_object[['WorkDirectory']], method_name, "_mainfit.Rdata")
    
    # save progress
    if(save_object){
      saveProject(opt_object)
    }
    
    if(return_optimization_result){
      return(optRes)
    }else{
      return(opt_object)
    }
  }
}

#' Run Multiple Optimization Methods
#'
#' Executes multiple optimization runs sequentially or in parallel, based on a predefined agenda.
#' Each optimization method and its control settings are specified in the \code{optimization_agenda}.
#' Results are stored in the project object and saved to disk after all runs complete.
#'
#' @import parallel
#' 
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}. Must include a valid population sample.
#' @param optimization_agenda A list of optimization method configurations. Each element must be a named list with:
#' \describe{
#'   \item{\code{method}}{A character string specifying the optimization method (e.g., \code{"PSO-LBFGSB"}, \code{"GA"}, \code{"BFGS"}).}
#'   \item{\code{control}}{A named list of control parameters specific to the method. See \code{runOptimization()} for details.}
#' }
#' @param kurtosis Numeric. A penalty tuning parameter passed to the objective function. Represents the kurtosis parameter in a logistic function. Default is 50.
#' @param verbose Integer (0–3). Level of verbosity during optimization. Default is 3.
#' @param n_cores Integer. Number of CPU cores to use. If \code{n_cores > 1}, optimization is run in parallel using \code{mclapply()}. Default is 1.
#'
#' @return The updated \code{opt_object} after running all specified optimization methods.
#' @export
runOptimization_Multiple <- function(opt_object, optimization_agenda,
                                     kurtosis=50, verbose=3, n_cores=1){
  # run optimization
  if(n_cores==1){
    opt_objects <- lapply(optimization_agenda, function(x){
      runOptimization(opt_object, method_name=x[['method']],
                      kurtosis=kurtosis, verbose=verbose,
                      opt_control = x[['control']], save_object=F)
    })  
  }else{
    # library(parallel)
    opt_objects <- mclapply(optimization_agenda, function(x){
      runOptimization(opt_object, method_name=x[['method']],
                      kurtosis=kurtosis, verbose=verbose,
                      opt_control = x[['control']], save_object=F)
    }, mc.cores=n_cores)
  }
  
  # add marker to object
  opt_object[['latest_optimization_run']] <- paste(unlist(lapply(optimization_agenda, function(x) x[['method']])), collapse=" | ")
  
  # save object
  saveProject(opt_object)
  
  return(opt_object)
}

#' Run PKPD Simulation
#'
#' Executes a PKPD simulation for a given set of treatment parameters using the model and settings
#' stored in the optimization task object. The simulation includes event table generation,
#' optional dose scaling based on population characteristics, and optional post-processing of results.
#'
#' @import dplyr
#' @import rxode2
#' 
#' @param opt_object An optimization task object as created by \code{initiateOptimizationTask}. Must contain a valid population sample, model, and simulation settings.
#' @param treatment_parameters A named numeric vector of treatment parameter values. Names must match the row names of the \code{searchSpace} component in \code{opt_object}.
#'
#' @return A data frame containing the simulated PKPD results. If a post-processing function is defined in \code{opt_object}, it is applied before returning.
#' @export
runSimulation <- function(opt_object, treatment_parameters){
  # generate event matrix
  event_matrix <- opt_object[['fun_generateEventTable_general']](treatment_parameters, 
                                                                 rownames(opt_object[['searchSpace']]), 
                                                                 opt_object[['simulationSettings']]) %>% 
    eventMatrix_appendEtaDemographicsSample(opt_object[['populationSample']], 
                                            n_subjects = opt_object[['simulationSettings']]['n_subjects'],
                                            dose_scaling_function = opt_object[['fun_doseScaling']])
  
  # run simulation
  sim_res <- opt_object[['pkpdModel']]$solve(unlist(opt_object[['modelParameters']]), event_matrix) %>% data.frame()
  
  # apply post-processing function if given
  if(!is.null(opt_object[['fun_postProcessing']])){
    sim_res <- opt_object[['fun_postProcessing']](sim_res)
  }
  
  return(sim_res)
}

#NON-USER ---------------
# Optimization - units
activationFunction <- function(x, x_threshold, kurtosis = 50){
  1/(1 + exp(-kurtosis * (x/x_threshold - 1)))
}

#' @import dplyr
#' @import tidyr
evaluation_singleObjective <- function(sim_result, evaluation_column, time_evaluation, 
                                       threshold_max, threshold_min=NA,
                                       type = 'state', time_max = NA, 
                                       source_type = 'simulation',
                                       event_matrix = NA,
                                       gamma_value = 50,
                                       eval_format = 'continuous'){
  # use event matrix if prompted instead
  if(source_type=='event_matrix' & length(event_matrix)>1){
    sim_result <- event_matrix
    colnames(sim_result) <- tolower(colnames(sim_result))
    sim_result <- sim_result[,c('id', 'time', evaluation_column)] %>% 
      na.omit() %>% rename("sim.id"=id)
  }
  
  # select relevant columns
  sim_result <- sim_result[,c("sim.id", 'time', evaluation_column)] 
  colnames(sim_result)[colnames(sim_result)==evaluation_column] <- "evaluation_value"
  
  # subset relevant time period
  if(!is.na(time_max)){
    sim_result <- subset(sim_result, (time >= time_evaluation) & (time <= time_max))
  }else{
    if(!is.na(time_evaluation)){
      sim_result <- subset(sim_result, time == time_evaluation)
    }
  }
  
  # calculate AUC if necessary
  if(tolower(type)=="auc"){
    sim_result$idChange <- c(T, sim_result$sim.id[1:(nrow(sim_result)-1)]!=sim_result$sim.id[2:nrow(sim_result)])
    sim_result$fragment <- c(NA, 
                             (sim_result$evaluation_value[2:nrow(sim_result)] + 
                                sim_result$evaluation_value[1:(nrow(sim_result)-1)]) *
                               (sim_result$time[2:nrow(sim_result)] - 
                                  sim_result$time[1:(nrow(sim_result)-1)]) * 0.5)
    sim_result <- filter(sim_result, !idChange) %>% 
      dplyr::select(-evaluation_value) %>% group_by(sim.id) %>% 
      summarize(evaluation_value = sum(fragment), .groups='drop')
  }else if(tolower(type)=='maxima' | tolower(type)=='minima'){
    # search for maxima anyway
    sim_result <- group_by(sim_result, sim.id) %>% 
      summarize(evaluation_value = max(evaluation_value), .groups='drop')
  }else if(tolower(type)=='mean' | tolower(type)=='average'){
    sim_result <- group_by(sim_result, sim.id) %>% 
      summarize(evaluation_value = mean(evaluation_value), .groups='drop')
  }
  
  # define threshold
  if(!is.na(threshold_max) & !is.na(threshold_min)){
    # two-tail; convert to lower-tail test
    mean_threshold <- mean(c(threshold_max, threshold_min))
    value_threshold <- abs(threshold_max - mean_threshold)
    sim_result$evaluation_value <- abs(sim_result$evaluation_value - mean_threshold) # use absolute instead of square to keep reasonable scale during penalty calculation
    test_type <- 'upper_tail'
  }else{
    # value_threshold <- threshold_max
    if(is.na(threshold_max)){
      # is an upper-tailed test
      value_threshold <- threshold_min
      test_type <- 'lower_tail'
    }else{
      # is a lower-tailed test
      value_threshold <- threshold_max
      test_type <- 'upper_tail'
    }
  }
  
  # choose evaluation format (binary or continuous)
  if(eval_format=='binary'){
    # flip result and threshold value for lower-tailed test
    if(test_type=='lower_tail'){
      sim_result[,'evaluation_value'] <- -sim_result[,'evaluation_value']
      value_threshold <- -value_threshold
    }
    
    # perform evaluation
    crit_failed <- sim_result$evaluation_value > value_threshold
    
    # penalize NA values
    crit_failed[is.na(crit_failed)] <- 1
    
    # return mean penalty value
    return(mean(crit_failed))
  }else{
    # perform evaluation
    penalty <- activationFunction(sim_result$evaluation_value, value_threshold, kurtosis = gamma_value)
    
    # flip penalty value if test is a minima test
    if(test_type=='lower_tail'){penalty <- 1-penalty}
    
    # penalize NA values
    penalty[is.na(penalty)] <- 1
    
    # return mean penalty value
    return(mean(penalty))    
  }
}

#' @import dplyr
#' @import tidyr
evaluation_allObjectives <- function(simulation_result, evaluation_table, 
                                     event_matrix_eval,
                                     gamma_kurtosis = 50,
                                     evaluation_format = 'continuous'){
  
  # standardize input column caps
  colnames(simulation_result) <- tolower(colnames(simulation_result))
  evaluation_table$evaluation_column <- tolower(evaluation_table$evaluation_column)
  
  # id column name adjustment
  if("id" %in% colnames(simulation_result)){
    simulation_result <- rename(simulation_result, "sim.id"=id)
  }
  
  # calculate mean penalty value for each criteria
  penalty_values <- vapply(c(1:nrow(evaluation_table)), FUN.VALUE=1, function(xi){
    evaluation_singleObjective(simulation_result,
                               evaluation_column=evaluation_table$evaluation_column[xi],
                               time_evaluation=evaluation_table$time_evaluation[xi],
                               threshold_max=evaluation_table$threshold_max[xi],
                               threshold_min=evaluation_table$threshold_min[xi],
                               type=evaluation_table$type[xi],
                               time_max=evaluation_table$time_max[xi],
                               source=evaluation_table$source[xi],
                               event_matrix = event_matrix_eval,
                               gamma = gamma_kurtosis,
                               eval_format = evaluation_format)
  })
  
  # calculate weighted penalty value if evaluation is non-binary
  if(evaluation_format != 'binary'){
    weighted_penalty <- penalty_values * evaluation_table$weight
  }else{
    weighted_penalty <- penalty_values
  }
  
  # append names
  names(weighted_penalty) <- evaluation_table$evaluation_column
  
  return(weighted_penalty)
}

# loop event matrix operations
#' @import dplyr
#' @import tidyr
eventMatrix_appendEtaDemographicsSample <- function(event_matrix, demographics_eta_samples, 
                                                    dose_scaling_function=NULL, n_subjects=NA){
  # use less ID if prompted
  if(!is.na(n_subjects)){
    max_id_in_table <- max(demographics_eta_samples$id)
    if(n_subjects <= max_id_in_table){
      demographics_eta_samples <- subset(demographics_eta_samples, id <= n_subjects)
    }else{
      print("Number of subjects less than available in the demographics/eta sample table!")
      print(paste0("Max. number of IDs available on sample table applied (", max_id_in_table, ")"))
    }
  }
  
  # time-id matrix to expand time (regular definition by et() was slow)
  time_matrix <- expand_grid(id = unique(demographics_eta_samples$id),
                             time = unique(event_matrix$time))
  
  # combine expansion and append demographics & eta samples
  event_matrix <- dplyr::select(event_matrix, -id) %>% 
    left_join(time_matrix, by='time', relationship='many-to-many') %>% 
    left_join(demographics_eta_samples, by='id') %>% 
    arrange(time) %>% arrange(id)
  
  # check if dose needs to be scaled
  if(!is.null(dose_scaling_function)){
    if(is.function(dose_scaling_function)){
      event_matrix <- dose_scaling_function(event_matrix)
    }else{
      for(i in c(1:nrow(dose_scaling_function))){
        # get cmt number to scale
        cmt_number <- dose_scaling_function$cmt_number[i]
        
        # get scaling column name
        scaling_colname <- dose_scaling_function$scaler[i]
        
        # scale dose
        event_matrix$amt[event_matrix$cmt==cmt_number & !is.na(event_matrix$amt)] <- 
          event_matrix$amt[event_matrix$cmt==cmt_number & !is.na(event_matrix$amt)] *
          event_matrix[event_matrix$cmt==cmt_number & !is.na(event_matrix$amt),scaling_colname]
      }      
    }
  }
  
  return(event_matrix)
}

# Optimization - Complex
#' @import dplyr
#' @import tidyr
#' @import rxode2
objFn <- function(pars, par_names,
                  simulation_settings,
                  model, evaluation_matrix,
                  pkpd_pars, fun_generalEventMatrix, 
                  population_samples, n_sub, 
                  gamma = 50,
                  fun_postprocessing=NULL,
                  verbose_level=2,
                  write_file=NULL,
                  dose_scaling = NULL){
  # generate event matrix
  event_matrix <- fun_generalEventMatrix(pars, par_names, simulation_settings) %>% 
    eventMatrix_appendEtaDemographicsSample(population_samples, n_subjects = n_sub,
                                            dose_scaling_function = dose_scaling)
  
  # run simulation
  sim_res <- model$solve(pkpd_pars, event_matrix) %>% data.frame()
  
  # apply post-processing function if given
  if(!is.null(fun_postprocessing)){
    sim_res <- fun_postprocessing(sim_res)
  }
  
  # run evaluation
  penalty_value <- evaluation_allObjectives(sim_res, evaluation_matrix,
                                            gamma_kurtosis = gamma,
                                            event_matrix_eval = event_matrix) %>% unlist()
  sum_penalty <- sum(penalty_value)
  
  # reporting
  if(verbose_level>=2){
    # visual report
    names(penalty_value) <- paste(evaluation_matrix$evaluation_column, evaluation_matrix$type, sep="_")
    print(round(penalty_value, 3))
    
    # write to history file at verbose level = 3
    if(verbose_level == 3){
      write.table(data.frame(t(round(penalty_value, 3))), write_file,
                  sep=",", row.names=F, col.names=F, append=T)
    }
    
  }else if(verbose_level==1){
    print(sum_penalty)
  }
  
  # return sum of weighted penalty values
  if(verbose_level > 0){
    return(sum_penalty)
  }else{
    return(penalty_value)
  }
}