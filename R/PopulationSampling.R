#COMPLEX - USER-EXPOSED-------------
#' Generate a Population Sample
#'
#' This function generates a population sample by sampling inter-individual variability (IIV) parameters
#' and patient demographics. Continuous and categorical demographic values can be sampled from
#' provided ranges or distributions. The resulting sample is stored in the given optimization object,
#' written to file, and the project is saved.
#'
#' @import dplyr
#' 
#' @param opt_object mobTD object (a list object containing project configuration).
#' @param n_subject_population Integer. The number of individuals to include in the population sample.
#' @param random_effect_parameters A named numeric vector or list specifying variability parameters. The names must match the ETA identifiers used in the PKPD model included in the opt_object. If a parameter is provided as a percentage (%IIV), it will be internally converted to ETA scale for simulation. Additionally, any names containing 'iiv' will be automatically renamed to 'eta' to align with the model syntax.
#' @param random_effect_covariance Optional matrix. The covariance matrix of the random effects. If NULL, diagonal structure is assumed.
#' @param patient_demographics_continuous Optional named list. Each element defines a continuous demographic. If demographics_from_range = TRUE, each element should be a numeric vector of three values: median, lower bound, and upper bound. If demographics_from_range = FALSE, each element should be a numeric vector of two values: mean and standard deviation.
#' @param patient_demographics_categorical Optional named list. Each element defines a categorical demographic with level probabilities.
#' @param IIV_from_OMEGA Logical. If TRUE, values will be transformed back to OMEGA before resampling. Default is FALSE.
#' @param demographics_from_range Logical. If TRUE, continuous demographics will be sampled from specified ranges instead of empirical distributions. Default is FALSE.
#' @param force_demographics_inrange Logical. If TRUE and `demographics_from_range` is also TRUE, resampling will enforce that values stay within specified ranges. Default is FALSE.
#'
#' @return The updated \code{opt_object} list with the generated population sample. Population sample is also exported in the working directory.
#' @export
generatePopulationSample <- function(opt_object, n_subject_population,
                                     random_effect_parameters,
                                     random_effect_covariance = NULL,
                                     patient_demographics_continuous=NULL,
                                     patient_demographics_categorical = NULL,
                                     IIV_from_OMEGA = F,
                                     demographics_from_range = F,
                                     force_demographics_inrange = F){
  # A | get ETA samples
  populationSample <- sampleETAs(random_effect_parameters, n_subject_population, 
                                 random_effect_covariance,
                                 from_iiv = !IIV_from_OMEGA)
  
  # B | Get continuous patient demographics sample
  if(!is.null(patient_demographics_continuous)){
    demographicsSample <- sampleDemographics_continuous(patient_demographics_continuous, n_subject_population, 
                                                        from_range=demographics_from_range) %>% 
      mutate(id = as.numeric(id))
    
    if(demographics_from_range & force_demographics_inrange){
      # resample if any individuals are not within the intervals
      index_resample <- getResampleIndex(demographicsSample, patient_demographics_continuous)
      while(length(index_resample)>0){
        resampled_ids <- sampleDemographics_continuous(patient_demographics_continuous, 
                                                       length(index_resample), 
                                                       from_range=demographics_from_range) %>% 
          mutate(id = index_resample)
        demographicsSample <- subset(demographicsSample, !id %in% index_resample) %>% 
          rbind.data.frame(resampled_ids) %>% arrange(id)
        index_resample <- getResampleIndex(demographicsSample, patient_demographics_continuous)
      }
    }
  }
  
  # C | Get categorical patient demographics sample
  if(!is.null(patient_demographics_categorical)){
    demographicsSample_categorical <- lapply(patient_demographics_categorical, sampleDemographics_categorical, n_subject_population)
    if(length(demographicsSample_categorical)>1){
      demographicsSample_categorical <- reduce(demographicsSample_categorical, left_join, by = "id") # combine to one dataframe
    }else{
      demographicsSample_categorical <- demographicsSample_categorical[[1]]
    }
  }
  
  # D | Combine ETA and patient demographics samples
  if(!is.null(patient_demographics_continuous)){
    populationSample <- left_join(populationSample, demographicsSample, by='id') 
  }
  if(!is.null(patient_demographics_categorical)){
    populationSample <- left_join(populationSample, demographicsSample_categorical, by='id') 
  }
  
  # E | Move ID to first column
  id_column_index <- which(colnames(populationSample)=='id')
  oth_column_index <- which(colnames(populationSample)!='id')
  populationSample <- populationSample[,c(id_column_index, oth_column_index)]
  
  # bind to opt object
  opt_object[['populationSample']] <- populationSample
  opt_object[['populationSample_file']] <- paste0(opt_object[['WorkDirectory']], 'PopulationSample.csv')
  
  # write population sample
  write.csv(populationSample, opt_object[['populationSample_file']],
            row.names=F)
  
  # save progress
  saveProject(opt_object)
  
  return(opt_object)
}

#NON-USER-----------
# population sampling - supporting functions
getResampleIndex <- function(demographics_sample, patient_demographics){
  pass_metrics <- rep(T, nrow(demographics_sample))
  for(i in c(1:length(patient_demographics))){
    pass_metrics <- pass_metrics & 
      (demographics_sample[,names(patient_demographics)[i]] <= max(patient_demographics[[i]])) &
      (demographics_sample[,names(patient_demographics)[i]] > min(patient_demographics[[i]]))
  }
  ids_resample <- demographics_sample$id[!pass_metrics]
  return(ids_resample)
}
untransformIIVtoETA <- function(N){log((N)^2+1,base=exp(1))} # CHECK #

# population sampling - units
#' @import dplyr
#' @importFrom reshape2 melt dcast
#' @importFrom MASS mvrnorm
sampleETAs <- function(iiv_pars, n_subjects, covariances=NULL, from_iiv=T){
  # convert IIV to omega matrix
  if(from_iiv){
    eta_pars <- untransformIIVtoETA(iiv_pars)
  }else{
    eta_pars <- iiv_pars
  }
  omega_matrix <- diag(eta_pars)
  rownames(omega_matrix) <- gsub("iiv", "eta", names(iiv_pars))
  colnames(omega_matrix) <- rownames(omega_matrix)
  
  # parse location of relevant covariances if given
  if(!is.null(covariances)){
    for(i in c(1:length(covariances))){
      # get applied parameter index
      applied_parameters <- strsplit(names(covariances)[i], split="[.]")[[1]][2:3]
      applied_parameters <- paste0("eta.", applied_parameters)
      applied_parameters <- sapply(applied_parameters, function(x) which(rownames(omega_matrix)==x))
      
      # append covariance
      omega_matrix[applied_parameters[1],applied_parameters[2]] <- covariances[i]
      omega_matrix[applied_parameters[2],applied_parameters[1]] <- covariances[i]
    }
  }
  
  # get random effect samples
  eta_samples <- mvrnorm(n=n_subjects, mu=rep(0, ncol(omega_matrix)), Sigma=omega_matrix) %>% 
    data.frame() %>% mutate(id = c(1:n_subjects))
  
  return(eta_samples)
} # CHECK #

#' @import dplyr
#' @import rlist
#' @import MASS
sampleDemographics_continuous <- function(patient_demographics, n_subjects, from_range=T){
  # estimate SD from range
  # method from: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-5-13
  if(from_range){
    patient_demographics <- lapply(patient_demographics, function(x){
      sd_value <- (max(x) - min(x))/4
      return(c(x[1], sd_value))
    })
  }
  
  # create patient demographics
  patient_demographics_matrix <- lapply(patient_demographics, function(x){
    rnorm(n_subjects, mean=x[1], sd=x[2])
  }) %>% list.rbind() %>% t() %>% data.frame()
  colnames(patient_demographics_matrix) <- names(patient_demographics)
  
  # append ID
  patient_demographics_matrix <- mutate(patient_demographics_matrix, id = c(1:n_subjects))
  return(patient_demographics_matrix)
}

#' @import dplyr
sampleDemographics_categorical <- function(dist_probabilities, n_subjects){
  # get probability thresholds for each range
  probability_thresholds <- dist_probabilities / sum(dist_probabilities)
  probability_ranges <- vapply(c(1:length(probability_thresholds)), FUN.VALUE=1, function(xi) sum(probability_thresholds[1:xi]))
  names(probability_ranges) <- names(probability_thresholds)
  
  # sample probability
  probability_diceroll <- runif(n_subjects, min=0, max=1)
  
  # create occurrence matrix from sampled probability
  occ_matrix <- matrix(0, nrow=n_subjects, ncol=length(probability_ranges))
  for(i in c(1:length(probability_ranges))){
    lower_bound <- ifelse(i==1, 0, probability_ranges[i-1])
    occ_matrix[,i] <- (probability_diceroll > lower_bound) & 
      (probability_diceroll <= probability_ranges[i])
  }
  
  # append ID and adjust column names
  occ_matrix <- data.frame(occ_matrix) %>% 
    mutate(id=c(1:n_subjects))
  colnames(occ_matrix)[1:(ncol(occ_matrix)-1)] <- names(probability_ranges)
  
  return(occ_matrix)
}