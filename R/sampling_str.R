#' @title A function to create the sample sets for calibration/validation in
#' cross-validation.
#' @description for internal use only! This is stratified sampling based on the
#' values of a continuous response variable (y). If group is provided, the
#' sampling is done based on the groups and the average of y per group. This
#' function is used to create groups forleave-group-out cross-validations (or
#' leave-group-of-groups-out cross-validation if group argument is provided).
#' @param y a matrix of one column with the response variable.
#' @param p the percentage of samples (or groups if group argument is used) to
#' retain in the validation_indices set
#' @param number the number of sample groups to be crated
#' @param group the labels for each sample in \code{y} indicating the group each
#' observation belongs to.
#' @param replacement A logical indicating sample replacements for the 
#' calibration set are required.
#' @return a list with two matrices (\code{validation_indices} and 
#' \code{calibration_indices}) giving the indices of the observations in each 
#' column. The number of columns represents the number of sampling repetitions.
#' @keywords internal

sample_str <- function(y, p, number, group = NULL, replacement = FALSE) {
  
  ## If the percentage of samples to build the hold_in subset is below 50% of 
  ## the total number of samples, the selection is based on the number of samples 
  ## to retain.
  ## On the other hand, if the percentage of samples to build the hold_in subset 
  ## is equal or above 50% of the total number of samples, the selection is based 
  ## on the number of samples to exclude. 
  if (p < 0.5) {
    p_to_sample <- p
    do_sampling_for <- "calibration"
  } else {
    p_to_sample <- 1 - p
    do_sampling_for <- "validation"
  }
  
  if (is.null(group)) {
    nv <- floor(p_to_sample * nrow(y))
    nv <- ifelse(nv < 1, 1, nv)
    
    if (p >= 0.5) {
      n_val <- nv
      if (replacement) {
        n_cal <- nrow(y)
      } else {
        n_cal <- nrow(y) - nv
      }
    } else {
      n_val <- nrow(y) - nv
      if (replacement) {
        n_cal <- nrow(y)
      } else {
        n_cal <- nrow(y) - n_val
      }
    }
    
    strata_category <- get_sample_strata(y = y, 
                                         n = nv)
    
    calibration_indices <- matrix(NA, n_cal, number)
    validation_indices <- matrix(NA, n_val, number)
    
    colnames(calibration_indices) <- colnames(validation_indices) <- paste0(
      "Resample_",
      seq(1:number)
    )
    rownames(calibration_indices) <- paste0("index_", seq(1:n_cal))
    rownames(validation_indices) <- paste0("index_", seq(1:nrow(validation_indices)))
    
    indcs <- 1:nrow(y)
    for (jj in 1:number) {
      strata_samples <- get_samples_from_strata(original_order = strata_category$original_order, 
                                                strata = strata_category$strata, 
                                                replacement = replacement,
                                                sampling_for = do_sampling_for)
      
      calibration_indices[, jj] <- strata_samples$calibration
      validation_indices[, jj] <- strata_samples$validation
    }
    
  } else {
    
    y_groups <- data.table(
      y = y,
      group = factor(group),
      original_order = 1:length(y)
    )
    gr_levels <- levels(y_groups$group)
    n_levels <- nlevels(y_groups$group)
    # Compute means per group and make the stratified sampling based on that
    # vector of means
    aggregated_y <- tapply(
      X = y_groups$y,
      INDEX = y_groups$group,
      FUN = mean
    )
    
    nv <- floor(p_to_sample * n_levels)
    nv <- ifelse(nv < 1, 1, nv)
    
    if (p >= 0.5) {
      n_val <- nv
      if (replacement) {
        n_cal <- nrow(y)
      } else {
        n_cal <- nrow(y) - nv
      }
    } else {
      n_val <- nrow(y) - nv
      if (replacement) {
        n_cal <- nrow(y)
      } else {
        n_cal <- nrow(y) - n_val
      }
    }
    
    strata_category <- get_sample_strata(y = aggregated_y, 
                                         n = nv)
    strata_category$strata <-  paste0("l_", strata_category$strata)
    
    calibration_indices <- matrix(0, nv, number)
    colnames(calibration_indices) <- paste0("Resample_", seq(1:number))
    rownames(calibration_indices) <- paste0("index_", seq(1:nv))
    
    calibration_groups_indices <- validation_groups_indices <- NULL
    for (jj in 1:number) {
      strata_samples <- get_samples_from_strata(original_order = strata_category$original_order, 
                                                strata = strata_category$strata, 
                                                replacement = replacement,
                                                sampling_for = do_sampling_for)
      
      if (resampling) {
        sel_sample_indices <- lapply(1:length(strata_samples$calibration), 
                                     FUN = function(ith, full_groups, cal_groups) {
                                       ## this equivalent to extract from 
                                       ## y_groups$original_order
                                       which(full_groups %in% cal_groups[ith])
                                     }, 
                                     full_groups = y_groups$group,
                                     cal_groups = gr_levels[strata_samples$calibration])
        sel_sample_indices <- do.call("c", sel_sample_indices)
      } else {
        sel_sample_indices <- y_groups$original_order[y_groups$group %in% gr_levels[strata_samples$calibration]]
      }
      
      calibration_groups_indices[[jj]] <- sel_sample_indices
      validation_groups_indices[[jj]] <- y_groups$original_order[y_groups$group %in% gr_levels[strata_samples$validation]]
    }
    lengths_list_calibration_indices <- lengths(calibration_groups_indices)
    lengths_list_validation_indices <- lengths(validation_groups_indices)
    calibration_indices <- matrix(NA, max(lengths_list_calibration_indices), number)
    validation_indices <- matrix(NA, max(lengths_list_validation_indices), number)
    colnames(calibration_indices) <- colnames(validation_indices) <- paste0(
      "Resample_",
      seq(1:number)
    )
    rownames(calibration_indices) <- paste0("index_", 1:nrow(calibration_indices))
    rownames(validation_indices) <- paste0("index_", 1:nrow(validation_indices))
    indcs <- 1:nrow(y)
    for (jj in 1:number) {
      ## ramdomly select the replacements for the missing indexes
      ## replacements are necessary!
      n_cal_missing <- max(lengths_list_calibration_indices) - lengths_list_calibration_indices[jj]
      jj_add_list_calibration_indices <- sample(calibration_groups_indices[[jj]],
                                                size = n_cal_missing,
                                                replace = TRUE
      )
      
      n_val_missing <- max(lengths_list_validation_indices) - lengths_list_validation_indices[jj]
      jj_add_list_validation_indices <- sample(validation_groups_indices[[jj]],
                                               size = n_val_missing,
                                               replace = TRUE
      )
      
      calibration_indices[1:lengths_list_calibration_indices[jj], jj] <- calibration_groups_indices[[jj]]
      validation_indices[1:lengths_list_validation_indices[jj], jj] <- validation_groups_indices[[jj]]
      calibration_indices[-(1:lengths_list_calibration_indices[jj]), jj] <- jj_add_list_calibration_indices
      validation_indices[-(1:lengths_list_validation_indices[jj]), jj] <- jj_add_list_validation_indices
    }
  }
  list(
    hold_in = calibration_indices,
    hold_out = validation_indices
  )
}

#' @title A function to assign values to sample distribution strata
#' @description for internal use only! This function takes a continuous variable,
#' creates n strata based on its distribution and assigns the corresponding starta
#' to every value. 
#' @param y a matrix of one column with the response variable.
#' @param n the number of strata.
#' @return a data table with the input \code{y} and the corresponding strata to 
#' every value.
#' @keywords internal
get_sample_strata <- function(y, n) {
  
  y_strata <- unique(quantile(y,
                              probs = seq(0, 1, length = (n + 1)),
                              names = FALSE))
  
  
  strata_labels <- 1:(length(y_strata) - 1)
  y_cuts <- cut(y,
                breaks = y_strata,
                labels = strata_labels,
                include.lowest = TRUE
  )
  
  strata_category <- data.table(
    original_order = 1:length(y),
    strata = y_cuts
  )
  strata_category
}


#' @title A function for stratified calibration/validation sampling 
#' @description for internal use only! This function selects samples
#' based on provided strata.
#' @param original_order a matrix of one column with the response variable.
#' @param starta the number of strata.
#' @param sampling_for sampling to select the calibration samples ("calibration") 
#' or sampling to select the validation samples ("validation").
#' @param replacement logical indicating if sampling with replacement must be 
#' done. 
#' @return a list with the indices of the calibration and validation samples.
#' @keywords internal
get_samples_from_strata <- function(original_order, 
                                    strata, 
                                    sampling_for = c("calibration", "validation"),
                                    replacement = FALSE) {
  
  if (replacement & sampling_for == "validation") {
    n_per_strata <- 2
  } else {
    n_per_strata <- 1
  }
  
  get_random_sample <- function(x, n = n_per_strata) {
    if (length(x) == 1) {
      # this is required to keep the name of the
      # strata, otherwise it fails
      x <- c(x, x)
    }
    sample(x, size = n)
  }
  
  ## for selecting the replacement samples in cases where a strata has only one 
  ## sample, the replacement sample is randomly selected from the data
  
  strata_samples <- tapply(
    X = original_order,
    FUN = get_random_sample,
    INDEX = strata, 
    simplify = FALSE
  )
  
  strata_samples <- do.call("rbind", strata_samples)
  
  if (sampling_for == "validation") {
    if (replacement) {
      unique_sample_strata <- levels(strata)[strata_samples[,1] == strata_samples[,2]]
      if (length(unique_sample_strata > 0)) {
        solve_replacement <- sample(original_order[-strata_samples[,2]], length(unique_sample_strata))
        strata_samples[unique_sample_strata, 2] <- solve_replacement
      }
      
      replacement_indices <- strata_samples[,2]
    } else {
      replacement_indices <- NULL
    }
    
    keep <- original_order[!original_order %in% strata_samples[,1]]
    exclude <- as.vector(sort(strata_samples[,1]))
  } 
  
  if (sampling_for == "calibration") {
    keep <- strata_samples[,1]
    exclude <- original_order[!original_order %in% keep]
    if (replacement) {
      replacement_indices <- sample(keep, length(original_order) - length(keep), replace = TRUE)        
    } else {
      replacement_indices <- NULL
    }
  }
  
  keep <- sort(as.vector(c(keep, replacement_indices)))
  
  strata_samples <- list(calibration = keep,
                         validation = exclude)
  strata_samples
}



