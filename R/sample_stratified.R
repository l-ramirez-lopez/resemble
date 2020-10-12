#' @title A function to create calibration and validation sample sets for
#' leave-group-out and bootstraping cross-validation
#' @description for internal use only! This is stratified sampling based on the
#' values of a continuous response variable (y). If group is provided, the
#' sampling is done based on the groups and the average of y per group. This
#' function is used to create calibration and validation groups for
#' leave-group-out cross-validations (or
#' leave-group-of-groups-out cross-validation if group argument is provided)
#' and/or bootstrapping.
#' @param y a matrix of one column with the response variable.
#' @param p the percentage of samples (or groups if group argument is used) to
#' retain in the validation_indices set
#' @param number the number of sample groups to be crated
#' @param group the labels for each sample in \code{y} indicating the group each
#' observation belongs to.
#' @param replacement A logical indicating sample replacements for the
#' calibration set are required.
#' @return a list with two matrices (\code{hold_in} and
#' \code{hold_out}) giving the indices of the observations in each
#' column. The number of columns represents the number of sampling repetitions.
#' @keywords internal

sample_stratified <- function(y, p, number, group = NULL, replacement = FALSE) {
  
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
    
    strata_category <- optim_sample_strata(
      y = y,
      n = nv
    )
    
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
      strata_samples <- get_samples_from_strata(
        original_order = strata_category$sample_strata$original_order,
        strata = strata_category$sample_strata$strata,
        samples_per_strata = strata_category$samples_to_get,
        replacement = replacement,
        sampling_for = do_sampling_for
      )
      
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
    
    strata_category <- optim_sample_strata(
      y = aggregated_y,
      n = nv
    )
    strata_category$strata <- paste0("l_", strata_category$strata)
    
    calibration_indices <- matrix(0, nv, number)
    colnames(calibration_indices) <- paste0("Resample_", seq(1:number))
    rownames(calibration_indices) <- paste0("index_", seq(1:nv))
    
    calibration_groups_indices <- validation_groups_indices <- NULL
    
    for (jj in 1:number) {
      
      strata_samples <- get_samples_from_strata(
        original_order = strata_category$sample_strata$original_order,
        strata = strata_category$sample_strata$strata,
        samples_per_strata = strata_category$samples_to_get,
        sampling_for = do_sampling_for,
        replacement = replacement
      )
      
      if (replacement) {
        sel_sample_indices <- lapply(1:length(strata_samples$calibration),
                                     FUN = function(ith, full_groups, cal_groups) {
                                       ## this equivalent to extract from
                                       ## y_groups$original_order
                                       which(full_groups %in% cal_groups[ith])
                                     },
                                     full_groups = y_groups$group,
                                     cal_groups = gr_levels[strata_samples$calibration]
        )
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
                              names = FALSE
  ))
  
  
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
}


optim_sample_strata <- function(y, n) {
  sample_strata <- get_sample_strata(y, n)
  n_strata <- length(unique(sample_strata$strata))
  iter <- 1
  table_strata <- table(sample_strata$strata)
  if (n_strata < n) {
    keep_going <- TRUE
    new_n <- n
    new_min_samples_per_strata <- 2 ## at least 2 samples per strata
    while (keep_going) {
      n_vec <- 1:new_n
      new_min_samples_per_strata <- new_min_samples_per_strata + 2
      s_missing <- n_vec[-as.numeric(names(table_strata))]
      v_missing <- rep(0, length(s_missing))
      names(v_missing) <- s_missing
      table_strata <- c(table_strata, v_missing)
      
      new_n <- ceiling(n / (iter + 1)) # ????????????????
      
      sample_strata <- get_sample_strata(y, new_n)
      table_strata <- table(sample_strata$strata)
      
      condition_1 <- new_min_samples_per_strata < min(table_strata)
      condition_2 <- length(table(sample_strata$strata)) == new_n
      condition_3 <- new_min_samples_per_strata >= n
      
      
      if ((condition_1 & condition_2) | condition_3) {
        break
      }
      iter <- iter + 1
    }
    
    samples_to_get_no_replacement <- rep(
      new_min_samples_per_strata / 2,
      nlevels(sample_strata$strata)
    )
    
    
    samples_to_get <- data.table(
      strata = levels(sample_strata$strata),
      samples_to_get = samples_to_get_no_replacement
    )
    total_samples <- new_n * (iter + 1)
    to_remove <- total_samples - n
    if (to_remove > 0) {
      highest_freqs <- names(sort(table_strata, decreasing = TRUE)[1:to_remove])
      samples_to_get$samples_to_get[samples_to_get$strata %in% highest_freqs] <-
        samples_to_get$samples_to_get[samples_to_get$strata %in% highest_freqs] - 1
    }
  } else {
    samples_to_get <- data.table(
      strata = levels(sample_strata$strata),
      samples_to_get = 1
    )
  }
  
  list(
    sample_strata = sample_strata,
    samples_to_get = samples_to_get
  )
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
                                    samples_per_strata,
                                    sampling_for = c("calibration", "validation"),
                                    replacement = FALSE) {
  n_samples <- sum(samples_per_strata$samples_to_get)
  with_replacement <- FALSE
  if (replacement & sampling_for == "validation") {
    samples_per_strata$samples_to_get <- 2 * samples_per_strata$samples_to_get
    with_replacement <- TRUE
  }
  
  get_random_sample <- function(x, ns) {
    if (length(x) == 1) {
      # this is required to keep the name of the
      # strata, otherwise it fails
      x <- c(x, x)
    }
    sample(x, size = ns)
  }
  
  ## for selecting the replacement samples in cases where a strata has only one
  ## sample, the replacement sample is randomly selected from the data
  max_samples <- max(samples_per_strata$samples_to_get)
  vec_samples <- rep(NA, max_samples)
  
  strata_samples <- lapply(levels(strata),
                           FUN = function(strata,
                                          original_order,
                                          samples_per_strata,
                                          vec_samples,
                                          replacement,
                                          ii) {
                             ith_n <- samples_per_strata$samples_to_get[samples_per_strata$strata == ii]
                             ith_set <- original_order[which(strata == ii)]
                             ith_sel <- get_random_sample(ith_set, ith_n)
                             if (replacement) {
                               ln <- (length(ith_sel) / 2)
                               vec_samples[1:ln] <- ith_sel[1:ln]
                               vec_samples[(length(vec_samples) - ln + 1):length(vec_samples)] <- ith_sel[-(1:ln)]
                             } else {
                               vec_samples[1:length(ith_sel)] <- ith_sel
                             }
                             
                             vec_samples
                           },
                           strata = strata,
                           original_order = original_order,
                           samples_per_strata = samples_per_strata,
                           vec_samples = vec_samples,
                           replacement = with_replacement
  )
  
  strata_samples <- do.call("rbind", strata_samples)
  if (replacement & sampling_for == "validation") {
    col_s <- 1:(ncol(strata_samples) / 2)
    col_replacement <- -col_s
    strata_samples <- cbind(
      sort(strata_samples[, col_s]),
      sort(strata_samples[, col_replacement])
    )
  } else {
    strata_samples <- as.matrix(sort(strata_samples))
  }
  
  
  if (sampling_for == "validation") {
    if (replacement) {
      unique_sample_strata <- levels(strata)[strata_samples[, 1] == strata_samples[, 2]]
      if (length(unique_sample_strata > 0)) {
        solve_replacement <- sample(original_order[-strata_samples[, 2]], length(unique_sample_strata))
        strata_samples[unique_sample_strata, 2] <- solve_replacement
      }
      
      replacement_indices <- strata_samples[, 2]
    } else {
      replacement_indices <- NULL
    }
    
    keep <- original_order[!original_order %in% strata_samples[, 1]]
    exclude <- as.vector(sort(strata_samples[, 1]))
  }
  
  if (sampling_for == "calibration") {
    keep <- strata_samples[, 1]
    exclude <- original_order[!original_order %in% keep]
    if (replacement) {
      replacement_indices <- sample(keep, length(original_order) - length(keep), replace = TRUE)
    } else {
      replacement_indices <- NULL
    }
  }
  
  keep <- sort(as.vector(c(keep, replacement_indices)))
  
  strata_samples <- list(
    calibration = keep,
    validation = exclude
  )
  strata_samples
}
