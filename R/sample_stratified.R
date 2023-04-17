#' @title A function to create calibration and validation sample sets for
#' leave-group-out cross-validation
#' @description for internal use only! This is stratified sampling based on the
#' values of a continuous response variable (y). If group is provided, the
#' sampling is done based on the groups and the average of y per group. This
#' function is used to create calibration and validation groups for
#' leave-group-out cross-validations (or
#' leave-group-of-groups-out cross-validation if group argument is provided).
#' @param y a matrix of one column with the response variable.
#' @param p the percentage of samples (or groups if group argument is used) to
#' retain in the validation_indices set
#' @param number the number of sample groups to be crated
#' @param group the labels for each sample in \code{y} indicating the group each
#' observation belongs to.
#' @param replacement A logical indicating sample replacements for the
#' calibration set are required.
#' @param seed an integer for random number generator (default \code{NULL}).
#' @return a list with two matrices (\code{hold_in} and
#' \code{hold_out}) giving the indices of the observations in each
#' column. The number of columns represents the number of sampling repetitions.
#' @keywords internal

sample_stratified <- function(y, p, number, group = NULL, replacement = FALSE, seed = NULL) {

  ## If the percentage of samples to build the hold_in subset is below 50% of
  ## the total number of samples, the selection is based on the number of samples
  ## to retain.
  ## On the other hand, if the percentage of samples to build the hold_in subset
  ## is equal or above 50% of the total number of samples, the selection is based
  ## on the number of samples to exclude.

  # Account for machine precision for calculating the number of samples to take
  # later. The reason is that the floor function can return wrong results, e.g.
  # floor((1 - 0.9) * 50) returns (incorrectly) 4. The following way to round
  # the number of samples to take should solve this issue
  machine_precision <- ifelse(
    is.null(.Machine$sizeof.longdouble) || .Machine$sizeof.longdouble == 0,
    8,
    .Machine$sizeof.longdouble - 2
  )

  if (p < 0.5) {
    p_to_sample <- round(p, machine_precision)
    do_sampling_for <- "calibration"
  } else {
    p_to_sample <- round(1 - p, machine_precision)
    do_sampling_for <- "validation"
  }

  # Stratified sampling if no group is given
  if (is.null(group)) {
    # Number of samples that must be sampled
    nv <- floor(p_to_sample * nrow(y))
    nv <- ifelse(nv < 1, 1, nv)

    # Calculate the number of samples in the validation and calibration sets,
    # respectively.
    n_val <- ifelse(p >= 0.5, nv, nrow(y) - nv)
    n_cal <- ifelse(replacement, nrow(y), nrow(y) - n_val)

    # Compute the strata, including how many samples must be taken in each of them.
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

    # The seed is set both here and if there is a group given. Alternatively,
    # one could set the seed before the if(is.null(group)) statement. However,
    # the way it is done now should ensure that future changes to the code do not
    # change the results of the sampling with the same seeds.
    if (!is.null(seed)) {
      set.seed(seed)
    }

    # Obtain the samples drawn from the calculated strata for a total number of
    # 'number' times.
    for (jj in 1:number) {
      strata_samples <- get_samples_from_strata(
        y = y,
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
    # Stratified sampling if a group is given.
    # The idea is to compute the means of each group and proceed to do the
    # stratified sampling based on the vector of means instead. In particular,
    # this ensures that all members of a groups will always be in the same
    # validation or calibration set.
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

    # Number of groups that must be sampled
    nv <- floor(p_to_sample * n_levels)
    nv <- ifelse(nv < 1, 1, nv)

    # Calculate the number of samples in the validation and calibration sets,
    # respectively.
    n_val <- ifelse(p >= 0.5, nv, nrow(y) - nv)
    n_cal <- ifelse(replacement, nrow(y), nrow(y) - n_val)

    # Compute the strata, including how many groups must be taken in each of them.
    strata_category <- optim_sample_strata(
      y = aggregated_y,
      n = nv
    )

    calibration_groups_indices <- validation_groups_indices <- NULL
    # The seed is set both here and if there is no group given. Alternatively,
    # one could set the seed before the if(is.null(group)) statement. However,
    # the way it is done now should ensure that future changes to the code do not
    # change the results of the sampling with the same seeds.
    if (!is.null(seed)) {
      set.seed(seed)
    }

    # Obtain the groups drawn from the calculated strata for a total number of
    # 'number' times.
    for (jj in 1:number) {
      strata_samples <- get_samples_from_strata(
        y = aggregated_y,
        original_order = strata_category$sample_strata$original_order,
        strata = strata_category$sample_strata$strata,
        samples_per_strata = strata_category$samples_to_get,
        sampling_for = do_sampling_for,
        replacement = replacement
      )

      # Since we sampled based on the group levels, we must convert these levels
      # back to the sample indices.
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
    # Put everything together into the calibration indices.
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
    if (!is.null(seed)) {
      set.seed(seed + 1)
    }
    # Groups can be of different sizes, so the
    for (jj in 1:number) {
      ## randomly select the replacements for the missing indexes
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
get_sample_strata <- function(y, n = NULL, probs = NULL) {
  if (!is.null(n) & !is.null(probs)) {
    stop("both n and probs have been passed to the function, only one of them can be accepted")
  }

  if (!is.null(n)) {
    probs <- seq(0, 1, length = (n + 1))
  }

  y_strata <- unique(
    quantile(
      y,
      probs = probs,
      names = FALSE
    )
  )

  strata_labels <- 1:max(length(y_strata) - 1, 1)
  y_cuts <- cut(
    y,
    breaks = y_strata,
    labels = strata_labels,
    include.lowest = TRUE,
    right = FALSE # Use left closed intervals for compatibility with Julia code
  )

  strata_category <- data.table(
    original_order = 1:length(y),
    strata = y_cuts
  )
  strata_category
}

#' @title A function to construct an optimal strata for the samples, based on
#' the distribution of the given y.
#' @description for internal use only! This function computes the optimal strata
#' from the distribution of the given y
#' @param y a matrix of one column with the response variable.
#' @param n number of samples that must be sampled.
#' @return a list with two \code{data.table} objects: \code{sample_strata} contains
#' the optimal strata, whereas \code{samples_to_get} contains information on how
#' many samples per stratum are supposed to be drawn.
#' @keywords internal
optim_sample_strata <- function(y, n) {
  sample_strata <- get_sample_strata(y, n)
  n_strata <- length(unique(sample_strata$strata))
  iter <- 1
  table_strata <- table(sample_strata$strata)

  # If the number of strata is lower than number of samples, we must correct the
  # strata such that we have the correct number of samples in each stratum.
  #
  # The idea of the following loop is to find the optimal number of strata to
  # properly collect random samples from each of them, as there might be some
  # strata that do not contain any samples. Therefore, we have to widen the
  # strata. As a consequence, the number of strata will be reduced.
  if (n_strata < n || min(table_strata, na.rm = TRUE) < 3) {
    keep_going <- TRUE
    new_n <- n
    new_min_samples_per_strata <- 2 # There must be at least 2 samples per strata
    while (keep_going) {
      # n_vec <- 1:new_n

      # We have to increase the minimum number of samples per strata as we decrease
      # the number of strata. We increase by two; one for the drawing the sample
      # set, one for drawing the replacement set.
      new_min_samples_per_strata <- new_min_samples_per_strata + 2

      # s_missing <- n_vec[-as.numeric(names(table_strata))]
      # v_missing <- rep(0, length(s_missing))
      # names(v_missing) <- s_missing
      # table_strata <- c(table_strata, v_missing)

      # Assign a new number of strata
      new_n <- ceiling(n / (iter + 1))

      # Get the new sample strata for the reduced number of strata
      sample_strata <- get_sample_strata(y, new_n)
      table_strata <- table(sample_strata$strata)

      # Check that each strata contains at least the previously set minimum
      # number of samples. Must be strictly less; otherwise eliminate randomness
      # from the sampling.
      condition_1 <- new_min_samples_per_strata < min(table_strata)
      # Check that all the strata contain samples
      condition_2 <- length(table_strata) == new_n
      # Check that we can further reduce the number of strata without issue
      condition_3 <- new_min_samples_per_strata >= n

      if ((condition_1 & condition_2) | condition_3) {
        break
      }
      iter <- iter + 1
    }

    # Describes how many samples to get for each strata. We start by assigning
    # the same number of samples in each stratum. For some strata, we then reduce
    # this number, as there will (in most cases) more total samples in all strata
    # than the total number of samples available.
    samples_to_get <- data.table(
      strata = levels(sample_strata$strata),
      samples_to_get = rep(
        new_min_samples_per_strata / 2,
        nlevels(sample_strata$strata)
      )
    )
    # Compute how many samples are currently assigned in samples_to_get.
    total_samples <- sum(samples_to_get$samples_to_get)
    # Derive the difference between the total samples in all strata and the number
    # of available samples. In case this number is larger then zero, we must remove
    # some samples to get from the strata.
    to_remove <- total_samples - n
    if (to_remove > 0) {
      # Remove one sample to be drawn from strata with least samples in them
      highest_freqs <- names(sort(table_strata, decreasing = FALSE)[1:to_remove])
      samples_to_get$samples_to_get[samples_to_get$strata %in% highest_freqs] <-
        samples_to_get$samples_to_get[samples_to_get$strata %in% highest_freqs] - 1
    }
  } else {
    # In case the strata already satisfies the above requirements of having exactly
    # n strata and at least 3 samples in each stratum, we do not have to correct
    # the strata further and can get exactly 1 sample in each stratum.
    samples_to_get <- data.table(
      strata = levels(sample_strata$strata),
      samples_to_get = 1
    )
  }
  # Return the strata and how many samples per stratum to get.
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
get_samples_from_strata <- function(y,
                                    original_order,
                                    strata,
                                    samples_per_strata,
                                    sampling_for = c("calibration", "validation"),
                                    replacement = FALSE) {
  # For validation, we double the samples to get per strata in case of replacement.
  # The reason is that we will sample two distinct sets of samples; one denotes
  # the validation indices, while the other describes the replacement indices
  # that will replace the drawn samples. This avoids drawing with replacement
  # for the validation set, which would result in over-optimistic estimates
  # of the RMSE.
  val_replacement <- replacement & sampling_for == "validation"
  if (val_replacement) {
    samples_per_strata$samples_to_get <- 2 * samples_per_strata$samples_to_get
  }

  # This functions draws samples from a given single stratum. We do not draw with
  # replacement here, as the replacements are done manually later.
  # For selecting the replacement samples in cases where a strata has only one
  # sample, the replacement sample is randomly selected from the data
  get_random_sample <- function(x, ns) {
    if (length(x) == 1) {
      # this is required to keep the name of the strata, otherwise it fails
      x <- c(x, x)
    }
    sample(x, size = ns)
  }
  max_samples <- max(samples_per_strata$samples_to_get)
  vec_samples <- rep(NA, max_samples)

  # Here, we draw samples from each strata. In case of sampling for validation
  # with replacement, we draw twice as many samples compared to the other methods.
  # These drawn samples will be split into two sets, one denoting the validation
  # indices, whereas the second one describes the replacement indices that replace
  # the drawn samples.
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
        # Split the drawn samples into the two above described sets. The sets are
        # both contained in vec_samples, but these vectors will be split in half
        # later. Note that we must take special care with the indexing here, to
        # ensure that the samples are put into the correct set.
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
    replacement = val_replacement
  )

  # Put all drawn samples into one single matrix
  strata_samples <- do.call("rbind", strata_samples)
  if (val_replacement) {
    # The samples in case of validation sampling with replacement are separated
    # here into validation and replacement sets. The rows of the matrix
    # samples are the above sampled indices. Hence, the columns correspond to
    # above constructed separation, where the first half are describing the
    # validations set and the rest describe the replacement set.
    col_s <- 1:(ncol(strata_samples) / 2)
    col_replacement <- -col_s
    strata_samples <- cbind(
      sort(strata_samples[, col_s]),
      sort(strata_samples[, col_replacement])
    )
  } else {
    # In any other case, we can simply put them into a matrix
    strata_samples <- as.matrix(sort(strata_samples))
  }


  if (sampling_for == "validation") {
    if (replacement) {
      # Ensure that the validation and replacement sets are distinct.
      unique_sample_strata <- which(strata_samples[, 1] == strata_samples[, 2])
      if (length(unique_sample_strata > 0)) {
        solve_replacement <- sample(original_order[-strata_samples[, 2]], length(unique_sample_strata))
        strata_samples[unique_sample_strata, 2] <- solve_replacement
      }
      replacement_indices <- strata_samples[, 2]
    } else {
      replacement_indices <- NULL
    }
    # These are the calibration and validation sets. The replacement indices
    # are added later to the keep parameter in case of validation sampling
    # with replacement.
    keep <- original_order[!original_order %in% strata_samples[, 1]]
    exclude <- as.vector(sort(strata_samples[, 1]))
  }


  if (sampling_for == "calibration") {
    # For calibration, we can already use the constructed calibration and validation
    # sets. In case of sampling with replacements, the replacement indices are
    # added afterwards
    keep <- strata_samples[, 1]
    exclude <- original_order[!original_order %in% keep]
    if (replacement) {
      # replacement_indices <- sample(keep, length(original_order) - length(keep), replace = TRUE)
      # We draw the replacement indices for calibration from a quartile strata,
      # based on the calibration set.
      replacement_indices <- get_quartile_samples(y[keep], n = length(original_order) - length(keep))
      replacement_indices <- keep[replacement_indices]
    } else {
      replacement_indices <- NULL
    }
  }
  # original
  # if (sampling_for == "calibration") {
  #   keep <- strata_samples[, 1]
  #   exclude <- original_order[!original_order %in% keep]
  #   if (replacement) {
  #     (length(original_order) +length(keep)) / (length(original_order) - length(keep))
  #     ip <-
  #      replacement_indices <-  sample_stratified(y = y[keep, ], p = ip, number = 1, seed = seed)
  #     replacement_indices <- sample(keep, length(original_order) - length(keep), replace = TRUE)
  #   } else {
  #     replacement_indices <- NULL
  #   }
  # }

  # Add replacement indices to the calibration set.
  keep <- sort(as.vector(c(keep, replacement_indices)))

  strata_samples <- list(
    calibration = keep,
    validation = exclude
  )
  strata_samples
}




get_quartile_samples <- function(y, n) {
  quartile_probs <- seq(0, 1, 0.25)
  if (length(y) > 1) {
    strata <- get_sample_strata(y, probs = quartile_probs)
  } else {
    return(NULL)
  }
  nn <- floor(n / length(quartile_probs))

  if (nn < 1) {
    nn <- 1
  }

  ## samples per stratum
  ssampling <- data.frame(
    stratum = sort(as.numeric(unique(strata$strata))),
    n = nn
  )

  # I think this is impossible to reach.
  # if (sum(ssampling$n) > n) {
  #   repl <- length(ssampling$stratum) < sum(ssampling$n) - n
  #   removes <- sample(ssampling$stratum,
  #                     sum(ssampling$n) - n,
  #                     replace = repl)
  #   removes <- table(removes)
  #   for (i in as.numeric(names(adds))) {
  #     if (!is.na(adds[i])) {
  #       i_add <- adds[i]
  #     } else {
  #       i_add <- 0
  #     }
  #     ssampling$n[i] <- ssampling$n[ssampling$stratum == i] - i_add
  #   }
  # }

  if (sum(ssampling$n) < n) {
    repl <- length(ssampling$stratum) < n - sum(ssampling$n)
    adds <- sample(ssampling$stratum, n - sum(ssampling$n), replace = repl)
    if (length(adds) > 1) {
      adds <- table(adds)
      for (i in as.numeric(names(adds))) {
        if (!is.na(adds[i])) {
          i_add <- adds[i]
        } else {
          i_add <- 0
        }
        ssampling$n[i] <- ssampling$n[ssampling$stratum == i] + i_add
      }
    } else {
      ssampling$n[ssampling$stratum == adds] <- ssampling$n[ssampling$stratum == adds] + 1
    }
  }

  strat_samples <- lapply(
    ssampling$stratum,
    FUN = function(stratum_data, stratum_n, i) {
      smp <- stratum_data$original_order[stratum_data$strata == i]
      i_n <- stratum_n$n[stratum_n$stratum == i]
      sample(smp,
        i_n,
        replace = ifelse(i_n > length(smp), TRUE, FALSE)
      )
    },
    stratum_data = strata,
    stratum_n = ssampling
  )

  strat_samples <- sort(unlist(strat_samples))

  if (length(strat_samples) > n) {
    overs <- length(strat_samples) - n
    strat_samples <- strat_samples[-sample(length(strat_samples), overs)]
  }

  if (length(strat_samples) < n) {
    unders <- n - length(strat_samples)
    strat_samples <- c(strat_samples, sort(sample(strat_samples, unders, replace = TRUE)))
  }
  strat_samples
}
