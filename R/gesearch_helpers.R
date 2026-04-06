#' @title iterates over elements to print
#' @description internal. used for printing
#' @keywords internal
cat_iter <- function(x) {
  iter_x <- iter(x, by = "cell", recycle = TRUE)
  
  next_e <- function() {
    nextElem(iter_x)
  }
  next_e
  class(next_e) <- c("isubset", "abstractiter", "iter")
  next_e
}

#' @title fits pls models at each \emph{gesearch} iteration
#' @description internal
#' @keywords internal
biter <- function(
    itersubs,
    Xu,
    Yu,
    y_lim_left,
    y_lim_right,
    iter_sequence,
    n,
    optimization,
    ncomp,
    tune,
    p,
    number,
    scale,
    max_iter,
    tol,
    seed = NULL,
    algorithm, 
    allow_parallel, 
    pchunks = 1
) {
  
  "%mydo%" <- get("%do%")
  if (allow_parallel & getDoParRegistered()) {
    "%mydo%" <- get("%dopar%")
  }
  
  if (is.null(Yu)) {
    n_yu <- 1   
  } else {
    n_yu <- ncol(Yu) 
  }
  
  soptions <- c(
    "rmse_response", 
    "rmse_reconstruction", 
    "score_dissimilarity",
    "deviation_from_boundaries",
    "residual_variance"
  )
  noptions <- length(soptions)
  
  results_temp <- matrix(NA, pchunks, n + (n_yu * noptions))
  
  if (!is.null(colnames(Yu))) {
    nms <- rep(colnames(Yu), each = noptions)
  } else {
    nms <- rep(1:n_yu, each = noptions)
  }
  
  nms <- paste0(
    soptions, 
    "_", 
    nms)
  
  colnames(results_temp) <- c(nms, paste0("sample_", 1:n))
  #   sample = its$idx
  its <- NULL
  rdf <- foreach(
    idx = iter_sequence,
    its = itersubs, 
    .inorder = FALSE,
    .noexport = c("Xr")
  ) %mydo% {
    
    ith_results <- results_temp[1:ncol(its$idx), , drop = FALSE]
    for (h in 1:ncol(its$idx)) {
      hsub_subset <- its$unique_idx %in% its$idx[, h]
      
      # #
      # # Step 2 - sample a training data set of size k from K without replacement
      # #
      # selected.idx <- sample(k_idx,
      #                        size = k,
      #                        replace = FALSE)
      #
      # # retrieve the spectra and dependent variable from the SL for the training
      # set selection
      # X <- Xr[selected.idx,]
      # Y <- Yr[selected.idx]
      
      # dataset <- cbind(X, Y)
      #
      # Step 3 calibrate a PLS model using the selected SL samples
      #
      for (i in  1:ncol(its$y)) {
        if (tune) {
          # perform a leave-group-out cross.validation of the selected SL (k)
          # observations to determine a suitable number of pls factors
          plsv <- pls_cv(
            x = its$x[hsub_subset, , drop = FALSE],
            y = its$y[hsub_subset, i, drop = FALSE],
            ncomp = ncomp,
            method = "pls",
            center = TRUE,
            scale = scale,
            min_component = 1,
            p = p,
            number = number,
            group = its$group,
            retrieve = FALSE,
            tune = TRUE,
            max_iter = max_iter,
            tol = tol, 
            seed = seed, 
            algorithm = algorithm
          )
          
          selected_pls_c <- which.min(plsv$cv_results$rmse_cv)[1]
        } else {
          selected_pls_c <- ncomp
        }
        gs_pls <- opls_gesearch(
          Xr = its$x[hsub_subset, , drop = FALSE],
          Yr = its$y[hsub_subset, i, drop = FALSE],     
          Xu = Xu, 
          ncomp = selected_pls_c,
          scale = scale,     
          response = "response" %in% optimization | "range" %in% optimization,
          reconstruction = "reconstruction" %in% optimization,
          similarity = "similarity" %in% optimization,
          fresponse = "fresponse" %in% optimization,
          algorithm =  algorithm
        )
        
        deviation_from_boundaries <- residual_variance <- rmse_reconstruction <- rmse_response <- score_dissimilarity <- NULL
        if ("similarity" %in% optimization) {
          score_dissimilarity <- gs_pls$score_dissimilarity[[1]]
        } else {
          score_dissimilarity <- NA
        }
        
        if ("response" %in% optimization) {
          # It is possible we get NA predictions, catch this case and ignore this
          # observation test
          y_pred <- gs_pls$pred_response[, 1]
          if (any(is.na(y_pred))) {
            
          } else {
            rmse_response <- sqrt(mean((y_pred - Yu[, i])^2, na.rm = TRUE))
          }
        } else {
          rmse_response <- NA
        }
        
        if ("reconstruction" %in% optimization) {
          rmse_reconstruction <- gs_pls$rmse_reconstruction[[1]]
        } else {
          rmse_reconstruction <- NA
        }
        
        if ("range" %in% optimization) {
          out_g_left <- gs_pls$pred_response[which(gs_pls$pred_response < y_lim_left), ]
          out_g_right <- gs_pls$pred_response[which(gs_pls$pred_response > y_lim_right), ]
          ws_left <- sum(y_lim_left - out_g_left) / length(gs_pls$pred_response)
          ws_right <- sum(out_g_right - y_lim_right) / length(gs_pls$pred_response)
          deviation_from_boundaries <- ws_out <- ws_left + ws_right
          
        } else {
          deviation_from_boundaries <- NA
        }
        
        if ("fresponse" %in% optimization) {
          residual_variance <- gs_pls$residual_variance[[1]]
        } else {
          residual_variance <- NA
        }
        
        
        # combine the validation results and the observations that were selected for this
        # model
        ith_results[h, (i * noptions) - 4] <- rmse_response
        ith_results[h, (i * noptions) - 3] <- rmse_reconstruction
        ith_results[h, (i * noptions) - 2] <- score_dissimilarity
        ith_results[h, (i * noptions) - 1] <- deviation_from_boundaries
        ith_results[h, (i * noptions)] <- residual_variance
      }
      ith_results[h, -c(1:(n_yu * noptions))] <- its$unique_idx[hsub_subset]
    }
    ith_results
  }
  
  # compile a data frame of all results of the sampling iterations
  rdf <- do.call("rbind", rdf)
  
  obj_cols <- NULL
  if ("response" %in% optimization) {
    obj_cols <- grep("rmse_response", colnames(rdf))
  } 
  
  if ("reconstruction" %in% optimization) {
    obj_cols <- sort(c(obj_cols, grep("rmse_reconstruction", colnames(rdf))))
  }
  
  if ("similarity" %in% optimization) {
    obj_cols <- sort(c(obj_cols, grep("score_dissimilarity", colnames(rdf))))
  }
  
  if ("range" %in% optimization) {
    obj_cols <- sort(c(obj_cols, grep("deviation_from_boundaries", colnames(rdf))))
  }
  
  
  if ("fresponse" %in% optimization) {
    obj_cols <- sort(c(obj_cols, grep("residual_variance", colnames(rdf))))
  }
  
  weakness_scores_subset <- rdf[, obj_cols, drop = FALSE]
  rdf <- rdf[, -grep("^rmse_|^score_|^deviation_|^residual_", colnames(rdf)), drop = FALSE]
  return(list(
    weakness_scores_subset = weakness_scores_subset,
    sampleidx = rdf
  ))
}

#' @title iterator that subsets a matrix based on input indices for pls modeling
#' at each \emph{gesearch} iteration
#' @description internal
#' @param x an input matrix of predictors
#' @param y a response variable
#' @param group a variable giving the groups/calsses to which the observations
#' belong to (used for avoiding pseudo-replication during validation)
#' @param indx the indices to retrieve
#' @param chunksize the number of elements of \code{indx} to return with each 
#' iteration.
#' @return an object of \code{class} iterator giving a \code{list} with (a) a
#' subset of the input matrix with rows re-arranged according to `indx`
#' @details this is designed to be used within (parallel) loops to avoid sending
#' entire matrices to each core. It just sends the subset of that is required.
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
ithrssubsets <- function(
    x,
    y,
    group = NULL,
    indx, 
    chunksize = 1L
) {
  it_indx <- iter(indx, by = "column", chunksize = chunksize)
  
  if (is.null(group)) {
    this_next_element <- function() {
      ss <- nextElem(it_indx)
      unique_ss <- unique(as.vector(ss))
      list(
        x = x[unique_ss, , drop = FALSE],
        y = y[unique_ss, , drop = FALSE],
        group = NULL,
        idx = ss, 
        unique_idx = unique_ss
      )
    }
  } else {
    nextEl <- function() {
      ss <- nextElem(it_indx)
      unique_ss <- unique(as.vector(ss))
      list(
        x = x[unique_ss, , drop = FALSE],
        y = y[unique_ss, , drop = FALSE],
        group = factor(group[unique_ss]),
        idx = ss, 
        unique_idx = unique_ss
      )
    }
  }
  obj <- list(nextElem = this_next_element)
  class(obj) <- c("isubset", "abstractiter", "iter")
  obj
}

#' @title drop the samples that consistently appear in high error models
#' @description internal
#' @keywords internal
drop_indices <- function(
    uu, 
    vv,
    r,
    obs_indices,
    resampling_indices, 
    weakness_subset,
    retain_by, 
    cull_quantity, 
    percentile_type, 
    max_ncomp
) {
  
  # vv <- as.vector(table(resampling_indices))
  # uu <- vv <- rep(0, n)
  
  for (i in 1:nrow(resampling_indices)) {
    uu[resampling_indices[i, ]] <- uu[resampling_indices[i, ]] + weakness_subset[i]
    vv[resampling_indices[i, ]] <- vv[resampling_indices[i, ]] + 1
  }
  
  vv[vv == 0] <- NA
  # normalise U using V
  uu <- uu / vv
  
  # uu <- sapply(obs_indices, 
  #              FUN = function(weakness, indices, ii) {
  #                mean(weakness[arrayInd(which(indices == ii), .dim = nrow(indices))])
  #              },
  #              weakness = weakness_subset, 
  #              indices = resampling_indices)
  # uu[is.nan(uu)] <- NA
  
  u_df <- data.frame(
    idx = obs_indices,
    weakness_score = uu
  )
  
  # now order / rank by weakness
  u_ordered <- u_df[order(u_df$weakness_score, decreasing = TRUE), ]
  u_ordered <- u_ordered[!is.na(u_ordered$weakness_score), ]
  
  uws <- length(unique(u_ordered$weakness_score))
  if (uws == 1) {
    results <- list(interrupted = FALSE,
                    s_to_drop = 0,
                    cutoff = 0,
                    max_weakness_score = 0,
                    obs_indices = obs_indices)
    return(results)
  } else {
    if (retain_by == "proportion") {
      # select the poorest performing SL samples, the amount is determined by
      # cull_quantity calculated earlier
      if (uws < cull_quantity) {
        cull_quantity <- uws
      }
      worst_reference_set <- u_ordered[1:cull_quantity, ]
      ok_reference_set <- u_ordered[-c(1:cull_quantity), ]
      cutoff <- min(worst_reference_set$weakness_score)
      pp <- 1 - nrow(ok_reference_set) / nrow(u_ordered)
    } else {
      ## Instead removing a given number of samples
      ## remove samples that are above a given cutoff prob
      cutoff <- quantile(uu, 1 - r, na.rm = TRUE, type = percentile_type)
      worst_reference_set <- u_ordered[u_ordered$weakness_score >= cutoff, ]
      ok_reference_set <- u_ordered[u_ordered$weakness_score < cutoff, ]
      if (nrow(ok_reference_set) < max_ncomp) {
        mss <- (paste0(
          "Iteration interrupted as the number of selected observations (",
          nrow(ok_reference_set),
          ") is lower than the number of pls factors (", max_ncomp, ")"
        ))
        return(list(interrupted = TRUE,
                    message = mss))
      }
      pp <- 1 - nrow(ok_reference_set) / nrow(u_ordered)
    }
    obs_indices[worst_reference_set$idx] <- NA
    results <- list(interrupted = FALSE,
                    s_to_drop = nrow(worst_reference_set),
                    cutoff = cutoff,
                    max_weakness_score = max(ok_reference_set$weakness_score),
                    obs_indices = obs_indices)
    return(results)
  }
}

#' @title combine two randomly selected individuals
#' @description internal
#' @keywords internal
combine_individuals <- function(x, k, k_idx) {
  # Select two individuals out of x (i.e. columns), randomly (without replacements)
  sel_individuals <- sample(1:ncol(x), 2, replace = FALSE)
  # Combine the two samples and remove duplicated entries and NAs
  combined_individuals <- unique(na.omit(c(x[, sel_individuals])))
  # If the combined individual has more than k samples, we draw k samples from 
  # the combined gene pool
  # Otherwise, do mutation: randomly selected genes from the global gene pool 
  # are used to fill up the number of genes to k
  if (length(combined_individuals) < k) {
    mutation_gene_pool <- setdiff(k_idx, combined_individuals)
    c(
      combined_individuals, 
      sample(
        mutation_gene_pool, 
        k - length(combined_individuals), 
        replace = FALSE
      )
    )
  } else {
    sample(combined_individuals, k, replace = FALSE)
  }
}
