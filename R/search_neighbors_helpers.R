#' @title From dissimilarity matrix to neighbors
#' @param diss_matrix a matrix representing the dissimilarities between 
#' observations in a matrix \code{Xu} and observations in another matrix 
#' \code{Xr}. \code{Xr} in rows \code{Xu} in columns.
#' @param k an integer value indicating the k-nearest neighbors of each 
#' observation in \code{Xu} that must be selected from \code{Xr}.
#' @param k_diss an integer value indicating a dissimilarity treshold. 
#' For each observation in \code{Xu}, its nearest neighbors in \code{Xr} 
#' are selected as those for which their dissimilarity to \code{Xu} is below 
#' this \code{k_diss} threshold. This treshold depends on the corresponding 
#' dissimilarity metric specified in \code{diss_method}. Either \code{k} or 
#' \code{k_diss} must be specified.
#' @param k_range an integer vector of length 2 which specifies the minimum 
#' (first value) and the maximum (second value) number of neighbors to be 
#' retained when the \code{k_diss} is given.
#' @param spike a vector of integers indicating what observations in \code{Xr} 
#' (and \code{Yr}) must be 'forced' to always be part of all the neighborhoods.
#' @param return_dissimilarity logical indicating if the input dissimilarity 
#' must be mirroed in the output.
#' @description internal
#' @keywords internal
diss_to_neighbors <- function(diss_matrix, 
                              k = NULL, k_diss = NULL, k_range = NULL,
                              spike = NULL,
                              return_dissimilarity = FALSE) {
  
  
  if (!is.null(spike)) {
    
    f_order_neigh <- function(x, s) {
      x <- order(x)
      c(s, x[!x %in% s])
    }
    
    neighbor_indcs <- apply(diss_matrix, 
                            MARGIN = 2, 
                            FUN = f_order_neigh, 
                            s = spike)
    stats <- seq(1, 
                 nrow(diss_matrix) * ncol(diss_matrix), 
                 by = nrow(diss_matrix)) - 1
    neighbors_diss <- sweep(neighbor_indcs, 
                            MARGIN = 2, 
                            STATS = stats, 
                            FUN = "+")
    neighbors_diss[] <- diss_matrix[as.vector(neighbors_diss)]
  } else {
    neighbor_indcs <- apply(diss_matrix, MARGIN = 2, FUN = order)
    neighbors_diss <- apply(diss_matrix, MARGIN = 2, FUN = sort)
  }
  
  
  if (!is.null(k)) {
    neighbor_indcs <- neighbor_indcs[1:k, , drop = FALSE]
    neighbors_diss <- neighbors_diss[1:k, , drop = FALSE]
  } else {
    k_min <- min(k_range)
    k_max <- max(k_range)
    neighbors <- neighbors_diss <= k_diss
    if (!is.null(spike)) {
      neighbors[1:length(spike),] <- TRUE
    }
    n_neighbors <- original_n_neighbors <- colSums(neighbors)
    n_neighbors[original_n_neighbors < k_min] <- k_min
    n_neighbors[original_n_neighbors > k_max] <- k_max
    neighbor_indcs <- neighbor_indcs[1:max(n_neighbors), , drop = FALSE]
    neighbors_diss <- neighbors_diss[1:max(n_neighbors), , drop = FALSE]
    
    
    neigh_list <- lapply(n_neighbors, FUN = function(x) 1:x)
    sts <- (max(n_neighbors) * (1:length(n_neighbors))) - max(n_neighbors)
    neighbor_vec_indcs <- sort(unlist(Map('+', neigh_list, sts)))
    neighbor_indcs[-neighbor_vec_indcs] <- NA
    neighbors_diss[-neighbor_vec_indcs] <- NA
  }
  rownames(neighbor_indcs) <- paste0("k_", 1:nrow(neighbor_indcs))
  rownames(neighbors_diss) <- rownames(neighbor_indcs)
  
  colnames(neighbor_indcs) <- paste0("Xu_", 1:ncol(neighbor_indcs))
  colnames(neighbors_diss) <- colnames(neighbor_indcs)
  
  
  results <- list(neighbors_diss = neighbors_diss,
                  neighbors = neighbor_indcs,
                  unique_neighbors = sort(unique(as.vector(neighbor_indcs))))
  
  if (!is.null(k_diss)) {
    results$k_diss_info <- data.table(Xu_index = 1:ncol(diss_matrix), 
                                      n_k = original_n_neighbors,
                                      final_n_k = n_neighbors)
  }
  
  if(return_dissimilarity) {
    results$dissimilarity <- diss_matrix
  }
  results
}
