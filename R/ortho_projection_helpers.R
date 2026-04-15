#' @title evaluation of multiple distances obtained with multiple PCs
#' @description internal
#' @keywords internal
#' @noRd
eval_multi_pc_diss <- function(scores,
                               side_info,
                               from = 1,
                               to = ncol(scores),
                               steps = 1,
                               method = c("pc", "pls"),
                               check_dims = TRUE) {
  if (ncol(side_info) == 1) {
    extract_sim_results <- function(x) {
      ## this takes either rmsd or kappa as only one will be in the matrix
      measure_names <- colnames(x$eval)[colnames(x$eval) %in% c("rmsd", "kappa")]
      return(list(
        result = x$eval[, measure_names],
        measure_names = measure_names
      ))
    }
    n_cols_results <- 1
  } else {
    extract_sim_results <- function(x) {
      measure_names <- c(
        paste0("rmsd_", names(x$eval[, "rmsd"])),
        "mean_standardized_rmsd_Yr"
      )
      return(list(
        result = c(
          x$eval[, "rmsd"],
          x$global_eval[, "mean_standardized_rmsd"]
        ),
        measure_names = measure_names
      ))
    }
    n_cols_results <- ncol(side_info) + 1
  }

  eval_pcs <- seq(from = 1, to = to, by = steps)
  results <- matrix(NA, length(eval_pcs), n_cols_results)

  std_scores <- sweep(scores,
    MARGIN = 2,
    STATS = get_column_sds(scores),
    FUN = "/"
  )

  if (check_dims) {
    if (nrow(side_info) != nrow(scores)) {
      stop(paste0(
        "Number of observations in 'scores' do not match with the ",
        "number of observations in 'side_info'"
      ))
    }
  } else {
    std_scores <- std_scores[1:nrow(side_info), , drop = FALSE]
  }

  d <- 0
  for (i in eval_pcs) {
    d <- d + fast_diss_vector(
      std_scores[, i, drop = TRUE]
    )
    tmp <- diss_evaluate(
      diss = d,
      side_info = side_info
    )
    ith_result <- extract_sim_results(tmp)
    results[i, 1:n_cols_results] <- unlist(ith_result$result)
  }

  colnames(results) <- ith_result$measure_names
  eval_pcs <- matrix(eval_pcs, dimnames = list(eval_pcs, method))
  results <- cbind(pc = eval_pcs, results)

  if (ncol(side_info) == 1) {
    colnames(results) <- gsub("rmsd", "rmsd_Yr", colnames(results))
  }

  if ("kappa" %in% ith_result$measure_names) {
    best_pc <- eval_pcs[which.max(results[, "kappa"])]
  }
  if ("mean_standardized_rmsd_Yr" %in% ith_result$measure_names) {
    best_pc <- eval_pcs[which.min(results[, "mean_standardized_rmsd_Yr"])]
  }
  if ("rmsd" %in% ith_result$measure_names) {
    best_pc <- eval_pcs[which.min(results[, "rmsd_Yr"])]
  }

  return(list(results = results, best_pc = best_pc))
}
