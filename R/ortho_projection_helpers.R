#' @title evaluation of multiple distances obtained with multiple PCs
#' @description internal
#' @keywords internal
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
    tmp <- sim_eval(
      d = d,
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

#' @title checks the pc_selection argument
#' @description internal
#' @keywords internal
check_pc_arguments <- function(n_rows_x, n_cols_x, pc_selection,
                               default_max_comp = 40,
                               default_max_cumvar = 0.99,
                               default_max_var = 0.01) {
  pc_selection_method <- pc_selection[[1]]

  if (pc_selection_method %in% c("opc", "manual")) {
    if (length(pc_selection) == 1) {
      treshold_comp <- min(n_rows_x, n_cols_x)
      treshold_comp <- if_else(treshold_comp > default_max_comp,
        default_max_comp, treshold_comp
      )

      pc_selection_checked <- list(
        method = pc_selection_method,
        value = treshold_comp
      )

      message(paste0(
        "Missing value in 'pc_selection', maximum components to be ",
        "tested automatically set to ", treshold_comp
      ))
    } else {
      if (!is.list(pc_selection)) {
        stop(paste0(
          "The 'pc_selection' argument must be a list in which the ",
          "first object indicates the selection method and the second ",
          "object indicates the parameter value of the method. ",
          "Optionally, a character string specifiying only the ",
          "method can be used, in this case the parameter value is ",
          "set automatically"
        ))
      }
      pc_selection_checked <- list(
        method = pc_selection[[1]],
        value = floor(pc_selection[[2]])
      )
      if (!is.numeric(pc_selection_checked$value)) {
        stop("The second value in pc_selection must be an integer value")
      }
      if (pc_selection_checked$value < 2 | pc_selection_checked$value > min(n_rows_x, n_cols_x)) {
        stop(paste(
          "The maximum number of principal components must be a value ",
          " between 2 and", min(n_rows_x, n_cols_x)
        ))
      }
    }
    max_comp <- pc_selection_checked$value
  }

  if (pc_selection_method %in% c("cumvar", "var")) {
    if (length(pc_selection) == 1) {
      if (pc_selection_method == "cumvar") {
        pc_selection_checked <- list(
          method = "cumvar",
          value = default_max_cumvar
        )
        message(paste(
          "Missing value in 'pc_selection', amount of cumulative ",
          "variance to be retained automatically set to 0.99 (99%)"
        ))
      } else {
        pc_selection_checked <- list(method = "var", value = default_max_var)
        message(paste0(
          "Since the value of the pc_selection argument is missing, ",
          "retaining components that expain at least 0.01 (1%) ",
          "of the original variance"
        ))
      }
    } else {
      pc_selection_checked <- list(
        method = pc_selection[[1]],
        value = pc_selection[[2]]
      )

      if (!is.numeric(pc_selection_checked$value)) {
        stop("The second element in 'pc_selection' must be a numeric value")
      }
      if (pc_selection_checked$value > 1 | pc_selection_checked$value <= 0) {
        stop(paste0(
          "When the method for 'pc_selection' is either 'var' or ",
          "'cumvar' the value in pc_selection must be a number larger ",
          "than 0 and below or equal to 1"
        ))
      }
    }
    max_comp <- min(n_rows_x, n_cols_x) - 1
  }

  list(
    pc_selection_checked = pc_selection_checked,
    max_comp = max_comp
  )
}
