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

#' @title Check and Process `pc_selection` Argument
#' @description
#' Internal utility to validate and normalize the `pc_selection` argument,
#' used for controlling how the number of components is selected.
#'
#' @details
#' Supports the following selection methods:
#' - `"opc"` or `"manual"`: manual control of component count (defaults to 40)
#' - `"cumvar"`: based on cumulative explained variance (defaults to 0.99)
#' - `"var"`: based on per-component variance (defaults to 0.01)
#'
#' Handles automatic defaulting when values are missing, and validates
#' against the shape of the input data (`n_rows_x`, `n_cols_x`).
#'
#' @param n_rows_x Number of rows in the reference dataset
#' @param n_cols_x Number of columns in the reference dataset
#' @param pc_selection List or character. See details.
#' @param default_max_comp Default maximum components to test (default: 40)
#' @param default_max_cumvar Default cumulative variance threshold (default: 0.99)
#' @param default_max_var Default single component variance threshold (default: 0.01)
#'
#' @return A list with:
#' \describe{
#'   \item{pc_selection_checked}{Validated list with `method` and `value`}
#'   \item{max_comp}{Maximum number of components to consider}
#' }
#' 
#' @keywords internal
#' @noRd
check_pc_arguments <- function(
    n_rows_x, n_cols_x, pc_selection,
    default_max_comp = 40,
    default_max_cumvar = 0.99,
    default_max_var = 0.01
) {
  # Validate pc_selection structure
  if (!is.list(pc_selection)) {
    stop(
      "The 'pc_selection' argument must be a list with two elements: ",
      "method and value. Optionally, a character string with only the ",
      "method can be used (value will be set automatically)."
    )
  }
  
  pc_method <- pc_selection[[1]]
  pc_value <- if (length(pc_selection) > 1) pc_selection[[2]] else NULL
  max_dim <- min(n_rows_x, n_cols_x)
  
  if (pc_method %in% c("opc", "manual")) {
    # Assign default if value is missing
    if (is.null(pc_value)) {
      value <- min(default_max_comp, max_dim)
      message(
        "Missing value in 'pc_selection': maximum components automatically set to ", 
        value
      )
    } else {
      if (!is.numeric(pc_value)) {
        stop("The 'value' in 'pc_selection' must be numeric for method '", pc_method, "'.")
      }
      value <- floor(pc_value)
      if (value < 1 || value > max_dim) {
        stop("The number of components must be between 1 and ", max_dim, ".")
      }
    }
    max_comp <- value
  } else if (pc_method %in% c("cumvar", "var")) {
    # Assign default if value is missing
    if (is.null(pc_value)) {
      value <- if (pc_method == "cumvar") default_max_cumvar else default_max_var
      msg <- if (pc_method == "cumvar") {
        "Missing value in 'pc_selection': cumulative variance automatically set to 0.99 (99%)"
      } else {
        "Missing value in 'pc_selection': minimum component variance automatically set to 0.01 (1%)"
      }
      message(msg)
    } else {
      if (!is.numeric(pc_value)) {
        stop("The 'value' in 'pc_selection' must be numeric for method '", pc_method, "'.")
      }
      if (pc_value <= 0 || pc_value > 1) {
        stop("For method '", pc_method, "', 'value' must be > 0 and ≤ 1.")
      }
      value <- pc_value
    }
    max_comp <- max_dim - 1
  } else {
    stop("Invalid 'pc_selection' method: must be one of 'opc', 'manual', 'cumvar', or 'var'.")
  }
  
  list(
    pc_selection_checked = list(method = pc_method, value = value),
    max_comp = max_comp
  )
}
