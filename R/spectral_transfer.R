#' Build a Spectral Transfer Model
#'
#' Creates a model to transfer spectra from a source instrument (e.g., low-cost)
#' to match a target instrument (e.g., high-end). Supports Piecewise Direct
#' Standardization (PDS) and Subtraction methods.
#'
#' @param source_spectra A numeric matrix of spectra from the source instrument
#'   (samples as rows).
#' @param target_spectra A numeric matrix of spectra from the target instrument.
#'   Must have the same dimensions and sample order as `source_spectra`.
#' @param method The transfer method to use. Can be either `"pds"` (the default)
#'   or `"subtraction"`.
#' @param ... Additional arguments passed to the specific method. For PDS, these
#'   are `n_components` and `window_size`. For Subtraction, this is `outlier_wl`
#'   (the wavelength used for outlier detection).
#'
#' @return A `transfer_model` object containing the necessary components to apply
#'   the transfer to new data.
#'
#' @importFrom stats median quantile IQR lm coef
#' @importFrom pls plsr
#' @export
#'
#' @examples
#' \dontrun{
#' # Build a PDS model with 2 components and a window size of 5
#' pds_model <- build_transfer_model(isc_spectra, asd_spectra, method = "pds",
#'                                   n_components = 2, window_size = 5)
#'
#' # Build a Subtraction model, using wavelength 1415 for outlier removal
#' sub_model <- build_transfer_model(isc_spectra, asd_spectra, method = "subtraction",
#'                                   outlier_wl = "1415")
#' }
build_transfer_model <- function(source_spectra, target_spectra, method = "pds", ...) {
  # --- Input validation ---
  if (!identical(dim(source_spectra), dim(target_spectra))) {
    stop("Source and target spectra must have the same dimensions.", call. = FALSE)
  }

  # Convert to matrix
  source_matrix <- as.matrix(source_spectra)
  target_matrix <- as.matrix(target_spectra)
  args <- list(...)

  if (method == "pds") {
    # --- PDS Method ---
    n_comp <- args$n_components %||% 3
    window_size <- args$window_size %||% 5
    k <- floor(window_size / 2)
    num_wavs <- ncol(target_matrix)

    transfer_matrix <- matrix(0, nrow = num_wavs, ncol = num_wavs)
    intercepts <- numeric(num_wavs)

    for (i in (k + 1):(num_wavs - k)) {
      # Define the local window
      window_indices <- (i - k):(i + k)

      # Perform local PLS regression
      fit <- pls::plsr(target_matrix[, i] ~ source_matrix[, window_indices],
                       ncomp = n_comp, scale = FALSE)

      # Extract coefficients
      pls_coefs <- as.numeric(stats::coef(fit, ncomp = n_comp, intercept = TRUE))

      # Store in transfer matrix and intercept vector
      intercepts[i] <- pls_coefs[1]
      transfer_matrix[window_indices, i] <- pls_coefs[-1]
    }

    model <- list(
      transfer_matrix = transfer_matrix,
      intercepts = intercepts,
      method = "pds"
    )

  } else if (method == "subtraction") {
    # --- Subtraction Method ---
    outlier_wl <- args$outlier_wl %||% "1415"

    if (!outlier_wl %in% colnames(source_matrix)) {
      stop("The specified outlier_wl is not a valid column name.", call. = FALSE)
    }

    differences <- target_matrix - source_matrix

    # Outlier removal based on the specified wavelength
    quartiles <- stats::quantile(differences[, outlier_wl], probs = c(0.25, 0.75), na.rm = TRUE)
    iqr_val <- stats::IQR(differences[, outlier_wl], na.rm = TRUE)
    lower_bound <- quartiles[1] - 1.5 * iqr_val
    upper_bound <- quartiles[2] + 1.5 * iqr_val

    is_outlier <- differences[, outlier_wl] < lower_bound | differences[, outlier_wl] > upper_bound

    # Calculate median correction factor from non-outlier samples
    correction_vector <- apply(differences[!is_outlier, ], 2, stats::median, na.rm = TRUE)

    model <- list(
      correction_vector = correction_vector,
      method = "subtraction"
    )

  } else {
    stop("Invalid method specified. Choose 'pds' or 'subtraction'.", call. = FALSE)
  }

  class(model) <- "transfer_model"
  return(model)
}


#' Apply a Spectral Transfer Model
#'
#' Applies a pre-built transfer model to new source spectra to make them match
#' the target instrument's characteristics.
#'
#' @param new_source_spectra A numeric matrix of new spectra from the source
#'   instrument.
#' @param transfer_model A `transfer_model` object created by the
#'   `build_transfer_model` function.
#'
#' @return A matrix containing the transformed spectra.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'pds_model' was built earlier
#' corrected_new_spectra <- apply_transfer(new_isc_samples, pds_model)
#' }
apply_transfer <- function(new_source_spectra, transfer_model) {
  if (!inherits(transfer_model, "transfer_model")) {
    stop("The provided model is not a valid 'transfer_model' object.", call. = FALSE)
  }

  source_matrix <- as.matrix(new_source_spectra)

  if (transfer_model$method == "pds") {
    # Apply PDS: matrix multiplication and add intercept
    transferred <- source_matrix %*% transfer_model$transfer_matrix
    transferred <- sweep(transferred, 2, transfer_model$intercepts, "+")

  } else if (transfer_model$method == "subtraction") {
    # Apply Subtraction: add the correction vector to each row
    transferred <- sweep(source_matrix, 2, transfer_model$correction_vector, "+")
  }

  return(transferred)
}


# Internal helper function to provide a default value if an argument is NULL
# (equivalent to the future `??` operator)
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
