#' Split Data via cLHS
#'
#' This function splits a spectral dataset into calibration and validation sets
#' using the Conditioned Latin Hypercube Sampling (cLHS) method on the
#' principal components of the data.
#'
#' @param spectra A numeric matrix where rows are samples and columns are wavelengths.
#' @param split_ratio A numeric value for the proportion of data in the calibration set. Defaults to 0.7.
#' @param pc_selection A numeric value for the cumulative variance to be explained by PCA. Defaults to 0.999.
#' @return A list containing two integer vectors: `calibration_idx` and `validation_idx`.
#' @export
#' @importFrom resemble pc_projection
#' @importFrom clhs clhs
#' @examples
#' \dontrun{
#' # Create a dummy matrix
#' spec_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10)
#' colnames(spec_matrix) <- 1:10
#' # Get split indices
#' splits <- split_clhs(spec_matrix)
#' # Create calibration and validation sets
#' cal_set <- spec_matrix[splits$calibration_idx, ]
#' val_set <- spec_matrix[splits$validation_idx, ]
#' }
split_clhs <- function(spectra, split_ratio = 0.7, pc_selection = 0.999) {
  spectra_matrix <- as.matrix(spectra)

  # Perform PCA to reduce dimensionality for cLHS
  pc_proj <- resemble::pc_projection(
    Xr = spectra_matrix,
    pc_selection = list("cumvar", pc_selection),
    center = TRUE,
    scale = FALSE
  )

  # Perform cLHS on the PC scores
  clhs_result <- clhs::clhs(
    x = as.data.frame(pc_proj$scores),
    size = round(split_ratio * nrow(spectra_matrix)),
    iter = 1000,
    simple = FALSE,
    progress = FALSE
  )

  # Return a list of indices
  list(
    calibration_idx = clhs_result$index_samples,
    validation_idx = setdiff(1:nrow(spectra_matrix), clhs_result$index_samples)
  )
}


#' Train a Predictive Soil Model
#'
#' A wrapper around `caret::train` to build a model for predicting a soil
#' property from spectral data.
#'
#' @param spectra A numeric matrix or data frame of spectra (samples as rows).
#' @param response A numeric vector of the soil property to be predicted. Must
#'   have the same length as the number of rows in `spectra`.
#' @param method The modeling method to be used by `caret`. Defaults to "cubist".
#' @param train_control A `trainControl` object from the `caret` package. A
#'   default is provided if none is specified.
#'
#' @return A trained model object of class `train`.
#'
#' @importFrom caret train trainControl
#' @export
#'
#' @examples
#' \dontrun{
#' # Train a Cubist model using default repeated cross-validation
#' soil_model <- train_soil_model(cal_spectra, cal_soil_data$tc.perc)
#'
#' # Predict new values
#' predictions <- predict(soil_model, val_spectra)
#' }
train_soil_model <- function(spectra, response, method = "cubist", train_control = NULL) {
  if (is.null(train_control)) {
    # Provide a sensible default trainControl object if the user doesn't
    train_control <- caret::trainControl(
      method = 'repeatedcv',
      number = 10,
      repeats = 10,
      search = 'random'
    )
  }

  # Train the model using caret
  model <- caret::train(
    x = spectra,
    y = response,
    method = method,
    trControl = train_control,
    importance = TRUE
  )

  return(model)
}


#' Calculate Goodness-of-Fit Statistics
#'
#' Computes a comprehensive set of performance metrics for evaluating a model's
#' predictions against observed values.
#'
#' @param observed A numeric vector of the true, observed values.
#' @param predicted A numeric vector of the model's predicted values.
#' @param type The type of output required. Use `"spec"` for detailed
#'   spectrometric model evaluation or `"DSM"` for a simpler set.
#'
#' @return A data frame containing the calculated performance metrics.
#'
#' @importFrom stats lm sd quantile var
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate performance metrics
#' performance_metrics <- gof(observed = validation_data$tc.perc, predicted = predictions)
#' print(performance_metrics)
#' }
gof <- function(observed, predicted, type = "spec") {
  # Coefficient of determination (from a linear model)
  rLM <- stats::lm(predicted ~ observed)
  MEC <- as.matrix(summary(rLM)$adj.r.squared)

  # Mean Squared Error and Root Mean Squared Error
  MSE <- mean((observed - predicted)^2)
  RMSE <- sqrt(MSE)

  # Bias
  bias <- mean(predicted) - mean(observed)

  # Bias-corrected MSE and RMSE
  MSEc <- sum(((predicted - bias - observed)^2) / length(observed))
  RMSEc <- sqrt(MSEc)

  # Ratio of Performance to Deviation (RPD)
  RPD <- stats::sd(observed) / RMSE

  # Ratio of Performance to Interquartile distance (RPIQ)
  IQ <- diff(stats::quantile(observed, probs = c(0.25, 0.75)))
  RPIQ <- IQ / RMSE

  # Concordance Correlation Coefficient (CCC)
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- stats::var(observed)
  s2y <- stats::var(predicted)
  sxy <- mean((observed - mx) * (predicted - my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)

  if (type == "DSM") {
    gf <- data.frame(MEC = MEC, concordance = ccc, MSE = MSE, RMSE = RMSE, bias = bias, row.names = NULL)
  } else if (type == "spec") {
    gf <- data.frame(MEC = MEC, concordance = ccc, MSE = MSE, RMSE = RMSE, bias = bias,
                     MSEc = MSEc, RMSEc = RMSEc, RPD = RPD, RPIQ = RPIQ, row.names = NULL)
  } else {
    stop("ERROR: Revise the type of output you require. Select from either 'DSM' or 'spec'.", call. = FALSE)
  }

  return(gf)
}
