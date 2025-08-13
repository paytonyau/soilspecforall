#' Apply Splice Correction to ASD Spectra
#'
#' Corrects for instrumental artifacts at detector joins in ASD spectra. This
#' function is a wrapper around `prospectr::spliceCorrection`.
#'
#' @param spectra A numeric matrix or data frame of spectra (samples as rows,
#'   wavelengths as columns).
#' @param wavelengths A numeric vector of the corresponding wavelengths.
#' @param splice_points A numeric vector indicating the wavelengths where
#'   detectors are joined (e.g., `c(1000, 1800)`).
#' @param interpol_bands The number of bands to use for interpolation around the
#'   splice points. Defaults to 10.
#'
#' @return A matrix with the splice-corrected spectra.
#'
#' @importFrom prospectr spliceCorrection
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'asd_spectra' is a matrix and 'wavs' is a vector of wavelengths
#' corrected_spectra <- correct_splice(asd_spectra, wavs, splice_points = c(1000, 1800))
#' }
correct_splice <- function(spectra, wavelengths, splice_points, interpol_bands = 10) {
  # Ensure input is a matrix for prospectr
  spectra_matrix <- as.matrix(spectra)

  corrected <- prospectr::spliceCorrection(
    X = spectra_matrix,
    wav = wavelengths,
    splice = splice_points,
    interpol.bands = interpol_bands
  )
  return(corrected)
}


#' Apply Savitzky-Golay Smoothing or Derivatives
#'
#' A wrapper around `prospectr::savitzkyGolay` to smooth spectra or compute
#' derivatives.
#'
#' @param spectra A numeric matrix or data frame of spectra.
#' @param derivative_order The order of the derivative to compute (m). Use 0 for
#'   smoothing only. Defaults to 0.
#' @param polynomial_order The order of the polynomial to fit (p). Defaults to 2.
#' @param window_size The size of the moving window (w). Must be an odd number.
#'   Defaults to 11.
#'
#' @return A matrix with the transformed spectra.
#'
#' @importFrom prospectr savitzkyGolay
#' @export
#'
#' @examples
#' \dontrun{
#' # Smooth the spectra
#' smoothed_spectra <- smooth_sgolay(my_spectra, window_size = 11)
#'
#' # Calculate the first derivative
#' first_derivative <- smooth_sgolay(my_spectra, derivative_order = 1)
#' }
smooth_sgolay <- function(spectra, derivative_order = 0, polynomial_order = 2, window_size = 11) {
  spectra_matrix <- as.matrix(spectra)

  smoothed <- prospectr::savitzkyGolay(
    X = spectra_matrix,
    m = derivative_order,
    p = polynomial_order,
    w = window_size
  )
  return(smoothed)
}


#' Resample Spectra to a New Wavelength Resolution
#'
#' A wrapper around `prospectr::resample` to standardize spectral resolution
#' using spline interpolation.
#'
#' @param spectra A numeric matrix or data frame of spectra.
#' @param current_wavelengths A numeric vector of the original wavelengths.
#' @param new_wavelengths A numeric vector of the desired new wavelengths.
#'
#' @return A matrix with the resampled spectra. Column names are set to the
#'   new wavelengths.
#'
#' @importFrom prospectr resample
#' @export
#'
#' @examples
#' \dontrun{
#' # Resample to a 5 nm interval between 950 and 1650 nm
#' new_wavs <- seq(950, 1650, by = 5)
#' resampled_spectra <- resample_spectra(my_spectra, old_wavs, new_wavs)
#' }
resample_spectra <- function(spectra, current_wavelengths, new_wavelengths) {
  spectra_matrix <- as.matrix(spectra)

  resampled <- prospectr::resample(
    X = spectra_matrix,
    wav = current_wavelengths,
    new.wav = new_wavelengths,
    interpol = "spline"
  )
  # Set column names for clarity
  colnames(resampled) <- new_wavelengths
  return(resampled)
}


#' Apply Standard Normal Variate (SNV) Transformation
#'
#' A wrapper around `prospectr::standardNormalVariate` to correct for
#' multiplicative scatter effects.
#'
#' @param spectra A numeric matrix or data frame of spectra.
#'
#' @return A matrix with the SNV-transformed spectra.
#'
#' @importFrom prospectr standardNormalVariate
#' @export
#'
#' @examples
#' \dontrun{
#' # Apply SNV after smoothing
#' snv_spectra <- my_spectra %>%
#'   smooth_sgolay() %>%
#'   transform_snv()
#' }
transform_snv <- function(spectra) {
  spectra_matrix <- as.matrix(spectra)
  snv_transformed <- prospectr::standardNormalVariate(X = spectra_matrix)
  return(snv_transformed)
}
