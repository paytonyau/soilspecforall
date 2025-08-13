# SoilSpec4All

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/paytonyau/SoilSpec4All/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paytonyau/SoilSpec4All/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `SoilSpec4All` R package provides a complete, feasible, and low-cost workflow for the environmentally sustainable analysis of soils. The project's objective was to develop a mathematical function for normalizing spectral data from a portable NIR spectrometer (ISC NIRvascan, 900-1700 nm) to match the quality of a full-range, research-grade instrument (ASD AgriSpec).

This allows users to leverage powerful predictive models built on high-quality reference data using only new data collected from the affordable portable sensor. The package implements a full workflow for predicting soil properties such as Total Organic Carbon (TOC), Total Nitrogen (TN), and pH, including:

* Data import and cleaning
* Spectral pre-processing (smoothing, resampling, SNV)
* Spectral transfer model building (PDS, Subtraction)
* Application of transfer models to new data
* Training and validation of predictive models

## Installation

You can install the development version of `SoilSpec4All` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("paytonyau/SoilSpec4All")
```

## A Complete Workflow Example

This example demonstrates the full power of the package, from pre-processing raw data to building a transfer model and making a final prediction.

```r
library(SoilSpec4All)
library(dplyr) # For the pipe operator %>%

# --- 1. Simulate Data ---
# Create realistic "paired" data for 70 samples, as if measured by both instruments.
set.seed(42)
n_samples <- 70
wavelengths <- seq(950, 1650, 5)
true_signal <- 0.5 - 0.4 * exp(-0.008 * (wavelengths - 1400)^2)

# Simulate high-quality ASD spectra
asd_spectra <- t(sapply(1:n_samples, function(i) true_signal + rnorm(length(wavelengths), 0, 0.005)))
colnames(asd_spectra) <- wavelengths

# Simulate low-cost ISC spectra (with more noise and a systematic offset)
isc_spectra <- t(sapply(1:n_samples, function(i) (true_signal + rnorm(length(wavelengths), 0, 0.01)) - 0.05))
colnames(isc_spectra) <- wavelengths

# Simulate corresponding soil carbon data
soil_data <- data.frame(tc_perc = 2.0 + (0.5 - asd_spectra[, "1400"]) * 8 + rnorm(n_samples, 0, 0.1))


# --- 2. Pre-process Spectra ---
# Apply the same pre-processing to both sets of spectra
asd_processed <- asd_spectra %>% smooth_sgolay() %>% transform_snv()
isc_processed <- isc_spectra %>% smooth_sgolay() %>% transform_snv()


# --- 3. Build the Transfer Model ---
# Use the paired, processed spectra to learn the transformation from ISC to ASD
transfer_model <- build_transfer_model(
  source_spectra = isc_processed,
  target_spectra = asd_processed,
  method = "pds"
)


# --- 4. Train the Predictive Model ---
# Train a "gold standard" model using only the high-quality ASD data
soil_model <- train_soil_model(
  spectra = asd_processed,
  response = soil_data$tc_perc,
  method = "cubist"
)


# --- 5. Use the Workflow on a New Sample ---
# Simulate a new sample measured ONLY with the low-cost ISC sensor
new_sample_isc <- matrix(rnorm(length(wavelengths), mean = 0.3, sd = 0.05), nrow = 1)
colnames(new_sample_isc) <- wavelengths

# Apply the same workflow: pre-process -> transfer -> predict
corrected_spectrum <- new_sample_isc %>%
  smooth_sgolay() %>%
  transform_snv() %>%
  apply_transfer(transfer_model)

# Re-apply column names lost during PDS matrix multiplication
colnames(corrected_spectrum) <- colnames(new_sample_isc)

# Get the final prediction
predicted_tc <- predict(soil_model, corrected_spectrum)

print(paste("Predicted Total Carbon (%) for new sample:", round(predicted_tc, 2)))
#> [1] "Predicted Total Carbon (%): 2.54"
```

## Team and Acknowledgements

This project was a cross-continental collaboration between Germany, Wales, Scotland, and the USA.

* **Principal Investigator:** Dr. Wanderson Mendes (Manaaki Whenua – Landcare Research, New Zealand)
* **Team Members:** Dr. Leandro Vieira-Filho (Clemson University, USA), Dr. Payton Yau (Scotland’s Rural College, Scotland), Dr. Kirst Elliott (Wardell Armstrong Inc., Wales).

We extend our gratitude to Prof. Dr. Michael Sommer (Leibniz Centre for Agricultural Landscape Research, Germany) and Prof. Dr. Maria Silveira (University of Florida, USA) for providing soil samples, and to Prof. Dr. Eyal Ben-Dor (Tel Aviv University, Israel) for providing a standard soil sample.

## Citation

If you use this package in your research, please cite the project and the associated publication (pending).

## Funding

This project was proudly supported by the **WCSS 2022 Early Career Interdisciplinary Grant** from the British Society of Soil Science.
