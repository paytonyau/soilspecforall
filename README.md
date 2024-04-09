# SoilSpec4All

## Introduction
SoilSpec4All is an innovative project aimed at enhancing soil testing services by utilising affordable, portable, and miniaturised NIR sensors. Our mission is to provide efficient and cost-effective soil property analysis, contributing to sustainable agriculture and natural resource management.

## Project Overview
We focus on normalising spectral data from the NIRvascan Portable Smart Spectrometer (https://www.alliedscientificpro.com/nirvascan) to match the quality of the research-grade ASD - AgriSpec sensor. This normalisation process is crucial for predicting soil properties such as carbon, nitrogen, and phosphorus levels across various soil types.

## Contents
- **Data Import & Cleaning**: [`1A) Importing and claning.R`](https://github.com/paytonyau/soilspecforall/blob/main/01_isc_importing_and_cleaning.R) | [`1B) Rmarkdown`](https://github.com/paytonyau/soilspecforall/blob/main/01_isc_importing_and_cleaning.Rmd)
- **Spectral Data Transfer**: [`2A) transfer_asd_to_isc_nir.R`](https://github.com/paytonyau/soilspecforall/blob/main/02_transfer_asd_to_isc_nir.R) | [`2B) Rmarkdown`](https://github.com/paytonyau/soilspecforall/blob/main/02_transfer_asd_to_isc_nir.Rmd)
- **Low-Cost NIR Modelling**: [`3A) proj_low_cost_nir_modelling.R`](https://github.com/paytonyau/soilspecforall/blob/main/03_proj_low_cost_nir_modelling.R) | [`3B) Rmarkdown`](https://github.com/paytonyau/soilspecforall/blob/main/03_proj_low_cost_nir_modelling.Rmd)
- **Free ASD Library Modelling**: [`4A) modelling_all_free_asd_library.R`](https://github.com/paytonyau/soilspecforall/blob/main/04_modelling_all_free_asd_library.R) | [`4B) Rmarkdown`](https://github.com/paytonyau/soilspecforall/blob/main/04_modelling_all_free_asd_library.Rmd)
- **Visualisation**: [`5) Quick graphs.R`](https://github.com/paytonyau/soilspecforall/blob/main/05_quick_graphs.R)

## References
- [Preparing a soil spectral library using the Internal Soil Standard (ISS) method: Influence of extreme different humidity laboratory conditions](https://www.sciencedirect.com/science/article/pii/S0016706118323619)
- [Reflectance measurements of soils in the laboratory: Standards and protocols](https://www.sciencedirect.com/science/article/pii/S0016706115000038)
- [Evaluating low-cost portable near infrared sensors for rapid analysis of soils from South Eastern Australia](https://www.sciencedirect.com/science/article/pii/S2352009419302391)

## Support
This project is proudly supported by the WCSS 2022 INTERDISCIPLINARY GRANT from the British Society of Soil Science.
