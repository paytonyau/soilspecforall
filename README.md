# SoilSpec4All

Over the past few decades, soil spectral libraries (SSLs) have been created on scales ranging from local farms to the entire globe using visible-near-infrared (Vis-NIR) sensors. As the availability of affordable, portable, and miniaturised NIR sensors with specific wavelength ranges increases, we are keen to investigate their effectiveness in estimating soil properties.

However, due to variations in technical specifications such as spectral resolution, range, and energy intensity, different NIR sensors may yield varying performances. This necessitates the normalisation of the spectral information retrieved from different sensors, enabling end-users and researchers to expand soil testing services within a more manageable budget.

In this context, we are working on a mathematical function to normalise the spectral data obtained from the cost-effective NIRvascan Portable Smart Spectrometer (https://www.alliedscientificpro.com/nirvascan). The aim is to make it comparable to the data from the research-grade spectral ASD - AgriSpec sensor. Our goal is to build a reliable model that can predict the levels of carbon, nitrogen, and phosphorus in any type of soil using a miniaturised NIR sensor.

This approach will facilitate more efficient, regular, and current soil information on soil properties at affordable costs. It will aid in the implementation of policies for sustainable agriculture and natural resource management.

**Table of contents**

(1A) [Importing and claning.R](https://github.com/paytonyau/soilspecforall/blob/main/01_isc_importing_and_cleaning.R)
(1B) [Importing and claning - Rmarkdown](https://github.com/paytonyau/soilspecforall/blob/main/01_isc_importing_and_cleaning.Rmd)
(2A) [transfer_asd_to_isc_nir.R](https://github.com/paytonyau/soilspecforall/blob/main/02_transfer_asd_to_isc_nir.R)
(2B)[transfer_asd_to_isc_nir - Rmarkdown](https://github.com/paytonyau/soilspecforall/blob/main/02_transfer_asd_to_isc_nir.Rmd)
(3A)[proj_low_cost_nir_modelling.R](https://github.com/paytonyau/soilspecforall/blob/main/03_proj_low_cost_nir_modelling.R)
(3B)[proj_low_cost_nir_modelling - Rmarkdown](https://github.com/paytonyau/soilspecforall/blob/main/03_proj_low_cost_nir_modelling.Rmd)
(4A)[modelling_all_free_asd_library.R](https://github.com/paytonyau/soilspecforall/blob/main/04_modelling_all_free_asd_library.R)
(4B)[modelling_all_free_asd_library - Rmarkdown](https://github.com/paytonyau/soilspecforall/blob/main/04_modelling_all_free_asd_library.Rmd)
(5)[Quick graphs.R](https://github.com/paytonyau/soilspecforall/blob/main/05_quick_graphs.R)


**References**
(1) [Preparing a soil spectral library using the Internal Soil Standard (ISS) method: Influence of extreme different humidity laboratory conditions - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0016706118323619)
(2) [Reflectance measurements of soils in the laboratory: Standards and protocols - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0016706115000038)
(3) [Evaluating low-cost portable near infrared sensors for rapid analysis of soils from South Eastern Australia - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2352009419302391)

The package is under development

This project is supported by **WCSS 2022 INTERDISCIPLINARY GRANT** from the British Society of Soil Science.