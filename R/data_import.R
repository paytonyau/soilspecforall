#' Read and Combine ISC Spectra Files
#'
#' This function reads all raw CSV files from an ISC spectrometer located in a
#' specified directory. It combines them into a single, tidy data frame and
#' extracts a sample ID from each filename.
#'
#' @param path The full file path to the directory containing the raw spectra files.
#' @param pattern A character string containing a regular expression to identify
#'   the correct files. Defaults to "_r.csv".
#' @param id_regex A regular expression used to extract the sample ID from the
#'   filename. The default pattern `(?<=_)(\\d+)(?=\\_)` looks for digits
#'   enclosed by underscores.
#' @param id_prefix A character string to prepend to the extracted numeric ID.
#'   Defaults to "isc-".
#' @param skip_lines The number of header lines to skip at the top of each
#'   CSV file. Defaults to 28.
#'
#' @return A tidy tibble (a modern data frame) with three columns: `id`,
#'   `wavelength`, and `reflectance`.
#'
#' @importFrom readr read_csv
#' @importFrom purrr map_dfr
#' @importFrom dplyr mutate select rename
#' @importFrom stringr str_extract
#' @importFrom tools file_path_sans_ext
#' @importFrom rlang .data
#' @import magrittr
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define the path to the folder containing your raw spectral files
#' spec_dir <- "C:/BSSS_project_lowcost_nir/database/raw_data/spectra_usa_ge/isc/"
#'
#' # Run the function to import and clean the data
#' my_spectra <- read_isc_spectra(path = spec_dir)
#'
#' # Print the first few rows of the resulting data frame
#' print(head(my_spectra))
#' }
read_isc_spectra <- function(path, pattern = "_r.csv", id_regex = "(?<=_)(\\d+)(?=\\_)", id_prefix = "isc-", skip_lines = 28) {

  # Get the full list of file paths matching the pattern
  files <- list.files(path = path,
                      pattern = pattern,
                      full.names = TRUE,
                      recursive = TRUE)

  if (length(files) == 0) {
    stop("No files found in the specified path with the given pattern.", call. = FALSE)
  }

  combined_data <- purrr::map_dfr(
    files,
    ~ readr::read_csv(.x, skip = skip_lines, show_col_types = FALSE, progress = FALSE),
    .id = "source_file"
  )

  # Using the .data pronoun and explicit renaming to avoid CMD check NOTEs.
  processed_data <- combined_data %>%
    dplyr::mutate(
      id_extracted = stringr::str_extract(
        tools::file_path_sans_ext(basename(.data$source_file)),
        id_regex
      ),
      id = paste0(id_prefix, .data$id_extracted)
    ) %>%
    # Use explicit, quoted renaming for robustness
    dplyr::rename(
      wavelength = "Wavelength (nm)",
      reflectance = "Reflectance (R)"
    ) %>%
    dplyr::select(.data$id, .data$wavelength, .data$reflectance)

  return(processed_data)
}
