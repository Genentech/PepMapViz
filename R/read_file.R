#' Combine CSV and TXT Files from a Folder
#'
#' This function reads all CSV and TXT files from a specified folder and combines them
#' into a single data.table.
#'
#' @param folder_path The path to the folder containing the CSV or TSV files.
#' @return A data.table containing the combined data from all files.
#' @import data.table
#'
#' @examples
#' folder_path <- ""
#' combined_df <- combine_files_from_folder(folder_path)
#' print(combined_df)
#'
#' @export
combine_files_from_folder <- function(folder_path) {
  # List all relevant files (case insensitive)
  file_list <- list.files(folder_path,
                          pattern = "\\.csv$|\\.tsv$|\\.txt$|\\.mzid$|\\.mztab$",
                          full.names = TRUE,
                          ignore.case = TRUE)

  # Function to read files based on extension
  read_file <- function(file_path) {
    ext <- tolower(tools::file_ext(file_path))

    if (ext %in% c("csv", "tsv", "txt")) {
      return(data.table::fread(file_path))
    } else if (ext == "mzid") {
      if (!requireNamespace("mzID", quietly = TRUE))
        stop("Package 'mzID' required for mzID files. Install with BiocManager::install('mzID')")
      return(mzID::flatten(mzID::mzID(file_path)))
    } else if (ext == "mztab") {
      if (!requireNamespace("MSnbase", quietly = TRUE))
        stop("Package 'MSnbase' required for mzTab files. Install with BiocManager::install('MSnbase')")
      return(MSnbase::fData(MSnbase::readMzTabData(file_path, "PSM")))
    }
  }

  # Read and combine files with error handling
  df_list <- lapply(file_list, function(file) {
    tryCatch(
      read_file(file),
      error = function(e) {
        message("Error processing ", basename(file), ": ", e$message)
        NULL
      }
    )
  })

  # Remove failed reads and combine
  combined_df <- data.table::rbindlist(Filter(Negate(is.null), df_list), fill = TRUE)

  return(combined_df)
}
