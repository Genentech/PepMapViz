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
  # List all CSV and TXT files in the folder
  file_list <- list.files(folder_path, pattern = "\\.csv$|\\.tsv$", full.names = TRUE)

  # Read and combine files using data.table's fread and rbindlist functions
  combined_df <- rbindlist(lapply(file_list, fread))

  return(combined_df)
}
