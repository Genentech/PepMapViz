#' Match peptide sequence with provided sequence and calculate positions
#'
#' This function matches peptide sequences from the 'peptide_data' data frame to
#' corresponding provided sequences in the 'whole_seq' data frame. It calculates
#' the start and end positions of the matched sequences and returns a data frame
#' with information about the matching positions.
#'
#' @param peptide_data A data frame containing peptide sequence information to match.
#' @param column The name of the column in peptide_data containing the peptide
#'              sequences to be matched.
#' @param whole_seq A data frame containing details about antibody sequence
#'              information including the domain and region information.
#'              'Region_Sequence' column is required for the sequence information.
#'              Change the column name if it is different than 'Region_Sequence'.
#' @param match_columns A character vector of column names that must exist in
#' both peptide_data and whole_seq. When searching for peptide sequence matches,
#' the function will only consider rows in whole_seq where the values in all
#' columns specified here exactly match the corresponding values in the current
#' row of peptide_data.
#' @param column_keep (Optional) The name of the columns in peptide_data to
#'                  keep in result data frame.
#' @param sequence_length (Optional) The sequence length range of peptide that
#'                we want to keep in the result. (e.g. c(1, 5) will include
#'                peptide sequence length from 1 to 5.)
#' @return A data frame with columns from 'peptide_data' and 'whole_seq'
#'        indicating the matched positions and related information.
#'
#' @importFrom stringr str_detect str_locate_all
#'
#' @examples
#' peptide_data <- data.frame(
#'   Sequence = c("AILNK", "BXLMR", "JJNXX", "DDEEF"),
#'   Condition_1 = c("Drug1", "Drug1", "Drug2", "Drug2"),
#'   Condition_2 = c("Donor1", "Donor2", "Donor1", "Donor2"),
#'   Region_1 = c("VH", "VL", "VH", "VL"),
#'   Region_2 = c("Arm_1", "Arm_2", "Arm_1", "Arm_2"),
#'   Area = c(100, 2, 4, NA)
#' )
#' whole_seq <- data.frame(
#'   Region_Sequence = c(
#'     "XYZAILNKPQR",
#'     "ABCBXLMRDEF",
#'     "GHIJJNXXKLM",
#'     "NOPDDEEFQRS",
#'     "AILXKPQR",
#'     "BNJLMRDEF",
#'     "ILNXXKLM",
#'     "DDEEXQRS",
#'     "XYZAAA",
#'     "XYZCCC",
#'     "XYZBBB",
#'     "XYZDDD",
#'     "XYZAAB",
#'     "XYZCCD",
#'     "XYZBBB",
#'     "XYZDDD"
#'   ),
#'   Condition_1 = c(
#'     "Drug1",
#'     "Drug1",
#'     "Drug2",
#'     "Drug2",
#'     "Drug1",
#'     "Drug1",
#'     "Drug2",
#'     "Drug2",
#'     "Drug1",
#'     "Drug1",
#'     "Drug2",
#'     "Drug2",
#'     "Drug1",
#'     "Drug1",
#'     "Drug2",
#'     "Drug2"
#'   ),
#'   Condition_2 = c(
#'     "Donor1",
#'     "Donor1",
#'     "Donor1",
#'     "Donor1",
#'     "Donor1",
#'     "Donor1",
#'     "Donor1",
#'     "Donor1",
#'     "Donor2",
#'     "Donor2",
#'     "Donor2",
#'     "Donor2",
#'     "Donor2",
#'     "Donor2",
#'     "Donor2",
#'     "Donor2"
#'   ),
#'   Region_1 = c(
#'     "VH",
#'     "VL",
#'     "VH",
#'     "VL",
#'     "VH",
#'     "VL",
#'     "VH",
#'     "VL",
#'     "VH",
#'     "VL",
#'     "VH",
#'     "VL",
#'     "VH",
#'     "VL",
#'     "VH",
#'     "VL"
#'   ),
#'   Region_2 = c(
#'     "Arm_1",
#'     "Arm_1",
#'     "Arm_1",
#'     "Arm_1",
#'     "Arm_2",
#'     "Arm_2",
#'     "Arm_2",
#'     "Arm_2",
#'     "Arm_1",
#'     "Arm_1",
#'     "Arm_1",
#'     "Arm_1",
#'     "Arm_2",
#'     "Arm_2",
#'     "Arm_2",
#'     "Arm_2"
#'   )
#' )
#' match_columns <- c("Condition_1", "Condition_2", "Region_1")
#' column_keep <- c("Region_2")
#' sequence_length <- c(1, 5)
#' column <- "Sequence"
#' matching_result <- match_and_calculate_positions(peptide_data,
#'                                                  column,
#'                                                  whole_seq,
#'                                                  match_columns,
#'                                                  sequence_length,
#'                                                  column_keep)
#'
#' @export
#'
match_and_calculate_positions <- function(peptide_data, column, whole_seq, match_columns, sequence_length = NULL, column_keep = NULL) {
  peptide_data <- as.data.frame(peptide_data)
  # capture the unquoted column name and then convert it to a character
  excluded_columns <- c(match_columns, "Region_Sequence")

  result <- list()

  for (i in 1:nrow(peptide_data)) {
    peptide <- as.character(peptide_data[[column]][i])
    peptide_regex <- convert_to_regex_pattern(peptide)
    matching_index <- integer(0)

    # Check if sequence_length is specified and the peptide length is within the range
    if (is.null(sequence_length) || (nchar(peptide) >= sequence_length[1] && nchar(peptide) <= sequence_length[2])) {
      matching_index <- which(apply(whole_seq[match_columns], 1, function(row) all(row == peptide_data[match_columns][i,])))
    }

    # Initialize empty result
    result_single_peptide <- data.frame()
    if (length(matching_index) > 0) {
      for (j in matching_index) {
        match_positions_list <- str_locate_all(as.character(whole_seq$Region_Sequence[j]), regex(peptide_regex))
        if (length(match_positions_list) <= 0) {
          next
        }
        # Loop over all match_positions
        for (match_positions in match_positions_list) {
          if (nrow(match_positions) == 0) {
            next
          }
          for (k in 1:nrow(match_positions)) {
            start_pos <- match_positions[k, 1]
            end_pos <- match_positions[k, 2]
            if (!is.na(start_pos) & !is.na(end_pos)) {
              remaining_columns <- setdiff(names(whole_seq), excluded_columns)
              peptide_row <- as.matrix(peptide_data[i, c(column, match_columns, column_keep), drop = FALSE])
              whole_seq_row <- as.matrix(whole_seq[j, remaining_columns, drop = FALSE])
              new_row <- cbind(peptide_row,
                               whole_seq_row,
                               start_position = start_pos,
                               end_position = end_pos)
              result_single_peptide <- rbind(result_single_peptide, new_row)
            }
          }
        }
      }
    }

    result[[i]] <- result_single_peptide
  }

  result <- do.call(rbind, result)
  return(result)
}

#' Convert Peptide Sequence to Regex Pattern
#'
#' This function converts a peptide sequence into a regular expression pattern
#' that accounts for ambiguous amino acids. Each amino acid is replaced by a
#' character class that includes itself, 'X', and any specific ambiguities.
#'
#' @param peptide A character string representing the peptide sequence.
#' @return A character string containing the regex pattern for matching.
#' @examples
#' # Convert a peptide sequence to a regex pattern
#' peptide <- "NDEQIL"
#' regex_pattern <- convert_to_regex_pattern(peptide)
#' print(regex_pattern) # Output: "[NBX][DBX][EZX][QZX][ILX][ILX]"
#' @export
#'
convert_to_regex_pattern <- function(peptide) {
  # Create a list of amino acid substitutions
  substitutions <- list(
    A = "[AX]", C = "[CX]", D = "[BDX]", E = "[EZX]", F = "[FX]",
    G = "[GX]", H = "[HX]", I = "[ILJX]", K = "[KX]", L = "[ILJX]",
    M = "[MX]", N = "[BNX]", P = "[PX]", Q = "[QZX]", R = "[RX]",
    S = "[SX]", T = "[TX]", V = "[VX]", W = "[WX]", Y = "[YX]",
    B = "[BDNX]", Z = "[ZEQX]", J = "[ILJX]", O = "[OX]", U = "[UX]"
  )

  # Split the peptide into individual characters
  peptide_chars <- strsplit(peptide, "")[[1]]

  # Replace each character with its regex equivalent
  regex_chars <- sapply(peptide_chars, function(char) {
    if (char %in% names(substitutions)) {
      substitutions[[char]]
    } else {
      paste0("[", char, "X]")  # For any unspecified character, match itself or X
    }
  })

  # Join the characters back into a single string
  regex_pattern <- paste(regex_chars, collapse = "")

  return(regex_pattern)
}
