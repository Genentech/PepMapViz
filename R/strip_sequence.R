#' Strip peptide sequences based on the specified data type
#'
#' This function takes outputs from multiple platform, a data frame with a
#' column containing peptide sequences to be stripped,
#' and a column where the stripped sequences will be stored. The function chooses
#' the appropriate stripping method based on the specified data type ('PEAKS',
#' 'Spectronaut', 'MSFragger', 'Comet', 'DIANN', 'Skyline' or 'Maxquant').
#'
#' @param data A data frame with the peptide sequences.
#' @param column The name of the column containing the peptide sequences to be stripped.
#' @param convert_column The name of the column where the stripped sequences will be stored.
#' @param type A character string specifying the data type (e.g. 'Skyline' or 'Maxquant').
#'
#' @return A data frame with the specified column containing stripped sequences.
#'
#' @examples
#' library(data.table)
#' data_skyline <- data.table(
#'   'Peptide Modified Sequence' = c(
#'     "AGLC[+57]QTFVYGGC[+57]R",
#'     "AAAASAAEAGIATTGTEDSDDALLK",
#'     "IVGGWEC[+57]EK"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' data_maxquant <- data.table(
#'   'Modified sequence' = c(
#'     "_(ac)AAAAELRLLEK_",
#'     "_EAAENSLVAYK_",
#'     "_AADTIGYPVM(ox)IRSAYALGGLGSGICPNK_"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#'
#' converted_data_skyline <- strip_sequence(data_skyline,
#'                                          'Peptide Modified Sequence',
#'                                          'Sequence',
#'                                          "Skyline")
#' converted_data_maxquant <- strip_sequence(data_maxquant, 'Modified sequence',
#'                                           'Sequence', "Maxquant")
#'
#' @import data.table
#'
#' @export
strip_sequence <- function(data, column, convert_column, type) {
  if (type == "PEAKS") {
    data <- strip_sequence_PEAKS(data, column, convert_column)
  } else if (type == "Spectronaut") {
    data <- strip_sequence_Spectronaut(data, column, convert_column)
  } else if (type == "MSFragger") {
    data <- strip_sequence_MSFragger(data, column, convert_column)
  } else if (type == "Comet") {
    data <- strip_sequence_Comet(data, column, convert_column)
  } else if (type == "DIANN") {
    data <- strip_sequence_Maxquant(data, column, convert_column)
  } else if (type == "Skyline") {
    data <- strip_sequence_Skyline(data, column, convert_column)
  } else if (type == "Maxquant") {
    data <- strip_sequence_Maxquant(data, column, convert_column)
  } else {
    stop(
      "Invalid type. Supported types are 'PEAKS', 'Spectronaut', 'MSFragger',
         'Comet', 'DIANN', 'Skyline' and 'Maxquant'."
    )
  }
  return(data)
}
#' Strip sequence from PEAKS outputs
#'
#' This function takes PEAKS output containing a column with peptide sequences to be stripped
#' and converts it into a new dataframe with the stripped sequence
#'
#' @param data A dataframe with a column containing peptide sequences to be stripped
#' @param column The name of the column containing the peptide sequences to be stripped.
#' @param convert_column The name of the column where the stripped sequences will be stored.
#' @return A dataframe with a column containing stripped sequence
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   Peptide = c(
#'     "AAN(+0.98)Q(-0.98)RGSLYQCDYSTGSC(+57.02)EPIR",
#'     "K.AAQQTGKLVHANFGT.K",
#'     "K.(+0.98)AATVTGKLVHANFGT.K"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' column <- "Peptide"
#' convert_column <- "Sequence"
#' converted_data <- strip_sequence_PEAKS(data, column, convert_column)
#'
#' @import data.table
#'
#' @export
strip_sequence_PEAKS <- function(data, column, convert_column) {
  # Remove characters before the first dot and after the last dot
  data[, (convert_column) := gsub("^[A-Za-z]\\.|\\.[A-Za-z]$", "", get(column))]

  # Remove characters within parentheses
  data[, (convert_column) := gsub("\\([^)]+\\)", "", get(convert_column))]

  return(data)
}

#' Strip sequence from Spectronaut outputs
#'
#' This function takes Spectronaut output containing a column with peptide sequences to be stripped
#' and converts it into a new dataframe with the stripped sequence
#'
#' @param data A dataframe with a column containing peptide sequences to be stripped
#' @param column The name of the column containing the peptide sequences to be stripped.
#' @param convert_column The name of the column where the stripped sequences will be stored.
#' @return A dataframe with a column containing stripped sequence
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   EG.IntPIMID = c(
#'     "_[+42]M[-16]DDREDLVYQAK_",
#'     "_EAAENSLVAYK_",
#'     "_IEAELQDIC[+57]NDVLELLDK_"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' converted_data <- strip_sequence_Spectronaut(data, 'EG.IntPIMID', 'Sequence')
#'
#' @import data.table
#'
#' @export
strip_sequence_Spectronaut <- function(data, column, convert_column) {
  data[, (convert_column) := gsub("^_+|_+$|\\[.*?\\]", "", get(column))]
  return(data)
}

#' Strip sequence from MSFragger outputs
#'
#' This function takes MSFragger output containing a column with peptide sequences to be stripped
#' and converts it into a new dataframe with the stripped sequence
#'
#' @param data A dataframe with a column containing peptide sequences to be stripped
#' @param column The name of the column containing the peptide sequences to be stripped.
#' @param convert_column The name of the column where the stripped sequences will be stored.
#' @return A dataframe with a column containing stripped sequence
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   'Modified Peptide' = c(
#'     "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
#'     "K.AAQQTGKLVHANFGT.K",
#'     "K.[0.98]AATVTGKLVHANFGT.K"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' column <- 'Modified Peptide'
#' convert_column <- 'Sequence'
#' converted_data <- strip_sequence_MSFragger(data, 'Modified Peptide', 'Sequence')
#'
#' @import data.table
#'
#' @export
strip_sequence_MSFragger <- function(data, column, convert_column) {
  # Remove characters before the first dot and after the last dot
  data[, (convert_column) := gsub("^[A-Za-z]\\.|\\.[A-Za-z]$", "", get(column))]

  # Remove characters within square brackets
  data[, (convert_column) := gsub("^_+|_+$|\\[.*?\\]", "", get(convert_column))]

  return(data)
}

#' Strip sequence from Comet outputs
#'
#' This function takes Comet output containing a column with peptide sequences to be stripped
#' and converts it into a new dataframe with the stripped sequence
#'
#' @param data A dataframe with a column containing peptide sequences to be stripped
#' @param column The name of the column containing the peptide sequences to be stripped.
#' @param convert_column The name of the column where the stripped sequences will be stored.
#' @return A dataframe with a column containing stripped sequence
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   modified_peptide = c(
#'     "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
#'     "K.AAQQTGKLVHANFGT.K",
#'     "K.[0.98]AATVTGKLVHANFGT.K"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' column <- 'modified_peptide'
#' convert_column <- 'Sequence'
#' converted_data <- strip_sequence_Comet(data, column, convert_column)
#'
#' @import data.table
#'
#' @export
strip_sequence_Comet <- function(data, column, convert_column) {
  # Remove characters before the first dot and after the last dot
  data[, (convert_column) := gsub("^[A-Za-z]\\.|\\.[A-Za-z]$", "", get(column))]

  # Remove characters within square brackets
  data[, (convert_column) := gsub("^_+|_+$|\\[.*?\\]", "", get(convert_column))]
  return(data)
}

#' Strip sequence from DIANN outputs
#'
#' This function takes DIANN output containing a column with peptide sequences to be stripped
#' and converts it into a new dataframe with the stripped sequence
#'
#' @param data A dataframe with a column containing peptide sequences to be stripped
#' @param column The name of the column containing the peptide sequences to be stripped.
#' @param convert_column The name of the column where the stripped sequences will be stored.
#' @return A dataframe with a column containing stripped sequence
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   Modified.Sequence = c(
#'     "AAAAGPGAALS(UniMod:21)PRPC(UniMod:4)DSDPATPGAQSPK",
#'     "AAAASAAEAGIATTGTEDSDDALLK",
#'     "AAAAALSGSPPQTEKPT(UniMod:21)HYR"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' column <- 'Modified.Sequence'
#' convert_column <- 'Sequence'
#' converted_data <- strip_sequence_DIANN(data, column, convert_column)
#'
#' @import data.table
#'
#' @export
strip_sequence_DIANN <- function(data, column, convert_column) {
  data[, (convert_column) := gsub("\\([^)]+\\)", "", get(column))]
  return(data)
}

#' Strip sequence from Skyline outputs
#'
#' This function takes Skyline output containing a column with peptide sequences to be stripped
#' and converts it into a new dataframe with the stripped sequence
#'
#' @param data A dataframe with a column containing peptide sequences to be stripped
#' @param column The name of the column containing the peptide sequences to be stripped.
#' @param convert_column The name of the column where the stripped sequences will be stored.
#' @return A dataframe with a column containing stripped sequence
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   'Peptide Modified Sequence' = c(
#'     "AGLC[+57]QTFVYGGC[+57]R",
#'     "AAAASAAEAGIATTGTEDSDDALLK",
#'     "IVGGWEC[+57]EK"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' column <- 'Peptide Modified Sequence'
#' convert_column <- 'Sequence'
#' converted_data <- strip_sequence_Skyline(data, column, convert_column)
#'
#' @import data.table
#'
#' @export
strip_sequence_Skyline <- function(data, column, convert_column) {
  data[, (convert_column) := gsub("^_+|_+$|\\[.*?\\]", "", get(column))]
  return(data)
}

#' Strip sequence from Maxquant outputs
#'
#' This function takes Maxquant output containing a column with peptide sequences to be stripped
#' and converts it into a new dataframe with the stripped sequence
#'
#' @param data A dataframe with a column containing peptide sequences to be stripped
#' @param column The name of the column containing the peptide sequences to be stripped.
#' @param convert_column The name of the column where the stripped sequences will be stored.
#' @return A dataframe with a column containing stripped sequence
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   'Modified sequence' = c(
#'     "_(ac)AA(ox)AAELRLLEK_",
#'     "_EAAENSLVAYK_",
#'     "_AADTIGYPVM(ox)IRSAYALGGLGSGICPNK_"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' column <- 'Modified sequence'
#' convert_column <- 'Sequence'
#' converted_data <- strip_sequence_Maxquant(data, column, convert_column)
#'
#' @import data.table
#'
#' @export
strip_sequence_Maxquant <- function(data, column, convert_column) {
  # Remove leading and trailing underscores and the characters within square brackets
  data[, (convert_column) := gsub("^_+|_+$|\\(.*?\\)", "", get(column))]
}
