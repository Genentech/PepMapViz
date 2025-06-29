#' Obtain post translational modification(PTM) information from Peptide data
#' based on the specified data type
#'
#' This function takes outputs from multiple platform, a data frame with column
#' containing modified peptide sequence with the detailed post translational
#' modification(PTM) information and converts it into a new dataframe with the desired format of peptide
#' sequences and associated PTM information. Due to the flexibility of outputs from
#' multiple platform, the PTM mass to type table needs to be provided if convertion to PTM_type is needed.
#' The result includes 'Peptide', 'PTM_position', 'PTM_type' and 'PTM_mass' columns.The function chooses
#' the appropriate converting method based on the specified data type ('PEAKS',
#' 'Spectronaut', 'MSFragger', 'Comet', 'DIANN', 'Skyline', 'Maxquant', 'mzIdenML' or 'mzTab'),
#' allowing you to convert the data into a consistent format for further analysis.
#'
#' @param data A data frame with the peptide sequences.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param seq_column (Optional) The name of the column containing peptide sequences for MSFragger, mzid and mzTab.
#'                     This parameter is required for the "MSFragger", "mzIdenML" and "mzTab" type and can be omitted
#'                    for other types.
#' @param PTM_table A data frame with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param type A character string specifying the data type (e.g. 'Skyline' or 'Maxquant').
#' @param PTM_mass_column The name of the column containing the PTM mass information.
#'
#' @return A data.table with 'PTM_position', 'PTM_type', 'PTM_mass', 'reps', and other columns.
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
#' PTM_table <- data.table(
#'   PTM_mass = c(57.02, -0.98, 15.9949),
#'   PTM_type = c("Cam", "Amid", "Ox")
#' )
#' converted_data_skyline <- obtain_mod(
#'   data_skyline,
#'   'Peptide Modified Sequence',
#'   'Skyline',
#'   seq_column = NULL,
#'   PTM_table,
#'   PTM_annotation = TRUE,
#'   PTM_mass_column = "PTM_mass"
#' )
#'
#' data_maxquant <- data.table(
#'   'Modified sequence' = c(
#'     "_(ac)AAAAELRLLEK_",
#'     "_EAAENSLVAYK_",
#'     "_AADTIGYPVM(ox)IRSAYALGGLGSGICPNK_"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' PTM_table <- data.table(
#'   PTM_mass = c('Phospho (STY)', 'Oxidation (M)'),
#'   PTM_type = c("Phos", "Ox")
#' )
#' converted_data_maxquant <- obtain_mod(
#'   data_maxquant,
#'   'Modified sequence',
#'   'Maxquant',
#'   seq_column = NULL,
#'   PTM_table,
#'   PTM_annotation = TRUE,
#'   PTM_mass_column = "PTM_mass"
#' )
#'
#'
#' @import data.table
#' @import stringr
#'
#' @export

# Define the wrap-up function
obtain_mod <- function(data,
                       mod_column,
                       type,
                       seq_column = NULL,
                       PTM_table = NULL,
                       PTM_annotation = FALSE,
                       PTM_mass_column) {
  if (type == "MSFragger") {
    if (is.null(seq_column)) {
      stop("seq_column is required for 'MSFragger'.")
    }
    result <- obtain_mod_MSFragger(data,
                                   mod_column,
                                   seq_column,
                                   PTM_table,
                                   PTM_annotation,
                                   PTM_mass_column)
  } else if (type == "mzIdenML") {
      if (is.null(seq_column)) {
        stop("seq_column is required for 'mzIdenML'.")
      }
      result <- obtain_mod_mzIdenML(data,
                                     mod_column,
                                     seq_column,
                                     PTM_table,
                                     PTM_annotation,
                                     PTM_mass_column)
  } else if (type == "mzTab") {
    if (is.null(seq_column)) {
      stop("seq_column is required for 'mzTab'.")
    }
    result <- obtain_mod_mzTab(data,
                                   mod_column,
                                   seq_column,
                                   PTM_table,
                                   PTM_annotation,
                                   PTM_mass_column)
  } else if (type == "Spectronaut") {
    result <- obtain_mod_Spectronaut(data, mod_column, PTM_table, PTM_annotation, PTM_mass_column)
  } else if (type == "PEAKS") {
    result <- obtain_mod_PEAKS(data, mod_column, PTM_table, PTM_annotation, PTM_mass_column)
  } else if (type == "Comet") {
    result <- obtain_mod_Comet(data, mod_column, PTM_table, PTM_annotation, PTM_mass_column)
  } else if (type == "DIANN") {
    result <- obtain_mod_DIANN(data, mod_column, PTM_table, PTM_annotation, PTM_mass_column)
  } else if (type == "Skyline") {
    result <- obtain_mod_Skyline(data, mod_column, PTM_table, PTM_annotation, PTM_mass_column)
  } else if (type == "Maxquant") {
    result <- obtain_mod_Maxquant(data, mod_column, PTM_table, PTM_annotation, PTM_mass_column)
  } else {
    stop(
      "Invalid type. Supported types are 'PEAKS', 'Spectronaut', 'MSFragger',
         'Comet', 'DIANN', 'Skyline', 'Maxquant', 'mzIdenML' or 'mzTab'"
    )
  }
  return(result)
}

#' Obtain modification information from Peptide data generated by PEAKS
#'
#' This function takes PEAKS output containing a column with modified peptide
#' sequences including PTM information and converts it into a new dataframe with the
#' desired format of peptide sequences and associated PTM information.
#'
#' @param data A dataframe with a column containing modified peptide sequences.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param PTM_table A dataframe with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param PTM_mass_column The name of the column containing the PTM mass information
#' @return A data.table with 'PTM_position', 'PTM_type', 'PTM_mass', 'reps', and other columns.
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   Peptide = c(
#'     "AAN(+42)Q(-0.98)RGSLYQCDYSTGSC(+57.02)EPIR",
#'     "K.AAQQTGKLVHANFGT.K",
#'     "K.(-0.98)AATVTGKLVHANFGT.K"
#'   ),
#'   Sequence = c(
#'     "AANQRGSLYQCDYSTGSCEPIR",
#'     "AAQQTGKLVHANFGT",
#'     "AATVTGKLVHANFGT"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' PTM_table <- data.table(PTM_mass = c(42, -0.98, 57.02),
#'                         PTM_type = c("Acet", "Amid", "Cam"))
#' mod_column <- "Peptide"
#' PTM_mass_column <- "PTM_mass"
#' converted_data <- obtain_mod_PEAKS(data, mod_column, PTM_table,
#' PTM_annotation = TRUE, PTM_mass_column)
#'
#' @import data.table
#'
#' @export
#'
obtain_mod_PEAKS <- function(data,
                             mod_column,
                             PTM_table = NULL,
                             PTM_annotation = FALSE,
                             PTM_mass_column) {
  # Ensure data is a data.table
  PTM_table <- as.data.table(PTM_table)

  # Extract PTM information using regular expression
  ptm_matches <- gregexpr("\\(([^)]+)\\)", data[[mod_column]])
  ptm_info <- regmatches(data[[mod_column]], ptm_matches)
  for (i in seq_along(ptm_info)) {
    ptm_info[[i]] <- gsub("\\+|\\(|\\)", "", ptm_info[[i]])
  }

  # Calculate PTM positions for each row
  ptm_positions <- lapply(seq_along(ptm_matches), function(i) {
    ptm_lengths <- attr(ptm_matches[[i]], "match.length")
    if (-1 %in% ptm_lengths) {
      ptm_positions <- NA

    } else {
      ptm_values <- unlist(ptm_matches[i])
      ptm_lengths <- attr(ptm_matches[[i]], "match.length")
      ptm_positions <- ptm_values - cumsum(c(1, ptm_lengths[-length(ptm_lengths)]))
    }
  })

  # Combine results for each row
  rep_values <- ifelse(sapply(ptm_info, length) == 0, 1, sapply(ptm_info, length))
  result <- data.table(
    PTM_position = unlist(ptm_positions),
    reps = rep(rep_values, rep_values)
  )

  # Dynamically add the PTM_mass column
  result[, (PTM_mass_column) := unlist(lapply(ptm_info, function(x)
    if (length(x) > 0)
      x
    else
      NA))]

  for (col in colnames(data)) {
    result[, (col) := unlist(lapply(seq_along(ptm_info), function(i)
      if (length(ptm_info[[i]]) > 0)
        rep(data[[col]][i], length(ptm_info[[i]]))
      else
        data[[col]][i]))]
  }

  if (PTM_annotation & !is.null(PTM_table)) {
    PTM_table[, (PTM_mass_column) := as.character(get(PTM_mass_column))]
    result <- merge(result, PTM_table, by = PTM_mass_column, all.x = TRUE)
  }

  return(result)
}

#' Obtain modification information from Peptide data generated by Spectronaut
#'
#' This function takes Spectronaut output containing a column with modified peptide sequences
#' including PTM information and converts it into a new dataframe with the desired format of peptide
#' sequences and associated PTM information.
#'
#' @param data A data.table with a column containing modified peptide sequences.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param PTM_table A data.table with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param PTM_mass_column The name of the column containing the PTM mass information
#' @return A data.table with 'PTM_position', 'PTM_type', 'reps', and other columns.
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   EG.ModifiedPeptide = c(
#'     "_[Acetyl (Protein N-term)]M[Oxidation (M)]DDREDLVYQAK_",
#'     "_EAAENSLVAYK_",
#'     "_IEAELQDIC[Carbamidomethyl (C)]NDVLELLDK_"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' PTM_table <- data.table(
#'   PTM_mass = c(
#'     'Acetyl (Protein N-term)',
#'     'Oxidation (M)',
#'     'Carbamidomethyl (C)'
#'   ),
#'   PTM_type = c("Acet", "Ox", "Cam")
#' )
#' converted_data <- obtain_mod_Spectronaut(data, 'EG.ModifiedPeptide',
#'                                          PTM_table, PTM_annotation = TRUE,
#'                                          PTM_mass_column = "PTM_mass")
#' data <- data.table(
#'   EG.IntPIMID = c(
#'     "_[+42]M[-0.98]DDREDLVYQAK_",
#'     "_EAAENSLVAYK_",
#'     "_IEAELQDIC[+57]NDVLELLDK_"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' PTM_table <- data.table(PTM_mass = c(42, -0.98, 57),
#'                         PTM_type = c("Acet", "Amid", "Cam"))
#' PTM_mass_column <- "PTM_mass"
#' converted_data <- obtain_mod_Spectronaut(data,
#'                                          'EG.IntPIMID',
#'                                          PTM_table,
#'                                          PTM_annotation = TRUE,
#'                                          PTM_mass_column)
#'
#' @import data.table
#'
#' @export


obtain_mod_Spectronaut <- function(data,
                                   mod_column,
                                   PTM_table = NULL,
                                   PTM_annotation = FALSE,
                                   PTM_mass_column) {
  # Ensure data is a data.table
  PTM_table <- as.data.table(PTM_table)

  # Extract PTM information using regular expression
  ptm_matches <- gregexpr("\\[([^]]+)\\]", data[[mod_column]])
  ptm_info <- regmatches(data[[mod_column]], ptm_matches)
  for (i in seq_along(ptm_info)) {
    ptm_info[[i]] <- gsub("\\+|\\[|\\]", "", ptm_info[[i]])
  }
  # Calculate PTM positions for each row
  ptm_positions <- lapply(seq_along(ptm_matches), function(i) {
    ptm_lengths <- attr(ptm_matches[[i]], "match.length")
    if (-1 %in% ptm_lengths) {
      ptm_positions <- NA
    } else {
      ptm_values <- unlist(ptm_matches[i])
      ptm_lengths <- attr(ptm_matches[[i]], "match.length")
      ptm_positions <- ptm_values - cumsum(c(1, ptm_lengths[-length(ptm_lengths)]))
    }
  })
  # Combine results for each row
  rep_values <- ifelse(sapply(ptm_info, length) == 0, 1, sapply(ptm_info, length))
  result <- data.table(
    PTM_position = unlist(ptm_positions),
    reps = rep(rep_values, rep_values)
  )

  # Dynamically add the PTM_mass column
  result[, (PTM_mass_column) := unlist(lapply(ptm_info, function(x)
    if (length(x) > 0)
      x
    else
      NA))]

  for (col in colnames(data)) {
    result[, (col) := unlist(lapply(seq_along(ptm_info), function(i)
      if (length(ptm_info[[i]]) > 0)
        rep(data[[col]][i], length(ptm_info[[i]]))
      else
        data[[col]][i]))]
  }

  if (PTM_annotation & !is.null(PTM_table)) {
    PTM_table[, (PTM_mass_column) := as.character(get(PTM_mass_column))]
    result <- merge(result, PTM_table, by = PTM_mass_column, all.x = TRUE)
  }

  return(result)
}

#' Obtain modification information from Peptide data generated by MSFragger
#'
#' This function takes MSFragger output containing a 'Assigned Modifications' column with
#' PTM information and converts it into a new dataframe with the desired format of peptide
#' sequences and associated PTM information.
#'
#' @param data A data.table with a column containing stripped sequence and a column containing PTM information.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param seq_column The name of the column containing peptide sequences for MSFragger.
#' @param PTM_table A data.table with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param PTM_mass_column The name of the column containing the PTM mass information
#' @return A data.table with 'PTM_position', 'PTM_type', 'reps', and other columns.
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   Peptide = c("DDREDMLVYQAK", "EAAENSLVAYK", "IEAELQDICNDVLELLDK"),
#'   `Assigned Modifications` = c("C-term(15.9949), 6M(-0.98)", "", "N-term(42.0106)"),
#'   Condition1 = c("A", "B", "B"),
#'   Condition2 = c("C", "C", "D")
#' )
#' PTM_table <- data.table(
#'   PTM_mass = c(42.0106, -0.98, 15.9949),
#'   PTM_type = c("Acet", "Amid", "Ox")
#' )
#' mod_column <- "Assigned Modifications"
#' seq_column <- "Peptide"
#' converted_data <- obtain_mod_MSFragger(
#'   data,
#'   mod_column,
#'   seq_column,
#'   PTM_table,
#'   PTM_annotation = TRUE,
#'   PTM_mass_column = "PTM_mass"
#' )
#'
#' @import data.table
#'
#' @export
obtain_mod_MSFragger <- function(data,
                                 mod_column,
                                 seq_column,
                                 PTM_table = NULL,
                                 PTM_annotation = FALSE,
                                 PTM_mass_column) {
  # Ensure data is a data.table
  PTM_table <- as.data.table(PTM_table)

  # Initialize empty lists to store the results
  ptm_positions_list <- list()
  ptm_mass_list <- list()
  reps_list <- list()

  # Iterate through each row
  for (i in seq_len(nrow(data))) {
    mods <- strsplit(data[[mod_column]][i], ", ")[[1]]

    # Handle the case when no modifications are present
    if (length(mods) == 0) {
      ptm_positions_list[[i]] <- NA
      ptm_mass_list[[i]] <- NA
      reps_list[[i]] <- 1
    } else {
      ptm_positions <- integer(length(mods))
      ptm_masses <- numeric(length(mods))
      for (j in seq_along(mods)) {
        if (grepl("C-term", mods[j])) {
          ptm_positions[j] <- nchar(data[[seq_column]][i])
        } else if (grepl("N-term", mods[j])) {
          ptm_positions[j] <- 0
        } else {
          ptm_positions[j] <- as.numeric(gsub("(\\d+)[A-Za-z]*.*", "\\1", mods[j]))
        }
        mass_value <- as.numeric(sub(".*\\((-?\\d+\\.\\d+)\\).*", "\\1", mods[j]))
        ptm_masses[j] <- mass_value
      }

      ptm_positions_list[[i]] <- ptm_positions
      ptm_mass_list[[i]] <- ptm_masses
      reps_list[[i]] <- length(mods)
    }
  }

  # Create new rows for multiple modifications
  new_rows <- lapply(seq_len(nrow(data)), function(i) {
    new_row <- data.table(PTM_position = ptm_positions_list[[i]],
                          reps = rep(reps_list[[i]], reps_list[[i]]),
                          data[i, names(data), with = FALSE])
    # Add the PTM mass column dynamically
    new_row[[PTM_mass_column]] <- ptm_mass_list[[i]]
    return(new_row)
  })

  # Combine the results into a single data.table
  result <- rbindlist(new_rows)

  if (PTM_annotation & !is.null(PTM_table)) {
    result <- merge(result, PTM_table, by = PTM_mass_column, all.x = TRUE)
  }

  return(result)
}

#' Obtain modification information from Peptide data generated by mzIdenML
#'
#' This function takes mzIdenML output containing a 'modification' column with
#' PTM information and converts it into a new dataframe with the desired format of peptide
#' sequences and associated PTM information.
#'
#' @param data A data.table with a column containing stripped sequence and a column containing PTM information.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param seq_column The name of the column containing peptide sequences for mzIdenML.
#' @param PTM_table A data.table with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param PTM_mass_column The name of the column containing the PTM mass information
#' @return A data.table with 'PTM_position', 'PTM_type', 'reps', and other columns.
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   pepseq = c("DDREDMLVYQAK", "EAAENSLVAYK", "IEAELQDICNDVLELLDK"),
#'   modification = c("-0.984016 (10), 15.994915 (13)", NA, "15.994915 (12)"),
#'   Condition1 = c("A", "B", "B"),
#'   Condition2 = c("C", "C", "D")
#' )
#' PTM_table <- data.table(
#'   PTM_mass = c(-0.984016, 15.994915),
#'   PTM_type = c("Amid", "Ox")
#' )
#' mod_column <- "modification"
#' seq_column <- "pepseq"
#' converted_data <- obtain_mod_mzIdenML(
#'   data,
#'   mod_column,
#'   seq_column,
#'   PTM_table,
#'   PTM_annotation = TRUE,
#'   PTM_mass_column = "PTM_mass"
#' )
#'
#' @import data.table
#'
#' @export
obtain_mod_mzIdenML <- function(data,
                                 mod_column,
                                 seq_column,
                                 PTM_table = NULL,
                                 PTM_annotation = FALSE,
                                 PTM_mass_column) {
  # Ensure data is a data.table
  PTM_table <- as.data.table(PTM_table)

  # Initialize empty lists to store the results
  ptm_positions_list <- list()
  ptm_mass_list <- list()
  reps_list <- list()

  # Iterate through each row
  for (i in seq_len(nrow(data))) {
    mods <- strsplit(data[[mod_column]][i], ", ")[[1]]

    # Handle the case when no modifications are present
    if (length(mods) == 0) {
      ptm_positions_list[[i]] <- NA
      ptm_mass_list[[i]] <- NA
      reps_list[[i]] <- 1
    } else {
      ptm_positions <- integer(length(mods))
      ptm_masses <- numeric(length(mods))
      for (j in seq_along(mods)) {
        ptm_positions[j] <- as.numeric(sub(".*\\((\\d+)\\).*", "\\1", mods[j]))
        mass_value <- as.numeric(sub("^(-?\\d+\\.\\d+).*", "\\1", mods[j]))
        ptm_masses[j] <- mass_value
      }

      ptm_positions_list[[i]] <- ptm_positions
      ptm_mass_list[[i]] <- ptm_masses
      reps_list[[i]] <- length(mods)
    }
  }

  # Create new rows for multiple modifications
  new_rows <- lapply(seq_len(nrow(data)), function(i) {
    new_row <- data.table(PTM_position = ptm_positions_list[[i]],
                          reps = rep(reps_list[[i]], reps_list[[i]]),
                          data[i, names(data), with = FALSE])
    # Add the PTM mass column dynamically
    new_row[[PTM_mass_column]] <- ptm_mass_list[[i]]
    return(new_row)
  })

  # Combine the results into a single data.table
  result <- rbindlist(new_rows)

  if (PTM_annotation & !is.null(PTM_table)) {
    result <- merge(result, PTM_table, by = PTM_mass_column, all.x = TRUE)
  }

  return(result)
}

#' Obtain modification information from Peptide data generated by mzTab
#'
#' This function takes mzTab output containing a 'modifications' column with
#' PTM information and converts it into a new dataframe with the desired format of peptide
#' sequences and associated PTM information.
#'
#' @param data A data.table with a column containing stripped sequence and a column containing PTM information.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param seq_column The name of the column containing peptide sequences for mzTab
#' @param PTM_table A data.table with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param PTM_mass_column The name of the column containing the PTM mass information
#' @return A data.table with 'PTM_position', 'PTM_type', 'reps', and other columns.
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   sequence = c("DDREDMLVYQAK", "EAAENSLVAYK", "IEAELQDICNDVLELLDK"),
#'   modifications = c("4-UNIMOD:7,10-UNIMOD:35", NA, "8-UNIMOD:7"),
#'   Condition1 = c("A", "B", "B"),
#'   Condition2 = c("C", "C", "D")
#' )
#' PTM_table <- data.table(
#'   PTM_mass = c("UNIMOD:7", "UNIMOD:35"),
#'   PTM_type = c("Amid", "Ox")
#' )
#' mod_column <- "modifications"
#' seq_column <- "sequence"
#' converted_data <- obtain_mod_mzTab(
#'   data,
#'   mod_column,
#'   seq_column,
#'   PTM_table,
#'   PTM_annotation = TRUE,
#'   PTM_mass_column = "PTM_mass"
#' )
#'
#' @import data.table
#'
#' @export
obtain_mod_mzTab <- function(data,
                                mod_column,
                                seq_column,
                                PTM_table = NULL,
                                PTM_annotation = FALSE,
                                PTM_mass_column) {
  # Ensure data is a data.table
  PTM_table <- as.data.table(PTM_table)

  # Initialize empty lists to store the results
  ptm_positions_list <- list()
  ptm_mass_list <- list()
  reps_list <- list()

  # Iterate through each row
  for (i in seq_len(nrow(data))) {
    mods <- strsplit(data[[mod_column]][i], ",")[[1]]

    # Handle the case when no modifications are present
    if (length(mods) == 0) {
      ptm_positions_list[[i]] <- NA
      ptm_mass_list[[i]] <- NA
      reps_list[[i]] <- 1
    } else {
      ptm_positions <- integer(length(mods))
      ptm_masses <- numeric(length(mods))
      for (j in seq_along(mods)) {
        ptm_positions[j] <- as.numeric(sub("^(\\d+)-.*", "\\1", mods[j]))
        mass_value <- sub("^\\d+-(.*)", "\\1", mods[j])
        ptm_masses[j] <- mass_value
      }

      ptm_positions_list[[i]] <- ptm_positions
      ptm_mass_list[[i]] <- ptm_masses
      reps_list[[i]] <- length(mods)
    }
  }

  # Create new rows for multiple modifications
  new_rows <- lapply(seq_len(nrow(data)), function(i) {
    new_row <- data.table(PTM_position = ptm_positions_list[[i]],
                          reps = rep(reps_list[[i]], reps_list[[i]]),
                          data[i, names(data), with = FALSE])
    # Add the PTM mass column dynamically
    new_row[[PTM_mass_column]] <- ptm_mass_list[[i]]
    return(new_row)
  })

  # Combine the results into a single data.table
  result <- rbindlist(new_rows)

  if (PTM_annotation & !is.null(PTM_table)) {
    result <- merge(result, PTM_table, by = PTM_mass_column, all.x = TRUE)
  }

  return(result)
}

#' Obtain modification information from Peptide data generated by Comet
#'
#' This function takes Comet output containing a column with modified peptide
#' sequences including PTM information and converts it into a new dataframe with the
#' desired format of peptide sequences and associated PTM information.
#'
#' @param data A data.table with a column containing PTM information.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param PTM_table A data.table with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param PTM_mass_column The name of the column containing the PTM mass information
#' @return A data.table with 'PTM_position', 'PTM_type', 'reps', and other columns.
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   modified_peptide = c(
#'     "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
#'     "K.AAQQTGKLVHANFGT.K",
#'     "K.[-0.98]AATVTGKLVHANFGT.K"
#'   ),
#'   plain_peptide = c(
#'     "AAMQRGSLYQCDYSTGSCEPIR",
#'     "AAQQTGKLVHANFGT",
#'     "AATVTGKLVHANFGT"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' PTM_table <- data.table(
#'   PTM_mass = c(57.02, -0.98, 15.9949),
#'   PTM_type = c("Cam", "Amid", "Ox")
#' )
#' mod_column <- 'modified_peptide'
#' PTM_mass_column <- "PTM_mass"
#' converted_data <- obtain_mod_Comet(data, mod_column, PTM_table,
#' PTM_annotation = TRUE, PTM_mass_column)
#'
#' @import data.table
#'
#' @export
obtain_mod_Comet <- function(data,
                             mod_column,
                             PTM_table = NULL,
                             PTM_annotation = FALSE,
                             PTM_mass_column) {
  # Ensure data is a data.table
  PTM_table <- as.data.table(PTM_table)

  # Remove characters before the first dot and after the second dot
  cleaned_sequences <- gsub("^[A-Za-z]\\.|\\.[A-Za-z]$", "", data[[mod_column]])

  # Extract PTM information using regular expression
  ptm_matches <- gregexpr("\\[([^]]+)\\]", cleaned_sequences)
  ptm_info <- regmatches(cleaned_sequences, ptm_matches)
  for (i in seq_along(ptm_info)) {
    ptm_value <- gsub("\\+|\\[|\\]", "", ptm_info[[i]])
    ptm_info[[i]] <- ptm_value
  }

  # Calculate PTM positions for each row
  ptm_positions <- lapply(seq_along(ptm_matches), function(i) {
    ptm_lengths <- attr(ptm_matches[[i]], "match.length")
    if (-1 %in% ptm_lengths) {
      ptm_positions <- NA

    } else {
      ptm_values <- unlist(ptm_matches[i])
      ptm_lengths <- attr(ptm_matches[[i]], "match.length")
      ptm_positions <- ptm_values - cumsum(c(1, ptm_lengths[-length(ptm_lengths)]))
    }
  })

  # Combine results for each row
  rep_values <- ifelse(sapply(ptm_info, length) == 0, 1, sapply(ptm_info, length))
  result <- data.table(
    PTM_position = unlist(ptm_positions),
    reps = rep(rep_values, rep_values)
  )

  # Dynamically add the PTM_mass column
  result[, (PTM_mass_column) := unlist(lapply(ptm_info, function(x)
    if (length(x) > 0)
      x
    else
      NA))]

  for (col in colnames(data)) {
    result[, (col) := unlist(lapply(seq_along(ptm_info), function(i)
      if (length(ptm_info[[i]]) > 0)
        rep(data[[col]][i], length(ptm_info[[i]]))
      else
        data[[col]][i]))]
  }

  if (PTM_annotation & !is.null(PTM_table)) {
    PTM_table[, (PTM_mass_column) := as.character(get(PTM_mass_column))]
    result <- merge(result, PTM_table, by = PTM_mass_column, all.x = TRUE)
  }

  return(result)
}

#' Obtain modification information from Peptide data generated by DIA-NN
#'
#' This function takes DIA-NN output containing a column with modified peptide
#' sequences including PTM information and converts it into a new dataframe with the
#' desired format of peptide sequences and associated PTM information.
#'
#' @param data A dataframe with 'Stripped.Sequence' column and 'Modified.Sequence' column containing modified peptide sequences.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param PTM_table A dataframe with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param PTM_mass_column The name of the column containing the PTM mass information
#' @return A dataframe with 'Peptide', 'PTM_position', and 'PTM_type' columns.
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   Modified.Sequence = c(
#'     "AAAAGPGAALS(UniMod:21)PRPC(UniMod:4)DSDPATPGAQSPK",
#'     "AAAASAAEAGIATTGTEDSDDALLK",
#'     "AAAAALSGSPPQTEKPT(UniMod:21)HYR"
#'   ),
#'   Stripped.Sequence = c(
#'     "AAAAGPGAALSPRPCDSDPATPGAQSPK",
#'     "AAAASAAEAGIATTGTEDSDDALLK",
#'     "AAAAALSGSPPQTEKPTHYR"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' PTM_table <- data.table(PTM_mass = c('UniMod:21', 'UniMod:4'),
#'                         PTM_type = c("Phos", "Cam"))
#' converted_data <- obtain_mod_DIANN(
#'   data,
#'   'Modified.Sequence',
#'   PTM_table,
#'   PTM_annotation = TRUE,
#'   PTM_mass_column = "PTM_mass"
#' )
#'
#' @import data.table
#'
#' @export
obtain_mod_DIANN <- function(data,
                             mod_column,
                             PTM_table = NULL,
                             PTM_annotation = FALSE,
                             PTM_mass_column) {
  # Ensure data is a data.table
  PTM_table <- as.data.table(PTM_table)

  # Extract PTM information using regular expression
  ptm_matches <- gregexpr("\\(([^)]+)\\)", data[[mod_column]])
  ptm_info <- regmatches(data[[mod_column]], ptm_matches)
  for (i in seq_along(ptm_info)) {
    ptm_info[[i]] <- gsub("\\+|\\(|\\)", "", ptm_info[[i]])
  }

  # Calculate PTM positions for each row
  ptm_positions <- lapply(seq_along(ptm_matches), function(i) {
    ptm_lengths <- attr(ptm_matches[[i]], "match.length")
    if (-1 %in% ptm_lengths) {
      ptm_positions <- NA

    } else {
      ptm_values <- unlist(ptm_matches[i])
      ptm_lengths <- attr(ptm_matches[[i]], "match.length")
      ptm_positions <- ptm_values - cumsum(c(1, ptm_lengths[-length(ptm_lengths)]))
    }
  })

  # Combine results for each row
  rep_values <- ifelse(sapply(ptm_info, length) == 0, 1, sapply(ptm_info, length))
  result <- data.table(
    PTM_position = unlist(ptm_positions),
    reps = rep(rep_values, rep_values)
  )

  # Dynamically add the PTM_mass column
  result[, (PTM_mass_column) := unlist(lapply(ptm_info, function(x)
    if (length(x) > 0)
      x
    else
      NA))]

  for (col in colnames(data)) {
    result[, (col) := unlist(lapply(seq_along(ptm_info), function(i)
      if (length(ptm_info[[i]]) > 0)
        rep(data[[col]][i], length(ptm_info[[i]]))
      else
        data[[col]][i]))]
  }

  if (PTM_annotation & !is.null(PTM_table)) {
    PTM_table[, (PTM_mass_column) := as.character(get(PTM_mass_column))]
    result <- merge(result, PTM_table, by = PTM_mass_column, all.x = TRUE)
  }

  return(result)
}

#' Obtain modification information from Peptide data generated by Skyline
#'
#' This function takes Skyline output containing a column with modified peptide
#' sequences including PTM information and converts it into a new dataframe with the
#' desired format of peptide sequences and associated PTM information.
#'
#' @param data A data.table with a column containing PTM information.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param PTM_table A data.table with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param PTM_mass_column The name of the column containing the PTM mass information
#' @return A data.table with 'PTM_position', 'PTM_type', 'reps', and other columns.
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   'Peptide Modified Sequence' = c(
#'     "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
#'     "AAQQTGKLVHANFGT",
#'     "[-0.98]AATVTGKLVHANFGT"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' PTM_table <- data.table(
#'   PTM_mass = c(57.02, -0.98, 15.9949),
#'   PTM_type = c("Cam", "Amid", "Ox")
#' )
#' converted_data <- obtain_mod_Skyline(
#'   data,
#'   'Peptide Modified Sequence',
#'   PTM_table,
#'   PTM_annotation = TRUE,
#'   PTM_mass_column = "PTM_mass"
#' )
#'
#' @import data.table
#'
#' @export
obtain_mod_Skyline <- function(data,
                               mod_column,
                               PTM_table,
                               PTM_annotation = FALSE,
                               PTM_mass_column) {
  # Ensure data is a data.table
  PTM_table <- as.data.table(PTM_table)

  # Extract PTM information using regular expression
  ptm_matches <- gregexpr("\\[([^]]+)\\]", data[[mod_column]])
  ptm_info <- regmatches(data[[mod_column]], ptm_matches)
  for (i in seq_along(ptm_info)) {
    ptm_value <- gsub("\\+|\\[|\\]", "", ptm_info[[i]])
    ptm_info[[i]] <- ptm_value
  }

  # Calculate PTM positions for each row
  ptm_positions <- lapply(seq_along(ptm_matches), function(i) {
    ptm_lengths <- attr(ptm_matches[[i]], "match.length")
    if (-1 %in% ptm_lengths) {
      ptm_positions <- NA

    } else {
      ptm_values <- unlist(ptm_matches[i])
      ptm_lengths <- attr(ptm_matches[[i]], "match.length")
      ptm_positions <- ptm_values - cumsum(c(1, ptm_lengths[-length(ptm_lengths)]))
    }
  })

  # Combine results for each row
  rep_values <- ifelse(sapply(ptm_info, length) == 0, 1, sapply(ptm_info, length))
  result <- data.table(
    PTM_position = unlist(ptm_positions),
    reps = rep(rep_values, rep_values)
  )

  # Dynamically add the PTM_mass column
  result[, (PTM_mass_column) := unlist(lapply(ptm_info, function(x)
    if (length(x) > 0)
      x
    else
      NA))]

  for (col in colnames(data)) {
    result[, (col) := unlist(lapply(seq_along(ptm_info), function(i)
      if (length(ptm_info[[i]]) > 0)
        rep(data[[col]][i], length(ptm_info[[i]]))
      else
        data[[col]][i]))]
  }

  if (PTM_annotation & !is.null(PTM_table)) {
    PTM_table[, (PTM_mass_column) := as.character(get(PTM_mass_column))]
    result <- merge(result, PTM_table, by = PTM_mass_column, all.x = TRUE)
  }

  return(result)
}

#' Obtain modification information from Peptide data generated by Maxquant
#'
#' This function takes Maxquant output containing a column with modified peptide sequences
#' including PTM information and converts it into a new dataframe with the desired format of peptide
#' sequences and associated PTM information.
#'
#' @param data A data.table with a column containing modified peptide sequences.
#' @param mod_column The name of the column containing the modified peptide sequences.
#' @param PTM_table A data.table with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.
#' @param PTM_annotation A logical value indicating whether to include PTM annotation information in the result.
#' @param PTM_mass_column The name of the column containing the PTM mass information
#' @return A data.table with 'PTM_position', 'PTM_type', 'reps', and other columns.
#'
#' @examples
#' library(data.table)
#' data <- data.table(
#'   'Modified sequence' = c(
#'     "_GLGPSPAGDGPS(Phospho (STY))GSGK_",
#'     "_HSSYPAGTEDDEGM(Oxidation (M))GEEPSPFR_",
#'     "_HSSYPAGTEDDEGM(Oxidation (M))GEEPS(Phospho (STY))PFR_"
#'   ),
#'   Condition = c("A", "B", "B")
#' )
#' PTM_table <- data.table(
#'   PTM_mass = c('Phospho (STY)', 'Oxidation (M)'),
#'   PTM_type = c("Phos", "Ox")
#' )
#' converted_data <- obtain_mod_Maxquant(
#'   data,
#'   'Modified sequence',
#'   PTM_table,
#'   PTM_annotation = TRUE,
#'   PTM_mass_column = "PTM_mass"
#' )
#'
#' @import data.table
#'
#' @export
obtain_mod_Maxquant <- function(data,
                                mod_column,
                                PTM_table = NULL,
                                PTM_annotation = FALSE,
                                PTM_mass_column) {
  # Ensure data is a data.table
  PTM_table <- as.data.table(PTM_table)

  # Extract PTM information using regular expression
  ptm_matches <- gregexpr("\\(([^]]+)\\)", data[[mod_column]])
  ptm_info <- regmatches(data[[mod_column]], ptm_matches)
  for (i in seq_along(ptm_info)) {
    ptm_info[[i]] <- gsub("\\+|\\(|\\)", "", ptm_info[[i]])
  }

  # Calculate PTM positions for each row
  ptm_positions <- lapply(seq_along(ptm_matches), function(i) {
    ptm_lengths <- attr(ptm_matches[[i]], "match.length")
    if (-1 %in% ptm_lengths) {
      ptm_positions <- NA

    } else {
      ptm_values <- unlist(ptm_matches[i])
      ptm_lengths <- attr(ptm_matches[[i]], "match.length")
      ptm_positions <- ptm_values - cumsum(c(1, ptm_lengths[-length(ptm_lengths)])) - 1
    }
  })

  # Combine results for each row
  rep_values <- ifelse(sapply(ptm_info, length) == 0, 1, sapply(ptm_info, length))
  result <- data.table(
    PTM_position = unlist(ptm_positions),
    reps = rep(rep_values, rep_values)
  )

  # Dynamically add the PTM_mass column
  result[, (PTM_mass_column) := unlist(lapply(ptm_info, function(x)
    if (length(x) > 0)
      x
    else
      NA))]

  # Add other columns to the result by recycling values
  for (col in colnames(data)) {
    result[, (col) := unlist(lapply(seq_along(ptm_info), function(i)
      if (length(ptm_info[[i]]) > 0)
        rep(data[[col]][i], length(ptm_info[[i]]))
      else
        data[[col]][i]))]
  }

  if (PTM_annotation & !is.null(PTM_table)) {
    PTM_table[, (PTM_mass_column) := as.character(get(PTM_mass_column))]
    result <- merge(result, PTM_table, by = PTM_mass_column, all.x = TRUE)
  }

  return(result)
}

