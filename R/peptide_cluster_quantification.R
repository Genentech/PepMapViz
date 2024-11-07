#' Peptide Quantification
#'
#' @param whole_seq A dataframe holding whole sequence information. 'Region_Sequence' column is required for the sequence information. Change the column name if it is different than 'Region_Sequence'.
#' @param matching_result The dataframe that contains the matched results and PTM information.
#' @param matching_columns Vector of column names that should match between each row of 'whole_seq' and the 'matching_result' dataframe.
#' @param distinct_columns Vector of column names that should be used to calculate PSM or Area separately for each unique combination of these columns.
#' @param quantify_method A string indicating the quantification method. It can be either "PSM" or "Area".
#' @param area_column The name of the column in 'matching_result' that contains the area/intensity information. Required if quantify_method is "Area".
#' @param with_PTM A boolean parameter indicating whether PTM should be considered during calculation. Default is `FALSE`.
#' @param reps A boolean parameter indicating whether the area/intensity should be divided by the number of replicates. Default is `FALSE`.
#'
#' @return Returns a dataframe containing the calculated PSM or Area for each record in 'whole_seq'.
#'
#' @examples
#' whole_seq <- data.frame(
#'   Region_Sequence = c(
#'     "XYZAAA",
#'     "XYZCCC",
#'     "XYZBBB",
#'     "XYZDDD",
#'     "XYZAAB",
#'     "XYZCCD",
#'     "XYZBBB",
#'     "XYZDDD",
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
#' matching_result <- data.frame(
#'   Sequence = c("AAA", "DDD", "DDD"),
#'   Condition_1 = c("Drug1", "Drug2", "Drug2"),
#'   Condition_2 = c("Donor1", "Donor2", "Donor2"),
#'   Region_1 = c("VH", "VL", "VL"),
#'   Region_2 = c("Arm_1", "Arm_2", "Arm_2"),
#'   Start_Position = c(4, 4, 4),
#'   End_Position = c(6, 6, 6),
#'   PTM_position = c(NA, 2, 0),
#'   PTM_type = c(NA, "O", "C"),
#'   Area = c(100, 200, 200),
#'   reps = c(1, 2, 2)
#' )
#' matching_columns <- c("Condition_1", "Region_2")
#' area_column <- "Area"
#' data_with_quantification <- peptide_quantification(
#'   whole_seq,
#'   matching_result,
#'   matching_columns,
#'   distinct_columns = c("Condition_2", "Region_1"),
#'   quantify_method = "Area",
#'   area_column = area_column,
#'   with_PTM = TRUE,
#'   reps = TRUE
#' )
#' @export
#'
peptide_quantification <- function(whole_seq, matching_result, matching_columns, distinct_columns, quantify_method, area_column = NULL, with_PTM = FALSE, reps = FALSE) {
  if (quantify_method == "PSM") {
    data_with_quantification <- calculate_all_PSM(whole_seq, matching_result, matching_columns, distinct_columns, with_PTM, reps)
  } else if (quantify_method == "Area") {
    if (is.null(area_column)) {
      stop("The area_column parameter is required when quantify_method is 'Area'")
    }
    data_with_quantification <- calculate_all_Area(whole_seq, matching_result, matching_columns, distinct_columns, area_column, with_PTM, reps)
  } else {
    stop("Invalid quantify_method. It should be either 'PSM' or 'Area'")
  }

  return(data_with_quantification)
}

#' Calculate Spectra Count (PSM) for one row of the input sequence dataframe
#'
#' @param row A row of dataframe containing the sequence for the 'Character' column in region_data.
#' @param matching_result The dataframe that contains the matched results and PTM information.
#' @param matching_columns Vector of column names that should match between the 'row' and 'matching_result' dataframes.
#' @param distinct_columns Vector of column names that should be used to calculate PSM separately for each unique combination of these columns.
#' @param with_PTM A boolean parameter indicating whether PTM should be considered. If `with_PTM = TRUE`,
#' the function will also add 'PTM' and 'PTM_type' to the result 'region_data' dataframe. Default is `FALSE`.
#' @param reps A boolean parameter indicating whether the area/intensity should be divided by the number of replicates. Default is `FALSE`.
#'
#' @return This function returns the modified `region_data` dataframe that includes the "PSM" column, and optionally "PTM" and "PTM_type" columns.
#' If the 'filter_conditions' do not match, an empty dataframe will be returned early.
#' An AttributeError is raised if 'PTM_position' and 'PTM_type' columns do not exist in the 'result' dataframe when 'with_PTM' is `TRUE`.
#'
#' @examples
#' row <- data.frame(
#'  Region_Sequence = c("XYZDDD"),
#'  Condition_1 = c("Drug2"),
#'  Region_1 = c("VL"),
#'  Region_2 = c("Arm_2")
#' )
#' matching_result <- data.frame(
#'   Sequence = c("AAA", "DDD", "DDD"),
#'   Condition_1 = c("Drug1", "Drug2", "Drug2"),
#'   Condition_2 = c("Donor1", "Donor2", "Donor2"),
#'   Region_1 = c("VH", "VL", "VL"),
#'   Region_2 = c("Arm_1", "Arm_2", "Arm_2"),
#'   Start_Position = c(4, 4, 4),
#'   End_Position = c(6, 6, 6),
#'   PTM_position = c(NA, 2, 0),
#'   PTM_type = c(NA,"O","C"),
#'   Area = c(100, 200, 200),
#'   reps = c(1, 2, 2)
#' )
#' matching_columns <- c("Condition_1", "Region_2")
#' result <- calculate_PSM(
#'   row,
#'   matching_result,
#'   matching_columns,
#'   distinct_columns = c("Condition_2", "Region_1"),
#'   with_PTM = TRUE,
#'   reps = TRUE
#' )
#'
#' @export
#'
calculate_PSM <- function(row, matching_result, matching_columns, distinct_columns, with_PTM = FALSE, reps = FALSE) {
  # Create a new dataframe for each character in Region_Sequence with its position information
  region_data <- data.frame(
    Character = strsplit(row$Region_Sequence, "")[[1]],
    Position = seq_along(strsplit(row$Region_Sequence, "")[[1]]),
    stringsAsFactors = FALSE
  )

  # Add matching_columns to region_data
  for (column in matching_columns) {
    region_data[[column]] = row[[column]]
  }

  # Initialize a logical vector to store filtering conditions
  filter_condition <- rep(TRUE, nrow(matching_result))

  # Ensure all conditions are met simultaneously
  for (column in matching_columns) {
    filter_condition <- filter_condition & matching_result[[column]] == row[[column]]
  }

  # If the filtered condition is all FALSE, return an empty dataframe
  if (all(!filter_condition)) {
    region_data$PSM <- 0
    if (with_PTM){
      region_data$PTM <- FALSE
      region_data$PTM_type <- NA
    }
    if (is.null(distinct_columns)) {
      return(region_data)
    } else {
      for (col in distinct_columns) {
        region_data[[col]] <- NA
      }
      return(region_data)
    }
  }

  # Subset result based on filtering condition
  result_filtered <- matching_result[filter_condition, ]
  # Apply as.numeric to the specified columns
  columns_to_convert <- c("Start_Position", "End_Position", "reps", "PTM_position")
  result_filtered[columns_to_convert] <- lapply(result_filtered[columns_to_convert], as.numeric)

  if (is.null(distinct_columns)) {
    region_data$PSM <- sapply(region_data$Position, function(pos) {
      if (reps) {
        sum((result_filtered$Start_Position <= pos & result_filtered$End_Position >= pos) / result_filtered$reps)
      } else {
        sum(result_filtered$Start_Position <= pos & result_filtered$End_Position >= pos)
      }
    })
  } else {
    # Get unique combinations of distinct columns
    distinct_combinations <- unique(result_filtered[, distinct_columns, drop = FALSE])

    combined_result <- data.frame()
    # Calculate PSM for each character in region for each distinct combination
    for (i in 1:nrow(distinct_combinations)) {
      distinct_result <- region_data
      combo <- distinct_combinations[i, , drop = FALSE]
      combo_condition <- rep(TRUE, nrow(result_filtered))
      for (col in distinct_columns) {
        combo_condition <- combo_condition & result_filtered[[col]] == combo[[col]]
      }
      result_combo <- result_filtered[combo_condition, ]

      distinct_result$PSM <- sapply(distinct_result$Position, function(pos) {
        if (reps) {
          sum((result_combo$Start_Position <= pos & result_combo$End_Position >= pos) / result_combo$reps)
        } else {
          sum(result_combo$Start_Position <= pos & result_combo$End_Position >= pos)
        }
      })
    for (col in distinct_columns) {
      distinct_result[[col]] <- combo[[col]]
    }
    combined_result <- rbind(combined_result, distinct_result)
    }
  }

  # Modify combined_result to include PTM information if with_PTM = TRUE
  if (with_PTM) {
    combined_result$PTM <- FALSE
    combined_result$PTM_type <- NA
    if (!all(c('PTM_position', 'PTM_type') %in% names(matching_result))) {
      stop("The result dataframe is required to have 'PTM_position' and 'PTM_type' columns when 'with_PTM' is TRUE")
    } else {
      for (i in 1:nrow(result_filtered)) {
        if (!is.na(result_filtered$PTM_position[i])) {
          PTM_pos <- result_filtered$Start_Position[i] + result_filtered$PTM_position[i] - 1

          if (!is.null(distinct_columns)) {
            condition <- combined_result$Position == PTM_pos
            for (col in distinct_columns) {
              condition <- condition & combined_result[[col]] == result_filtered[[col]][i]
            }
            combined_result$PTM[condition] <- TRUE
            combined_result$PTM_type[condition] <- result_filtered$PTM_type[i]
          } else {
            combined_result$PTM[combined_result$Position == PTM_pos] <- TRUE
            combined_result$PTM_type[combined_result$Position == PTM_pos] <- result_filtered$PTM_type[i]
          }
        }
      }
    }
  }

  return(combined_result)
}

#' Calculate Spectra Count (PSM) for the whole input sequence dataframe
#'
#' @param whole_seq A dataframe holding whole sequence information. 'Region_Sequence' column is required for the sequence information. Change the column name if it is different than 'Region_Sequence'.
#' @param matching_result The dataframe that contains the matched results and PTM information.
#' @param matching_columns Vector of column names that should match between each row of 'whole_seq' and the 'matching_result' dataframe.
#' @param distinct_columns Vector of column names that should be used to calculate PSM separately for each unique combination of these columns.
#' @param with_PTM A boolean parameter indicating whether PTM should be considered during calculation of PSM. Default is `FALSE`.
#' @param reps A boolean parameter indicating whether the area/intensity should be divided by the number of replicates. Default is `FALSE`.
#'
#' @return Returns `data_with_psm`, a dataframe contains calculated PSM for each record in 'whole_seq'.
#'
#' @examples
#' whole_seq <- data.frame(
#'   Region_Sequence = c(
#'     "XYZAAA",
#'     "XYZCCC",
#'     "XYZBBB",
#'     "XYZDDD",
#'     "XYZAAB",
#'     "XYZCCD",
#'     "XYZBBB",
#'     "XYZDDD",
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
#' matching_result <- data.frame(
#'   Sequence = c("AAA", "DDD", "DDD"),
#'   Condition_1 = c("Drug1", "Drug2", "Drug2"),
#'   Condition_2 = c("Donor1", "Donor2", "Donor2"),
#'   Region_1 = c("VH", "VL", "VL"),
#'   Region_2 = c("Arm_1", "Arm_2", "Arm_2"),
#'   Start_Position = c(4, 4, 4),
#'   End_Position = c(6, 6, 6),
#'   PTM_position = c(NA, 2, 0),
#'   PTM_type = c(NA, "O", "C"),
#'   Area = c(100, 200, 200),
#'   reps = c(1, 2, 2)
#' )
#' matching_columns <- c("Condition_1", "Region_2")
#' data_with_psm <- calculate_all_PSM(
#'   whole_seq,
#'   matching_result,
#'   matching_columns,
#'   distinct_columns = c("Condition_2", "Region_1"),
#'   with_PTM = TRUE,
#'   reps = TRUE
#' )
#' @export
#'
calculate_all_PSM <- function(whole_seq, matching_result, matching_columns, distinct_columns, with_PTM = FALSE, reps = FALSE) {
  # Apply the calculate_PSM function on each row of whole_seq
  if (!'Region_Sequence' %in% names(whole_seq)) {
    stop("The whole_seq dataframe is required to have 'Region_Sequence' column")
  } else {
  data_with_psm <- do.call(rbind, apply(whole_seq, 1, function(row) {
    calculate_PSM(as.data.frame(t(row), stringsAsFactors=FALSE), matching_result, matching_columns, distinct_columns, with_PTM, reps)
  }))
  }

  return(data_with_psm)
}


#' Calculate Area/Intensity for one row of the input sequence dataframe
#'
#' @param row A row of dataframe containing the sequence for the 'Character' column in region_data.
#' @param matching_result The dataframe that contains the matched results and PTM information.
#' @param matching_columns Vector of column names that should match between the 'row' and 'matching_result' dataframes.
#' @param distinct_columns Vector of column names that should be used to calculate Area separately for each unique combination of these columns.
#' @param area_column The name of the column in 'matching_result' that contains the area/intensity information.
#' @param with_PTM A boolean parameter indicating whether PTM should be considered. If `with_PTM = TRUE`,
#' the function will also add 'PTM' and 'PTM_type' to the result 'region_data' dataframe. Default is `FALSE`.
#' @param reps A boolean parameter indicating whether the area/intensity should be divided by the number of replicates. Default is `FALSE`.
#'
#' @return This function returns the modified `region_data` dataframe that includes the "Area" column, and optionally "PTM" and "PTM_type" columns.
#' If the 'filter_conditions' do not match, an empty dataframe will be returned early.
#' An AttributeError is raised if 'PTM_position' and 'PTM_type' columns do not exist in the 'result' dataframe when 'with_PTM' is `TRUE`.
#'
#' @examples
#' row <- data.frame(
#'  Region_Sequence = c("XYZAAA"),
#'  Condition_1 = c("Drug1"),
#'  Condition_2 = c("Donor1"),
#'  Region_1 = c("VH"),
#'  Region_2 = c("Arm_1")
#' )
#' matching_result <- data.frame(
#'   Sequence = c("AAA", "DDD", "DDD"),
#'   Condition_1 = c("Drug1", "Drug2", "Drug2"),
#'   Condition_2 = c("Donor1", "Donor2", "Donor2"),
#'   Region_1 = c("VH", "VL", "VL"),
#'   Region_2 = c("Arm_1", "Arm_2", "Arm_2"),
#'   Start_Position = c(4, 4, 4),
#'   End_Position = c(6, 6, 6),
#'   PTM_position = c(NA, 2, 0),
#'   PTM_type = c(NA,"O","C"),
#'   Area = c(100, 200, 200),
#'   reps = c(1, 2, 2)
#' )
#' matching_columns <- c("Condition_1", "Region_2")
#' area_column <- "Area"
#' data_with_area <- calculate_Area(
#'   row,
#'   matching_result,
#'   matching_columns,
#'   distinct_columns = c("Condition_2", "Region_1"),
#'   area_column,
#'   with_PTM = TRUE,
#'   reps = TRUE
#' )
#'
#' @export
#'
calculate_Area <- function(row, matching_result, matching_columns, distinct_columns = NULL, area_column, with_PTM = FALSE, reps=FALSE) {
  # Create a new dataframe for each character in Region_Sequence with its position information
  region_data <- data.frame(
    Character = strsplit(row$Region_Sequence, "")[[1]],
    Position = seq_along(strsplit(row$Region_Sequence, "")[[1]]),
    stringsAsFactors = FALSE
  )

  # Add matching_columns to region_data
  for (column in matching_columns) {
    region_data[[column]] = row[[column]]
  }

  # Initialize a logical vector to store filtering conditions
  filter_condition <- rep(TRUE, nrow(matching_result))

  # Ensure all conditions are met simultaneously
  for (column in matching_columns) {
    filter_condition <- filter_condition & matching_result[[column]] == row[[column]]
  }

  # If the filtered condition is all FALSE, return an empty dataframe
  if (all(!filter_condition)) {
    region_data$Area <- 0
    if (with_PTM){
      region_data$PTM <- FALSE
      region_data$PTM_type <- NA
    }
    if (is.null(distinct_columns)) {
      return(region_data)
    } else {
      for (col in distinct_columns) {
        region_data[[col]] <- NA
      }
      return(region_data)
    }
  }

  # Subset result based on filtering condition
  result_filtered <- matching_result[filter_condition, ]
  # Apply as.numeric to the specified columns
  columns_to_convert <- c("Start_Position", "End_Position", "reps", "PTM_position")
  columns_to_convert <- c(columns_to_convert, area_column)
  result_filtered[columns_to_convert] <- lapply(result_filtered[columns_to_convert], as.numeric)

  # Remove rows where the area column is NA
  result_filtered <- result_filtered[!is.na(result_filtered[[area_column]]), ]

  # For rows with the same value in the Sequence column and distinct columns,
  # keep only the row with the largest value in the area column
  result_filtered$unique_id <- apply(result_filtered[, c("Sequence", distinct_columns)], 1, paste, collapse = "_")
  split_data <- split(result_filtered, result_filtered$unique_id)
  result_filtered <- do.call(rbind, lapply(split_data, function(group) {
    group[which.max(group[[area_column]]),]
  }))
  result_filtered$unique_id <- NULL
  rownames(result_filtered) <- NULL

  # Calculate Area for each character in region
  if (is.null(distinct_columns)) {
    region_data$Area <- sapply(region_data$Position, function(pos) {
      if (reps) {
        sum(ifelse(result_filtered$Start_Position <= pos & result_filtered$End_Position >= pos,
                   result_filtered[[area_column]] / result_filtered$reps, 0))
      } else {
        sum(result_filtered$Start_Position <= pos & result_filtered$End_Position >= pos)
      }
    })
  } else {
    # Get unique combinations of distinct columns
    distinct_combinations <- unique(result_filtered[, distinct_columns, drop = FALSE])

    combined_result <- data.frame()
    # Calculate PSM for each character in region for each distinct combination
    for (i in 1:nrow(distinct_combinations)) {
      distinct_result <- region_data
      combo <- distinct_combinations[i, , drop = FALSE]
      combo_condition <- rep(TRUE, nrow(result_filtered))
      for (col in distinct_columns) {
        combo_condition <- combo_condition & result_filtered[[col]] == combo[[col]]
      }
      result_combo <- result_filtered[combo_condition, ]

      distinct_result$Area <- sapply(distinct_result$Position, function(pos) {
        if (reps) {
          sum(ifelse(result_combo$Start_Position <= pos & result_combo$End_Position >= pos,
                     result_combo[[area_column]] / result_combo$reps, 0))
        } else {
          sum(result_combo$Start_Position <= pos & result_combo$End_Position >= pos)
        }
      })
      for (col in distinct_columns) {
        distinct_result[[col]] <- combo[[col]]
      }
      combined_result <- rbind(combined_result, distinct_result)
    }
  }

  # Apply log2 transformation to the area column, keeping 0 values as 0
  combined_result[[area_column]] <- ifelse(combined_result[[area_column]] == 0,
                                           0,
                                           log2(combined_result[[area_column]]))

  # Modify combined_result to include PTM information if with_PTM = TRUE
  if (with_PTM) {
    combined_result$PTM <- FALSE
    combined_result$PTM_type <- NA
    if (!all(c('PTM_position', 'PTM_type') %in% names(matching_result))) {
      stop("The result dataframe is required to have 'PTM_position' and 'PTM_type' columns when 'with_PTM' is TRUE")
    } else {
      for (i in 1:nrow(result_filtered)) {
        if (!is.na(result_filtered$PTM_position[i])) {
          PTM_pos <- result_filtered$Start_Position[i] + result_filtered$PTM_position[i] - 1

          if (!is.null(distinct_columns)) {
            condition <- combined_result$Position == PTM_pos
            for (col in distinct_columns) {
              condition <- condition & combined_result[[col]] == result_filtered[[col]][i]
            }
            combined_result$PTM[condition] <- TRUE
            combined_result$PTM_type[condition] <- result_filtered$PTM_type[i]
          } else {
            combined_result$PTM[combined_result$Position == PTM_pos] <- TRUE
            combined_result$PTM_type[combined_result$Position == PTM_pos] <- result_filtered$PTM_type[i]
          }
        }
      }
    }
  }



  return(combined_result)
}

#' Calculate Area/Intensity for the whole input sequence dataframe
#'
#' @param whole_seq A dataframe holding whole sequence information. 'Region_Sequence' column is required for the sequence information. Change the column name if it is different than 'Region_Sequence'.
#' @param matching_result The dataframe that contains the matched results and PTM information.
#' @param matching_columns Vector of column names that should match between each row of 'whole_seq' and the 'matching_result' dataframe.
#' @param distinct_columns Vector of column names that should be used to calculate Area separately for each unique combination of these columns.
#' @param area_column The name of the column in 'matching_result' that contains the area/intensity information.
#' @param with_PTM A boolean parameter indicating whether PTM should be considered during calculation of Area. Default is `FALSE`.
#' @param reps A boolean parameter indicating whether the area/intensity should be divided by the number of replicates. Default is `FALSE`.
#'
#' @return Returns `data_with_area`, a dataframe contains calculated Area for each record in 'whole_seq'.
#'
#' @examples
#' whole_seq <- data.frame(
#'   Region_Sequence = c(
#'     "XYZAAA",
#'     "XYZCCC",
#'     "XYZBBB",
#'     "XYZDDD",
#'     "XYZAAB",
#'     "XYZCCD",
#'     "XYZBBB",
#'     "XYZDDD",
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
#' matching_result <- data.frame(
#'   Sequence = c("AAA", "DDD", "DDD"),
#'   Condition_1 = c("Drug1", "Drug2", "Drug2"),
#'   Condition_2 = c("Donor1", "Donor2", "Donor2"),
#'   Region_1 = c("VH", "VL", "VL"),
#'   Region_2 = c("Arm_1", "Arm_2", "Arm_2"),
#'   Start_Position = c(4, 4, 4),
#'   End_Position = c(6, 6, 6),
#'   PTM_position = c(NA, 2, 0),
#'   PTM_type = c(NA, "O", "C"),
#'   Area = c(100, 200, 200),
#'   reps = c(1, 2, 2)
#' )
#' matching_columns <- c("Condition_1", "Region_2")
#' area_column <- "Area"
#' data_with_area <- calculate_all_Area(
#'   whole_seq,
#'   matching_result,
#'   matching_columns,
#'   distinct_columns = c("Condition_2", "Region_1"),
#'   area_column,
#'   with_PTM = TRUE,
#'   reps = TRUE
#' )
#'
#' @export
#'
calculate_all_Area <- function(whole_seq, matching_result, matching_columns, distinct_columns, area_column, with_PTM = FALSE, reps = FALSE) {
  # Apply the calculate_PSM function on each row of whole_seq
  if (!'Region_Sequence' %in% names(whole_seq)) {
    stop("The whole_seq dataframe is required to have 'Region_Sequence' column")
  } else {
    data_with_area <- do.call(rbind, apply(whole_seq, 1, function(row) {
      calculate_Area(as.data.frame(t(row), stringsAsFactors=FALSE), matching_result, matching_columns, distinct_columns, area_column, with_PTM, reps)
    }))
  }

  return(data_with_area)
}
