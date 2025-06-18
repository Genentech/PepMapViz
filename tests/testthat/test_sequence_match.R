test_that("match_and_calculate_positions works correctly", {

  # Define the input data
  peptide_data <- data.frame(
    Sequence = c("AAA", "BBB", "CCC", "DDD"),
    Condition_1 = c("Drug1", "Drug1", "Drug2", "Drug2"),
    Condition_2 = c("Donor1", "Donor2", "Donor1", "Donor2"),
    Region_1 = c("VH", "VL", "VH", "VL"),
    Region_2 = c("Arm_1", "Arm_2", "Arm_1", "Arm_2"),
    Area = c(100, 2, 4, NA)
  )

  whole_seq <- data.frame(
    Region_Sequence = c(
      "XYZAAA",
      "XYZCCC",
      "XYZBBB",
      "XYZDDD",
      "XYZAAB",
      "XYZCCD",
      "XYZBBB",
      "XYZDDD",
      "XYZAAA",
      "XYZCCC",
      "XYZBBB",
      "XYZDDD",
      "XYZAAB",
      "XYZCCD",
      "XYZBBB",
      "XYZDDD"
    ),
    Condition_1 = c(
      "Drug1",
      "Drug1",
      "Drug2",
      "Drug2",
      "Drug1",
      "Drug1",
      "Drug2",
      "Drug2",
      "Drug1",
      "Drug1",
      "Drug2",
      "Drug2",
      "Drug1",
      "Drug1",
      "Drug2",
      "Drug2"
    ),
    Condition_2 = c(
      "Donor1",
      "Donor1",
      "Donor1",
      "Donor1",
      "Donor1",
      "Donor1",
      "Donor1",
      "Donor1",
      "Donor2",
      "Donor2",
      "Donor2",
      "Donor2",
      "Donor2",
      "Donor2",
      "Donor2",
      "Donor2"
    ),
    Region_1 = c(
      "VH",
      "VL",
      "VH",
      "VL",
      "VH",
      "VL",
      "VH",
      "VL",
      "VH",
      "VL",
      "VH",
      "VL",
      "VH",
      "VL",
      "VH",
      "VL"
    ),
    Region_2 = c(
      "Arm_1",
      "Arm_1",
      "Arm_1",
      "Arm_1",
      "Arm_2",
      "Arm_2",
      "Arm_2",
      "Arm_2",
      "Arm_1",
      "Arm_1",
      "Arm_1",
      "Arm_1",
      "Arm_2",
      "Arm_2",
      "Arm_2",
      "Arm_2"
    )
  )

  match_columns <- c("Condition_1", "Condition_2", "Region_1")
  column_keep <- c("Region_2")
  sequence_length <- c(1, 5)
  column <- "Sequence"

  # Call the function
  matching_result <- match_and_calculate_positions(peptide_data,
                                                   column,
                                                   whole_seq,
                                                   match_columns,
                                                   sequence_length,
                                                   column_keep)

  # Check the results
  expect_equal(nrow(matching_result), 3)
  expect_equal(matching_result$Sequence, c("AAA", "DDD", "DDD"))
  expect_equal(matching_result$start_position, c("4", "4", "4"))
  expect_equal(matching_result$end_position, c("6", "6", "6"))
  expect_equal(matching_result$Region_2, c("Arm_1", "Arm_2", "Arm_2"))
})
