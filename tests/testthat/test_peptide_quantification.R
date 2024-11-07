test_that("peptide_quantification works correctly for PSM and Area", {

  # Define the input data
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

  matching_result <- data.frame(
    Sequence = c("AAA", "DDD", "DDD"),
    Condition_1 = c("Drug1", "Drug2", "Drug2"),
    Condition_2 = c("Donor1", "Donor2", "Donor2"),
    Region_1 = c("VH", "VL", "VL"),
    Region_2 = c("Arm_1", "Arm_2", "Arm_2"),
    Start_Position = c(4, 4, 4),
    End_Position = c(6, 6, 6),
    PTM_position = c(NA, 2, 0),
    PTM_type = c(NA, "O", "C"),
    Area = c(100, 200, 200),
    reps = c(1, 2, 2)
  )

  matching_columns <- c("Condition_1", "Region_2")
  distinct_columns <- c("Condition_2", "Region_1")
  area_column <- "Area"

  # Test for PSM quantification
  data_with_psm <- peptide_quantification(
    whole_seq,
    matching_result,
    matching_columns,
    distinct_columns,
    quantify_method = "PSM",
    with_PTM = TRUE,
    reps = TRUE
  )

  expect_equal(nrow(data_with_psm), 16 * 6)  # 16 sequences * 6 positions each
  expect_true("PSM" %in% colnames(data_with_psm))
  expect_true("PTM" %in% colnames(data_with_psm))
  expect_true("PTM_type" %in% colnames(data_with_psm))

  # Test for Area quantification
  data_with_area <- peptide_quantification(
    whole_seq,
    matching_result,
    matching_columns,
    distinct_columns,
    quantify_method = "Area",
    area_column = area_column,
    with_PTM = TRUE,
    reps = TRUE
  )

  expect_equal(nrow(data_with_area), 16 * 6)  # 16 sequences * 6 positions each
  expect_true("Area" %in% colnames(data_with_area))
  expect_true("PTM" %in% colnames(data_with_area))
  expect_true("PTM_type" %in% colnames(data_with_area))
})
