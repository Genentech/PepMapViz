test_that("obtain_mod works correctly for different platforms", {

  # Test for Skyline
  data_skyline <- data.table(
    'Peptide Modified Sequence' = c(
      "AGLC[+57]QTFVYGGC[+57]R",
      "AAAASAAEAGIATTGTEDSDDALLK",
      "IVGGWEC[+57]EK"
    ),
    Condition = c("A", "B", "B")
  )
  PTM_table <- data.table(
    PTM_mass = c(57, -0.98, 15.9949),
    PTM_type = c("Cam", "Amid", "Ox")
  )
  converted_data_skyline <- obtain_mod(
    data_skyline,
    'Peptide Modified Sequence',
    'Skyline',
    seq_column = NULL,
    PTM_table,
    PTM_annotation = TRUE,
    PTM_mass_column = "PTM_mass"
  )

  expect_equal(nrow(converted_data_skyline), 4)
  expect_equal(converted_data_skyline$PTM_type, c(NA, "Cam", "Cam", "Cam"))

  # Test for Maxquant
  data_maxquant <- data.table(
    'Modified sequence' = c(
      "_(ac)AAAAELRLLEK_",
      "_EAAENSLVAYK_",
      "_AADTIGYPVM(ox)IRSAYALGGLGSGICPNK_"
    ),
    Condition = c("A", "B", "B")
  )
  PTM_table <- data.table(
    PTM_mass = c('ac', 'ox'),
    PTM_type = c("Acet", "Ox")
  )
  converted_data_maxquant <- obtain_mod(
    data_maxquant,
    'Modified sequence',
    'Maxquant',
    seq_column = NULL,
    PTM_table,
    PTM_annotation = TRUE,
    PTM_mass_column = "PTM_mass"
  )

  expect_equal(nrow(converted_data_maxquant), 3)
  expect_equal(converted_data_maxquant$PTM_type, c(NA, "Acet", "Ox"))

  # Test for PEAKS
  data_peaks <- data.table(
    Peptide = c(
      "AAN(+42)Q(-0.98)RGSLYQCDYSTGSC(+57.02)EPIR",
      "K.AAQQTGKLVHANFGT.K",
      "K.(-0.98)AATVTGKLVHANFGT.K"
    ),
    Sequence = c(
      "AANQRGSLYQCDYSTGSCEPIR",
      "AAQQTGKLVHANFGT",
      "AATVTGKLVHANFGT"
    ),
    Condition = c("A", "B", "B")
  )
  PTM_table <- data.table(PTM_mass = c(42, -0.98, 57.02),
                          PTM_type = c("Acet", "Amid", "Cam"))
  mod_column <- "Peptide"
  PTM_mass_column <- "PTM_mass"
  converted_data_peaks <- obtain_mod(
    data_peaks,
    mod_column,
    'PEAKS',
    seq_column = NULL,
    PTM_table,
    PTM_annotation = TRUE,
    PTM_mass_column
  )

  expect_equal(nrow(converted_data_peaks), 5)
  expect_equal(converted_data_peaks$PTM_type, c(NA, "Amid","Amid","Acet", "Cam"))
})
