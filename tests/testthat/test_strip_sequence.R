test_that("strip_sequence works for PEAKS data", {
  data <- data.table(
    Peptide = c(
      "AAN(+0.98)Q(-0.98)RGSLYQCDYSTGSC(+57.02)EPIR",
      "K.AAQQTGKLVHANFGT.K",
      "K.(+0.98)AATVTGKLVHANFGT.K"
    ),
    Condition = c("A", "B", "B")
  )
  expected <- data.table(
    Peptide = c(
      "AAN(+0.98)Q(-0.98)RGSLYQCDYSTGSC(+57.02)EPIR",
      "K.AAQQTGKLVHANFGT.K",
      "K.(+0.98)AATVTGKLVHANFGT.K"
    ),
    Condition = c("A", "B", "B"),
    Sequence = c(
      "AANQRGSLYQCDYSTGSCEPIR",
      "AAQQTGKLVHANFGT",
      "AATVTGKLVHANFGT"
    )
  )
  result <- strip_sequence(data, "Peptide", "Sequence", "PEAKS")
  expect_equal(result, expected)
})

test_that("strip_sequence works for Spectronaut data", {
  data <- data.table(
    EG.IntPIMID = c(
      "_[+42]M[-16]DDREDLVYQAK_",
      "_EAAENSLVAYK_",
      "_IEAELQDIC[+57]NDVLELLDK_"
    ),
    Condition = c("A", "B", "B")
  )
  expected <- data.table(
    EG.IntPIMID = c(
      "_[+42]M[-16]DDREDLVYQAK_",
      "_EAAENSLVAYK_",
      "_IEAELQDIC[+57]NDVLELLDK_"
    ),
    Condition = c("A", "B", "B"),
    Sequence = c(
      "MDDREDLVYQAK",
      "EAAENSLVAYK",
      "IEAELQDICNDVLELLDK"
    )
  )
  result <- strip_sequence(data, "EG.IntPIMID", "Sequence", "Spectronaut")
  expect_equal(result, expected)
})

test_that("strip_sequence works for MSFragger data", {
  data <- data.table(
    'Modified Peptide' = c(
      "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
      "K.AAQQTGKLVHANFGT.K",
      "K.[0.98]AATVTGKLVHANFGT.K"
    ),
    Condition = c("A", "B", "B")
  )
  expected <- data.table(
    'Modified Peptide' = c(
      "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
      "K.AAQQTGKLVHANFGT.K",
      "K.[0.98]AATVTGKLVHANFGT.K"
    ),
    Condition = c("A", "B", "B"),
    Sequence = c(
      "AAMQRGSLYQCDYSTGSCEPIR",
      "AAQQTGKLVHANFGT",
      "AATVTGKLVHANFGT"
    )
  )
  result <- strip_sequence(data, 'Modified Peptide', 'Sequence', "MSFragger")
  expect_equal(result, expected)
})

test_that("strip_sequence works for Comet data", {
  data <- data.table(
    modified_peptide = c(
      "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
      "K.AAQQTGKLVHANFGT.K",
      "K.[0.98]AATVTGKLVHANFGT.K"
    ),
    Condition = c("A", "B", "B")
  )
  expected <- data.table(
    modified_peptide = c(
      "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
      "K.AAQQTGKLVHANFGT.K",
      "K.[0.98]AATVTGKLVHANFGT.K"
    ),
    Condition = c("A", "B", "B"),
    Sequence = c(
      "AAMQRGSLYQCDYSTGSCEPIR",
      "AAQQTGKLVHANFGT",
      "AATVTGKLVHANFGT"
    )
  )
  result <- strip_sequence(data, 'modified_peptide', 'Sequence', "Comet")
  expect_equal(result, expected)
})

test_that("strip_sequence works for DIANN data", {
  data <- data.table(
    Modified.Sequence = c(
      "AAAAGPGAALS(UniMod:21)PRPC(UniMod:4)DSDPATPGAQSPK",
      "AAAASAAEAGIATTGTEDSDDALLK",
      "AAAAALSGSPPQTEKPT(UniMod:21)HYR"
    ),
    Condition = c("A", "B", "B")
  )
  expected <- data.table(
    Modified.Sequence = c(
      "AAAAGPGAALS(UniMod:21)PRPC(UniMod:4)DSDPATPGAQSPK",
      "AAAASAAEAGIATTGTEDSDDALLK",
      "AAAAALSGSPPQTEKPT(UniMod:21)HYR"
    ),
    Condition = c("A", "B", "B"),
    Sequence = c(
      "AAAAGPGAALSPRPCDSDPATPGAQSPK",
      "AAAASAAEAGIATTGTEDSDDALLK",
      "AAAAALSGSPPQTEKPTHYR"
    )
  )
  result <- strip_sequence(data, 'Modified.Sequence', 'Sequence', "DIANN")
  expect_equal(result, expected)
})

test_that("strip_sequence works for Skyline data", {
  data <- data.table(
    'Peptide Modified Sequence' = c(
      "AGLC[+57]QTFVYGGC[+57]R",
      "AAAASAAEAGIATTGTEDSDDALLK",
      "IVGGWEC[+57]EK"
    ),
    Condition = c("A", "B", "B")
  )
  expected <- data.table(
    'Peptide Modified Sequence' = c(
      "AGLC[+57]QTFVYGGC[+57]R",
      "AAAASAAEAGIATTGTEDSDDALLK",
      "IVGGWEC[+57]EK"
    ),
    Condition = c("A", "B", "B"),
    Sequence = c(
      "AGLCQTFVYGGCR",
      "AAAASAAEAGIATTGTEDSDDALLK",
      "IVGGWECEK"
    )
  )
  result <- strip_sequence(data, 'Peptide Modified Sequence', 'Sequence', "Skyline")
  expect_equal(result, expected)
})

test_that("strip_sequence works for Maxquant data", {
  data <- data.table(
    'Modified sequence' = c(
      "_(ac)AA(ox)AAELRLLEK_",
      "_EAAENSLVAYK_",
      "_AADTIGYPVM(ox)IRSAYALGGLGSGICPNK_"
    ),
    Condition = c("A", "B", "B")
  )
  expected <- data.table(
    'Modified sequence' = c(
      "_(ac)AA(ox)AAELRLLEK_",
      "_EAAENSLVAYK_",
      "_AADTIGYPVM(ox)IRSAYALGGLGSGICPNK_"
    ),
    Condition = c("A", "B", "B"),
    Sequence = c(
      "AAAAELRLLEK",
      "EAAENSLVAYK",
      "AADTIGYPVMIRSAYALGGLGSGICPNK"
    )
  )
  result <- strip_sequence(data, 'Modified sequence', 'Sequence', "Maxquant")
  expect_equal(result, expected)
})
