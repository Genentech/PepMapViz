% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_modification.R
\name{obtain_mod_DIANN}
\alias{obtain_mod_DIANN}
\title{Obtain modification information from Peptide data generated by DIA-NN}
\usage{
obtain_mod_DIANN(
  data,
  mod_column,
  PTM_table = NULL,
  PTM_annotation = FALSE,
  PTM_mass_column
)
}
\arguments{
\item{data}{A dataframe with 'Stripped.Sequence' column and 'Modified.Sequence' column containing modified peptide sequences.}

\item{mod_column}{The name of the column containing the modified peptide sequences.}

\item{PTM_table}{A dataframe with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.}

\item{PTM_annotation}{A logical value indicating whether to include PTM annotation information in the result.}

\item{PTM_mass_column}{The name of the column containing the PTM mass information}
}
\value{
A dataframe with 'Peptide', 'PTM_position', and 'PTM_type' columns.
}
\description{
This function takes DIA-NN output containing a column with modified peptide
sequences including PTM information and converts it into a new dataframe with the
desired format of peptide sequences and associated PTM information.
}
\examples{
library(data.table)
data <- data.table(
  Modified.Sequence = c(
    "AAAAGPGAALS(UniMod:21)PRPC(UniMod:4)DSDPATPGAQSPK",
    "AAAASAAEAGIATTGTEDSDDALLK",
    "AAAAALSGSPPQTEKPT(UniMod:21)HYR"
  ),
  Stripped.Sequence = c(
    "AAAAGPGAALSPRPCDSDPATPGAQSPK",
    "AAAASAAEAGIATTGTEDSDDALLK",
    "AAAAALSGSPPQTEKPTHYR"
  ),
  Condition = c("A", "B", "B")
)
PTM_table <- data.table(PTM_mass = c('UniMod:21', 'UniMod:4'),
                        PTM_type = c("Phos", "Cam"))
converted_data <- obtain_mod_DIANN(
  data,
  'Modified.Sequence',
  PTM_table,
  PTM_annotation = TRUE,
  PTM_mass_column = "PTM_mass"
)

}
