% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_modification.R
\name{obtain_mod_mzTab}
\alias{obtain_mod_mzTab}
\title{Obtain modification information from Peptide data generated by mzTab}
\usage{
obtain_mod_mzTab(
  data,
  mod_column,
  seq_column,
  PTM_table = NULL,
  PTM_annotation = FALSE,
  PTM_mass_column
)
}
\arguments{
\item{data}{A data.table with a column containing stripped sequence and a column containing PTM information.}

\item{mod_column}{The name of the column containing the modified peptide sequences.}

\item{seq_column}{The name of the column containing peptide sequences for mzTab}

\item{PTM_table}{A data.table with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.}

\item{PTM_annotation}{A logical value indicating whether to include PTM annotation information in the result.}

\item{PTM_mass_column}{The name of the column containing the PTM mass information}
}
\value{
A data.table with 'PTM_position', 'PTM_type', 'reps', and other columns.
}
\description{
This function takes mzTab output containing a 'modifications' column with
PTM information and converts it into a new dataframe with the desired format of peptide
sequences and associated PTM information.
}
\examples{
library(data.table)
data <- data.table(
  sequence = c("DDREDMLVYQAK", "EAAENSLVAYK", "IEAELQDICNDVLELLDK"),
  modifications = c("4-UNIMOD:7,10-UNIMOD:35", NA, "8-UNIMOD:7"),
  Condition1 = c("A", "B", "B"),
  Condition2 = c("C", "C", "D")
)
PTM_table <- data.table(
  PTM_mass = c("UNIMOD:7", "UNIMOD:35"),
  PTM_type = c("Amid", "Ox")
)
mod_column <- "modifications"
seq_column <- "sequence"
converted_data <- obtain_mod_mzTab(
  data,
  mod_column,
  seq_column,
  PTM_table,
  PTM_annotation = TRUE,
  PTM_mass_column = "PTM_mass"
)

}
