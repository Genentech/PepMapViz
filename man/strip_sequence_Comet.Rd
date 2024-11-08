% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strip_sequence.R
\name{strip_sequence_Comet}
\alias{strip_sequence_Comet}
\title{Strip sequence from Comet outputs}
\usage{
strip_sequence_Comet(data, column, convert_column)
}
\arguments{
\item{data}{A dataframe with a column containing peptide sequences to be stripped}

\item{column}{The name of the column containing the peptide sequences to be stripped.}

\item{convert_column}{The name of the column where the stripped sequences will be stored.}
}
\value{
A dataframe with a column containing stripped sequence
}
\description{
This function takes Comet output containing a column with peptide sequences to be stripped
and converts it into a new dataframe with the stripped sequence
}
\examples{
library(data.table)
data <- data.table(
  modified_peptide = c(
    "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
    "K.AAQQTGKLVHANFGT.K",
    "K.[0.98]AATVTGKLVHANFGT.K"
  ),
  Condition = c("A", "B", "B")
)
column <- 'modified_peptide'
convert_column <- 'Sequence'
converted_data <- strip_sequence_Comet(data, column, convert_column)

}
