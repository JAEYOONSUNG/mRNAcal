% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mRNAregion.R
\name{get_mRNA_region}
\alias{get_mRNA_region}
\title{Generate Regions and RNA Folding Results}
\usage{
get_mRNA_region(
  genbank_table = NULL,
  upstream = 500,
  downstream = 100,
  complete = FALSE
)
}
\arguments{
\item{genbank_table}{A data frame containing genomic data. If not provided, it searches for an object ending with "_total" in the global environment.}

\item{upstream}{Length of the upstream region to extract. Default is 500.}

\item{downstream}{Length of the downstream region to extract. Default is 100.}

\item{complete}{Logical, whether to handle circular genomes. Default is FALSE.}

\item{data}{A data frame containing sequence columns matching the regex pattern "seq_u\link{0-9}+_d\link{0-9}+$".}

\item{rnafold_path}{Path to the RNAfold executable. Default is "/Users/JaeYoon/miniconda3/bin/RNAfold".}
}
\value{
A processed data frame containing extracted regions, RNAfold results, and visualization data.
}
\description{
This function extracts specific genomic regions based on user-defined upstream and downstream lengths, processes sequences through RNAfold, and visualizes RNA structures and energy distributions.
}
\details{
This function integrates genomic region extraction, RNAfold execution, and result visualization. It handles circular genomes, computes reverse complements for negative strands, and processes sequences in batches using RNAfold.
}
\examples{
# Example usage:
genbank_table <- my_genbank_data
rnafold_path <- "/path/to/RNAfold"
result <- process_all_sequences(genbank_table, rnafold_path)
}
