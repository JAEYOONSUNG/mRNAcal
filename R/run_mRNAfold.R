#' Execute mRNAfold on a Given RNA Sequence
#'
#' This function runs RNAfold on a given RNA sequence and extracts the dot-bracket structure
#' and minimum free energy (ΔG) from the RNAfold output.
#'
#' @param sequence A character string representing the RNA sequence to fold.
#' @param rnafold_path A character string specifying the path to the RNAfold executable.
#'                     Default is "/usr/local/bin/RNAfold".
#' @return A list containing the following components:
#'   \describe{
#'     \item{mRNA_Sequence}{The RNA sequence as processed by RNAfold.}
#'     \item{Structure}{The dot-bracket structure of the RNA.}
#'     \item{FreeEnergy}{The minimum free energy (ΔG) value of the RNA structure.}
#'   }
#' @details
#' The function assumes that RNAfold is installed and available at the specified path.
#' It captures the RNAfold output and parses the dot-bracket structure and ΔG value.
#' If RNAfold output is malformed or incomplete, the function throws an error.
#'
#' @examples
#' # Example usage:
#' sequence <- "GCGCUUCGCCGAAGCGCUUCGCCGAAGCGCUUCGCCGAAGCGCUUCGCCGAAGC"
#' result <- run_rnafold(sequence)
#' print(result)
#'
#' @export

run_mRNAfold <- function(sequence, rnafold_path="/Users/JaeYoon/miniconda3/bin/RNAfold") {
  result <- system(paste(rnafold_path), input = sequence, intern = TRUE)
  if (length(result) < 2) stop("RNAfold output is malformed.")
  
  mRNA_sequnece <- result[1]
  structure_line <- result[2]
  free_energy <- as.numeric(gsub("[()]", "", sub(".*\\(([^)]+)\\)$", "\\1", structure_line)))
  structure <- sub("\\s*\\([^)]*\\)$", "", structure_line)
  
  list(
    mRNA_Sequence = mRNA_sequnece,
    Structure = structure,
    FreeEnergy = free_energy
  )
}
