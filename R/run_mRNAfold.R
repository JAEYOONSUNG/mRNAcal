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

run_mRNAfold <- function(sequence, rnafold_path = "/Users/JaeYoon/miniconda3/bin/RNAfold") {
  # Load genbank_table from the global environment
  if (!exists("genbank_table", envir = .GlobalEnv)) {
    stop("genbank_table not found in the global environment.")
  }
  genbank_table <- get("genbank_table", envir = .GlobalEnv)

  # Find the corresponding locus_tag for the given sequence
  locus_tag <- genbank_table %>%
    dplyr::filter(dplyr::if_any(dplyr::everything(), ~ . == sequence)) %>%
    dplyr::pull(locus_tag)

  if (length(locus_tag) == 0) {
    locus_tag <- "Unknown"
  }

  tryCatch({
    # Run RNAfold
    result <- system(paste(rnafold_path), input = sequence, intern = TRUE)

    # Check for errors in RNAfold output
    if (length(result) < 2) {
      warning(paste("RNAfold output is malformed for locus_tag:", locus_tag,
                    "\nSequence:", sequence,
                    "\nRNAfold output:", paste(result, collapse = "\n")))
      # Return NA for malformed output
      return(data.frame(
        mRNA_Sequence = sequence,
        Structure = NA,
        FreeEnergy = NA,
        locus_tag = locus_tag,
        stringsAsFactors = FALSE
      ))
    }

    # Handle multiple results by taking the first two lines
    if (length(result) > 2) {
      warning(paste("Multiple results detected for locus_tag:", locus_tag,
                    "\nUsing only the first result."))
      result <- result[1:2]  # Use only the first two lines
    }

    # Parse RNAfold results
    mRNA_sequence <- result[1]
    structure_line <- result[2]
    free_energy <- as.numeric(gsub("[()]", "", sub(".*\\(([^)]+)\\)$", "\\1", structure_line)))
    structure <- sub("\\s*\\([^)]*\\)$", "", structure_line)

    # Return parsed results
    data.frame(
      mRNA_Sequence = mRNA_sequence,
      Structure = structure,
      FreeEnergy = free_energy,
      locus_tag = locus_tag,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    # Capture and display the full error context
    warning(paste("Error processing RNAfold for locus_tag:", locus_tag,
                  "\nSequence:", sequence,
                  "\nError Message:", e$message))
    # Return NA for failed sequences
    data.frame(
      mRNA_Sequence = sequence,
      Structure = NA,
      FreeEnergy = NA,
      locus_tag = locus_tag,
      stringsAsFactors = FALSE
    )
  })
}
