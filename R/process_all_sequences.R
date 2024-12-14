#' Process RNAfold for Multiple Sequences
#'
#' This function processes multiple RNA sequences using RNAfold, extracts their structures
#' and minimum free energy (ΔG) values, and appends the results to the input data frame.
#'
#' @param data A data frame containing sequence columns that match the regex pattern "seq_u[0-9]+_d[0-9]+$".
#' @param rnafold_path A character string specifying the path to the RNAfold executable.
#'                     Default is "/Users/JaeYoon/miniconda3/bin/RNAfold".
#' @return A data frame with the original columns and appended RNAfold results, including:
#'   \describe{
#'     \item{Column_mRNA_Sequence}{The RNA sequence processed by RNAfold.}
#'     \item{Column_Structure}{The dot-bracket structure of the RNA.}
#'     \item{Column_FreeEnergy}{The minimum free energy (ΔG) value of the RNA structure.}
#'   }
#' @details
#' The function identifies sequence columns in the input data frame using the regex pattern
#' "seq_u[0-9]+_d[0-9]+$", applies RNAfold to each sequence, and parses the RNAfold output.
#' Results are stored as new columns in the data frame, with column names prefixed by the original
#' sequence column name.
#'
#' @examples
#' # Example usage:
#' data <- data.frame(
#'   seq_u100_d100 = c("GCGCUUCGCC", "UUCGGA"),
#'   seq_u200_d100 = c("AUGGCC", "GGUUACCCGG")
#' )
#' result <- process_all_sequences(data)
#' print(result)
#'
#' @export

process_all_sequences <- function(data, rnafold_path = "/Users/JaeYoon/miniconda3/bin/RNAfold") {
  # Extract the name of the input data object
  data_name <- deparse(substitute(data))

  # Identify columns matching the regex pattern for sequences
  sequence_cols <- grep("seq_u[0-9]+_d([0-9]+|GOI)$", colnames(data), value = TRUE)

  if (length(sequence_cols) == 0) stop("No matching columns found.")

  # Loop through each sequence column
  for (col in sequence_cols) {
    cat("Processing column:", col, "\n")  # Progress message

    # Apply run_rnafold to each sequence in the column
    results <- lapply(data[[col]], function(seq) {
      if (!is.na(seq) && nchar(seq) > 0) {  # Check for valid sequence
        as.data.frame(run_mRNAfold(seq, rnafold_path), stringsAsFactors = FALSE)
      } else {
        data.frame(mRNA_Sequence = NA, Structure = NA, FreeEnergy = NA, stringsAsFactors = FALSE)
      }
    })

    # Combine results into a single data frame
    results_df <- do.call(rbind, results) %>%
      dplyr::mutate(
        across(everything(), stringr::str_squish),  # Remove extra spaces
        FreeEnergy = as.numeric(FreeEnergy)  # Ensure FreeEnergy is numeric
      )

    # Update column names to include the original column name
    colnames(results_df) <- paste(col, colnames(results_df), sep = "_")

    # Append the results to the original data
    data <- cbind(data, results_df)
  }

  # Save the updated data to the global environment with the original name
  assign(data_name, data, envir = .GlobalEnv)
}
