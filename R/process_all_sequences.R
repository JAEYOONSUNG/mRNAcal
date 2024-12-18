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

process_all_sequences <- function(data = NULL, rnafold_path = "/Users/JaeYoon/miniconda3/bin/RNAfold") {
  # If data is NULL, attempt to get genbank_table from the global environment
  if (is.null(data)) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      data <- get("genbank_table", envir = .GlobalEnv)
      message("Using 'genbank_table' from the global environment as the input data.")
    } else {
      stop("Input data is NULL, and 'genbank_table' does not exist in the global environment.")
    }
  }

  # Extract the name of the input data object
  data_name <- deparse(substitute(data))

  # Identify sequence columns using regex
  sequence_cols <- grep("seq_u[0-9]+_d([0-9]+|GOI)$", colnames(data), value = TRUE)
  if (length(sequence_cols) == 0) stop("No matching sequence columns found in the input data.")

  # Initialize a list to store results for all sequence columns
  all_results <- list()

  # Loop through each sequence column
  for (col in sequence_cols) {
    cat("Processing column:", col, "\n")  # Progress message

    # Initialize a data frame to store RNAfold results
    fold_results <- data.frame(
      locus_tag = data$locus_tag,  # Keep locus_tag as a key
      mRNA_Sequence = NA,
      Structure = NA,
      FreeEnergy = NA,
      stringsAsFactors = FALSE
    )

    # Apply run_mRNAfold to each sequence
    for (i in seq_along(data[[col]])) {
      seq <- data[[col]][i]
      locus_tag <- data$locus_tag[i]  # Get locus_tag for the current row

      if (!is.na(seq) && nchar(seq) > 0) {
        tryCatch({
          result <- run_mRNAfold(seq, rnafold_path)
          fold_results$mRNA_Sequence[i] <- result$mRNA_Sequence
          fold_results$Structure[i] <- result$Structure
          fold_results$FreeEnergy[i] <- result$FreeEnergy
        }, error = function(e) {
          warning(paste("RNAfold failed for locus_tag:", locus_tag,
                        "\nSequence:", seq,
                        "\nError Message:", e$message))
        })
      }
    }

    # Rename columns to reflect the original sequence column
    colnames(fold_results)[2:4] <- paste(col, c("mRNA_Sequence", "Structure", "FreeEnergy"), sep = "_")

    # Append results to the list
    all_results[[col]] <- fold_results
  }

  # Combine all RNAfold results into a single data frame using left join
  final_results <- data
  for (result_df in all_results) {
    final_results <- dplyr::left_join(final_results, result_df, by = "locus_tag")
  }

  # Save the updated data to the global environment
  assign(data_name, final_results, envir = .GlobalEnv)
  message("RNAfold processing complete. Updated data saved as:", data_name)
}
