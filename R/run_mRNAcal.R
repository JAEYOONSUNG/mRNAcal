#' Run the mRNA Analysis Pipeline
#'
#' This function executes a series of steps in the GenBank analysis pipeline, including
#' organizing GenBank data, extracting sequences, and processing RNAfold results.
#' Progress messages are displayed at each step to provide updates.
#'
#' @param genbank_file A character string specifying the path to the GenBank file. Default is `NULL`.
#'                     If `NULL`, the function uses the preloaded `genbank_table` in the environment.
#' @param rnafold_path A character string specifying the path to the RNAfold executable. Default is
#'                     `"/Users/JaeYoon/miniconda3/bin/RNAfold"`.
#' @param upstream_values A numeric vector specifying the upstream lengths to extract for mRNA regions.
#'                        Default is `seq(100, 500, 100)`.
#' @param downstream A numeric value specifying the downstream length to extract for mRNA regions.
#'                   Default is `100`.
#' @param custom_downstream_seq A character string specifying a custom downstream sequence for mRNA regions.
#'                              If provided, it overrides the downstream parameter. Default is `NULL`.
#' @return The processed data frame with RNAfold results appended.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Organizes GenBank data using `Genbank_organizer`.
#'   \item Extracts FASTA sequences with `genbank_fna_extractor`.
#'   \item Retrieves basic GenBank information using `gb_info`.
#'   \item Extracts mRNA regions with `generate_mRNA_regions`.
#'   \item Processes RNAfold results using `process_all_sequences`.
#' }
#' Progress messages are displayed at each step to inform the user of the pipeline status.
#'
#' @examples
#' # Example usage
#' run_mRNAcal(
#'   genbank_file = "example.gbk",
#'   rnafold_path = "/path/to/RNAfold",
#'   upstream_values = seq(100, 500, 100),
#'   downstream = 100,
#'   custom_downstream_seq = NULL
#' )
#'
#' @export

run_mRNAcal <- function(
    genbank_file = NULL,
    rnafold_path = "/Users/JaeYoon/miniconda3/bin/RNAfold",
    upstream_values = seq(500, 500, 100),
    downstream = 100,
    custom_downstream_seq = NULL
) {
  # Check if RNAfold path exists
  if (!file.exists(rnafold_path)) {
    stop(paste("RNAfold executable not found at:", rnafold_path,
               "\nPlease provide a valid path to the RNAfold executable."))
  }
  message("RNAfold executable found at:", rnafold_path)

  message("Starting the mRNA Analysis Pipeline...")

  # Step 1: Organize GenBank data
  message("Step 1: Organizing GenBank data...")

  # Check if genbank_table exists; if not, run Genbank_organizer()
  if (!exists("genbank_table", envir = .GlobalEnv) || is.null(get("genbank_table", envir = .GlobalEnv))) {
    message("genbank_table not found in the environment. Running Genbank_organizer()...")
    Genbank_organizer()

    # Verify that Genbank_organizer() successfully created genbank_table
    if (!exists("genbank_table", envir = .GlobalEnv) || is.null(get("genbank_table", envir = .GlobalEnv))) {
      stop("Genbank_organizer() did not create a valid genbank_table.")
    }
  } else {
    message("Using pre-existing genbank_table in the environment.")
  }
  message("GenBank data organization complete.")

  # Step 2: Extract FASTA sequences
  message("Step 2: Extracting FASTA sequences...")
  genbank_fna_extractor()
  message("FASTA sequence extraction complete.")

  # Step 3: Retrieve basic GenBank information
  message("Step 3: Retrieving GenBank information...")
  gb_info()
  message("GenBank information retrieval complete.")

  # Step 4: Generate mRNA regions
  message("Step 4: Generating mRNA regions...")
  generate_mRNA_regions(
    upstream_values = upstream_values,
    downstream = downstream,
    custom_downstream_seq = custom_downstream_seq
  )
  message("mRNA region generation complete.")

  # Step 5: Process RNAfold results
  if (exists("genbank_table", envir = .GlobalEnv)) {
    genbank_table <- get("genbank_table", envir = .GlobalEnv)
    process_all_sequences(rnafold_path = rnafold_path)
    assign("genbank_table", genbank_table, envir = .GlobalEnv)
  } else {
    stop("genbank_table not found in the global environment.")
  }

  message("RNAfold processing complete.")

  message("Pipeline execution completed successfully!")
}
