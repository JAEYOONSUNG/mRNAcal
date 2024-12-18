#' Generate Regions and RNA Folding Results
#'
#' This function extracts specific genomic regions based on user-defined upstream and downstream lengths, processes sequences through RNAfold, and visualizes RNA structures and energy distributions.
#'
#' @param genbank_table A data frame containing genomic data. If not provided, it searches for an object ending with "_total" in the global environment.
#' @param upstream Length of the upstream region to extract. Default is 500.
#' @param downstream Length of the downstream region to extract. Default is 100.
#' @param complete Logical, whether to handle circular genomes. Default is FALSE.
#' @param data A data frame containing sequence columns matching the regex pattern "seq_u[0-9]+_d[0-9]+$".
#' @param rnafold_path Path to the RNAfold executable. Default is "/Users/JaeYoon/miniconda3/bin/RNAfold".
#'
#' @return A processed data frame containing extracted regions, RNAfold results, and visualization data.
#'
#' @details
#' This function integrates genomic region extraction, RNAfold execution, and result visualization. It handles circular genomes, computes reverse complements for negative strands, and processes sequences in batches using RNAfold.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' genbank_table <- my_genbank_data
#' rnafold_path <- "/path/to/RNAfold"
#' result <- process_all_sequences(genbank_table, rnafold_path)

get_mRNA_region <- function(genbank_table = NULL, upstream = 100, downstream = 100, complete = FALSE) {
  # Find "_total" object if genbank_table is not provided
  if (is.null(genbank_table)) {
    total_objs <- ls(pattern = "_total$|genbank_table", envir = .GlobalEnv)
    if (length(total_objs) == 0) {
      stop("No object ending with '_total' found in the environment.")
    } else if (length(total_objs) > 1) {
      stop("Multiple objects ending with '_total' found. Please ensure only one such object exists.")
    }
    total_obj_name <- total_objs[1]
    genbank_table <- get(total_obj_name, envir = .GlobalEnv)
  }

  # Ensure contig sequence exists
  if (!exists("contig_name", envir = .GlobalEnv)) {
    genbank_fna_extractor()
    } else{next}

  contig_number_val <- genbank_table$contig_number[[1]]
  direction_val <- genbank_table$direction[[1]]
  start_val <- as.numeric(genbank_table$start[[1]])
  end_val <- as.numeric(genbank_table$end[[1]])

  contig_name <- paste0("contig_", contig_number_val, "_seq")
  seq_dna <- get(contig_name, envir = .GlobalEnv)
  seq_length <- nchar(seq_dna)
  cds_start_pos <- if (direction_val == "+") start_val else end_val

  get_circular_segment <- function(dna_double, seq_len, base_pos, length_needed) {
    start_mod <- ((base_pos - 1) %% seq_len) + 1
    end_mod <- start_mod + length_needed - 1
    substring(dna_double, start_mod, end_mod)
  }

  if (complete) {
    seq_dna_double <- paste0(seq_dna, seq_dna)
    upstream_seq <- if (upstream > 0) get_circular_segment(seq_dna_double, seq_length, cds_start_pos - upstream, upstream) else ""
    downstream_seq <- get_circular_segment(seq_dna_double, seq_length, cds_start_pos, downstream)
    combined_seq <- paste0(upstream_seq, downstream_seq)
  } else {
    region_start_up <- max(1, cds_start_pos - upstream)
    region_end_up <- max(0, cds_start_pos - 1)
    region_start_down <- cds_start_pos
    region_end_down <- min(seq_length, cds_start_pos + downstream)

    up_seq <- if (region_end_up >= region_start_up) substring(seq_dna, region_start_up, region_end_up) else ""
    down_seq <- if (region_start_down <= seq_length) substring(seq_dna, region_start_down, region_end_down) else ""
    combined_seq <- paste0(up_seq, down_seq)
  }

  if (direction_val == "-") {
    combined_seq <- seqinr::c2s(rev(seqinr::comp(seqinr::s2c(toupper(combined_seq)))))
  } else {
    combined_seq <- toupper(combined_seq)
  }

  return(combined_seq)
}
