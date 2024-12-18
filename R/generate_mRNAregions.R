#' Generate mRNA Regions from a GenBank Table
#'
#' This function generates mRNA regions by extracting upstream and downstream sequences
#' based on user-specified parameters. It supports dynamic upstream sequence generation,
#' custom downstream sequences, and circular genome handling.
#'
#' @param genbank_table A data frame containing genomic information such as `locus_tag`, `start`, `end`, and `direction`.
#'                      If `NULL`, the function will attempt to find a `genbank_table` object in the global environment.
#' @param upstream_values A numeric vector specifying the upstream lengths to extract. Default is `seq(100, 500, 100)`.
#' @param downstream A numeric value specifying the downstream length to extract. Default is `100`.
#' @param custom_downstream_seq A character string specifying a custom downstream sequence. If provided, the downstream
#'                              parameter is ignored, and the custom sequence is used.
#'
#' @return A data frame with additional columns for each upstream and downstream combination.
#' @export

generate_mRNA_regions <- function(
    genbank_table = NULL,
    upstream_values = seq(100, 500, 100),
    downstream = 100,
    custom_downstream_seq = NULL
) {
  start_time <- Sys.time()

  # Validate genbank_table
  if (is.null(genbank_table)) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      genbank_table <- get("genbank_table", envir = .GlobalEnv)
    } else {
      stop("genbank_table is not provided and could not be found in the environment.")
    }
  }

  # Ensure numeric types for start and end
  genbank_table$start <- as.numeric(genbank_table$start)
  genbank_table$end <- as.numeric(genbank_table$end)

  if (!is.null(custom_downstream_seq)) {
    if (!is.character(custom_downstream_seq) || length(custom_downstream_seq) != 1) {
      stop("custom_downstream_seq must be a single character string.")
    }
    downstream <- nchar(custom_downstream_seq)  # Override downstream length
  }

  # Identify circular contigs
  genome_summary <- get("Genome_summary", envir = .GlobalEnv)
  circular_contigs <- genome_summary %>%
    dplyr::filter(grepl("circular", trimws(tolower(LOCUS)), ignore.case = TRUE)) %>%
    dplyr::pull(LOCUS)

  circular_contigs <- unique(circular_contigs)

  # Preload contig sequences into a named list
  contig_sequences <- lapply(unique(genbank_table$contig_number), function(contig) {
    seq_name <- paste0("contig_", contig, "_seq")
    if (exists(seq_name, envir = .GlobalEnv)) {
      get(seq_name, envir = .GlobalEnv)
    } else {
      stop(paste("Contig sequence", seq_name, "not found in the global environment."))
    }
  })
  names(contig_sequences) <- unique(genbank_table$contig_number)

  # Define helper function for circular genome handling
  get_circular_segment <- function(dna, seq_len, base_pos, length_needed) {
    dna_double <- paste0(dna, dna)
    start_mod <- ((base_pos - 1) %% seq_len) + 1
    substring(dna_double, start_mod, start_mod + length_needed - 1)
  }

  # Helper function to extract regions
  extract_region <- function(row, up, down) {
    contig_seq <- contig_sequences[[row$contig_number]]
    seq_length <- nchar(contig_seq)
    is_circular <- row$contig_number %in% circular_contigs

    if (row$direction == "+") {
      # Forward strand
      cds_start_pos <- row$start

      # Upstream
      upstream_seq <- if (up > 0) {
        if (is_circular && (cds_start_pos - up < 1)) {
          missing_length <- abs(cds_start_pos - up) + 1
          paste0(
            substring(contig_seq, seq_length - missing_length + 1, seq_length),
            substring(contig_seq, 1, cds_start_pos - 1)
          )
        } else {
          substring(contig_seq, max(1, cds_start_pos - up), cds_start_pos - 1)
        }
      } else {
        ""
      }

      # Downstream
      downstream_seq <- if (is_circular) {
        get_circular_segment(contig_seq, seq_length, cds_start_pos, down)
      } else {
        substring(contig_seq, cds_start_pos, min(seq_length, cds_start_pos + down - 1))
      }

      # Combine upstream and downstream
      combined_seq <- paste0(upstream_seq, downstream_seq)
      combined_seq <- gsub("T", "U", toupper(combined_seq))  # Convert to RNA
      return(combined_seq)

    } else {
      # Reverse strand
      cds_end_pos <- row$end

      # Downstream (upstream for reverse strand)
      upstream_seq <- if (up > 0) {
        if (is_circular && (cds_end_pos + up > seq_length)) {
          overflow_length <- (cds_end_pos + up) - seq_length
          paste0(
            substring(contig_seq, cds_end_pos + 1, seq_length),
            substring(contig_seq, 1, overflow_length)
          )
        } else {
          substring(contig_seq, cds_end_pos + 1, min(seq_length, cds_end_pos + up))
        }
      } else {
        ""
      }

      # Upstream (downstream for reverse strand)
      downstream_seq <- if (down > 0) {
        if (is_circular && (cds_end_pos - down < 1)) {
          missing_length <- abs(cds_end_pos - down) + 1
          paste0(
            substring(contig_seq, seq_length - missing_length + 1, seq_length),
            substring(contig_seq, 1, cds_end_pos)
          )
        } else {
          substring(contig_seq, max(1, cds_end_pos - down), cds_end_pos)
        }
      } else {
        ""
      }

      # Combine upstream and downstream in correct order for reverse strand
      combined_seq <- paste0(downstream_seq, upstream_seq)

      # Reverse complement and convert to RNA
      combined_seq <- seqinr::c2s(rev(seqinr::comp(seqinr::s2c(combined_seq))))
      combined_seq <- gsub("T", "U", toupper(combined_seq))  # Convert to RNA
      return(combined_seq)
    }
  }

  # Generate sequences for each upstream value
  for (up in upstream_values) {
    col_name <- paste0("seq_u", up, "_d", ifelse(is.null(custom_downstream_seq), downstream, "GOI"))

    genbank_table[[col_name]] <- purrr::pmap_chr(
      genbank_table,
      function(contig_number, start, end, direction, ...) {
        extract_region(
          list(
            contig_number = contig_number,
            start = start,
            end = end,
            direction = direction
          ),
          up,
          downstream
        )
      }
    )
  }

  # Add custom downstream sequence column if provided
  if (!is.null(custom_downstream_seq)) {
    genbank_table$seq_custom_downstream <- custom_downstream_seq
  }

  assign("genbank_table", genbank_table, envir = .GlobalEnv)
  end_time <- Sys.time()
  cat("Execution Time: ", round(end_time - start_time, 2), "seconds\n")
}
