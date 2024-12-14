#' Generate mRNA Regions from a GenBank Table
#'
#' This function generates mRNA regions by extracting upstream and downstream sequences
#' based on user-specified parameters. It supports dynamic upstream sequence generation
#' and custom downstream sequences.
#'
#' @param genbank_table A data frame containing genomic information such as `locus_tag`, `start`, `end`, and `strand`.
#'                      If `NULL`, the function will attempt to find a `genbank_table` object in the global environment.
#' @param upstream_values A numeric vector specifying the upstream lengths to extract. Default is `seq(100, 1000, 100)`.
#' @param downstream A numeric value specifying the downstream length to extract. Default is `100`.
#' @param custom_downstream_seq A character string specifying a custom downstream sequence. If provided, the downstream
#'                              parameter is ignored, and the custom sequence is used.
#'
#' @return A data frame with additional columns for each upstream and downstream combination. Columns are named as:
#' \describe{
#'   \item{`seq_u[upstream]_d[downstream]`}{Generated sequences for upstream and downstream lengths.}
#'   \item{`seq_u[upstream]_d_GOI`}{Generated sequences with the custom downstream sequence (if provided).}
#'   \item{`seq_custom_downstream`}{The custom downstream sequence for all rows (if provided).}
#' }
#'
#' @details
#' This function processes each row of the input `genbank_table` to extract mRNA sequences
#' using the `get_mRNA_region` function. It dynamically generates column names based on
#' the specified upstream and downstream lengths. If a custom downstream sequence is provided,
#' it overrides the downstream length and generates columns with the `_GOI` suffix.
#'
#' @examples
#' # Example GenBank table
#' genbank_table <- data.frame(
#'   locus_tag = c("gene1", "gene2"),
#'   start = c(100, 200),
#'   end = c(500, 600),
#'   strand = c("+", "-")
#' )
#'
#' # Generate regions with default downstream length
#' result <- generate_mRNA_regions(genbank_table)
#'
#' # Generate regions with a custom downstream sequence
#' custom_result <- generate_mRNA_regions(
#'   genbank_table,
#'   upstream_values = seq(100, 500, 100),
#'   custom_downstream_seq = "AGGCTTAA"
#' )
#'
#' @export

generate_mRNA_regions <- function(genbank_table = NULL, upstream_values = seq(100, 500, 100), downstream = 100, custom_downstream_seq = NULL) {
  # If genbank_table is NULL, find it in the environment
  if (is.null(genbank_table)) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      genbank_table <- get("genbank_table", envir = .GlobalEnv)
    } else {
      stop("genbank_table is not provided and could not be found in the environment.")
    }
  }

  # Validate custom_downstream_seq input
  if (!is.null(custom_downstream_seq)) {
    if (!is.character(custom_downstream_seq) || length(custom_downstream_seq) != 1) {
      stop("custom_downstream_seq must be a single character string.")
    }
  }

  # Process each upstream value
  for (up in upstream_values) {
    col_name <- paste0("seq_u", up, "_d", ifelse(is.null(custom_downstream_seq), downstream, "GOI"))
    genbank_table <- genbank_table %>%
      dplyr::rowwise() %>%
      dplyr::mutate(!!col_name := get_mRNA_region(
        genbank_table = dplyr::cur_data(),
        upstream = up,
        downstream = if (is.null(custom_downstream_seq)) downstream else nchar(custom_downstream_seq)
      )) %>%
      dplyr::ungroup()
  }

  # Add custom downstream sequence if provided
  if (!is.null(custom_downstream_seq)) {
    genbank_table <- genbank_table %>%
      dplyr::mutate(seq_custom_downstream = custom_downstream_seq)
  }
  return(genbank_table)
}
