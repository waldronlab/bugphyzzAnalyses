#' Get counts for attribute
#'
#' \code{get_attr_counts} gets counts of attributes per rank.
#'
#' @param x Attribute dataframe imported from bugphyzz.
#'
#' @return A vector of counts per ranks.
#' @importFrom magrittr %>%
#' @export
#'
get_attr_counts <- function(x) {

    NCBI_ID <- Taxon_name <- Rank <- NULL

    select_cols <- c('NCBI_ID', 'Taxon_name', 'Rank')

    if( !any(select_cols %in% colnames(x)))
        stop('Columns missing.', call. = FALSE)

    valid_ranks <- c(
        'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
        'species', 'strain'
    )

    attr_counts_df <- x %>%
        dplyr::select(NCBI_ID, Taxon_name, Rank) %>%
        dplyr::distinct() %>%
        dplyr::count(Rank)

    attr_counts <- attr_counts_df$n
    names(attr_counts) <- attr_counts_df$Rank
    attr_counts[which(names(attr_counts) %in% valid_ranks)]

}
