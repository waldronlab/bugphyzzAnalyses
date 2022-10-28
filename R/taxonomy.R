
## Functions for dealing with the NCBI taxonomy

#' Get taxonomy classification
#'
#' @param x A vector of taxonomy names or taxids.
#'
#' @return The classification of the first non-NA element of the input vector.
#' @keywords internal
#'
.getTaxonClassification <- function(x) {
    x <- x[length(x):1]
    output <- vector('list', length(x))
    for (i in seq_along(output)) {
        if (is.na(x[i])) next
        res <- tryCatch(
            {
                taxize::classification(x[i], db = 'ncbi')
            },
            error = function(e) NA
        )
        if (all(is.na(res[[1]])) || all(is.null(res[[1]]))) {
            next
        } else {
            names(output)[i] <- x[i]
            output[[i]] <- res[[1]]
            break
        }
    }
    output <- output[which(!is.na(names(output)))]
    output[[1]]
}

#' Valid ranks
#'
#' @return A character vector.
#' @export
#'
validRanks <- function() {
    c(
        'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family',
        'genus', 'species', 'strain'
    )
}

#' Taxonomic names to taxIDS
#'
#' @param df A data frame.
#' @param names_from Name of the column with the names that should be used as
#' rownames. Default is to extract the rownames from the dataframe.
#'
#' @return A data frame with taxIDs instead of names.
#' @export
#'
taxNames2TaxIDs <- function(df, names_from) {

    colnames(df) <- tolower(colnames(df))
    names_from <- tolower(names_from)

    if (missing(names_from))
        stop('Missing argument "names_from"', call. = FALSE)

    if(!names_from %in% colnames(df))
        stop(names_from, ' was not found in the data.', call. = FALSE)

    col_names <- colnames(df)[which(colnames(df) %in% validRanks())]
    names <- df[[names_from]]

    tax_list <- df |>
        dplyr::select(tidyselect::all_of(col_names)) |>
        t() |>
        as.data.frame() |>
        as.list() |>
        magrittr::set_names(names)

    taxa_list_unique <- unique(tax_list)
    names(taxa_list_unique) <- vapply(
        taxa_list_unique, function(x) paste0(x, collapse = '|'), character(1)
    )

    list(
        row_names = names,
        key = vapply(
            tax_list, function(x) paste0(x, collapse = '|'), character(1),
            USE.NAMES = FALSE
        ),
        taxa_list = taxa_list_unique
    )

}





