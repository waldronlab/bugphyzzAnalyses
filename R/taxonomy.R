
## Functions for dealing with the NCBI taxonomy

#' Valid ranks
#'
#' \code{validRanks} get valid ranks.
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

#' Get taxonomy classification
#'
#' \code{.getTaxonClassification} Gets the taxonomy clasification using
#' the `taxize:classification` function of the first non-NA element of an
#' input vector of taxids.
#'
#' @param x A vector of taxonomy names or taxids.
#'
#' @return A data frame with taxonomic classification (as in the taxize format).
#' @keywords internal
#'
.getTaxonClassification <- function(x) {
    if (all(is.na(x))) {

        output <- data.frame(
            name = x,
            rank = names(x),
            id = x
        )
        return(output)
    }

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

#' Taxize classification to data frame
#'
#' \code{classif2Table} converts from classification format (`taxize` packatge)
#' to a data frame.
#'
#' @param x Output from the `taxize::classification` function.
#' @param ranks Ranks to be selected. Default is given by `validRanks()`.
#'
#' @return A data frame.
#' @export
#'
classif2Table <- function(x, ranks) {

    if (missing(ranks)) {
        valid_ranks <- validRanks()
    } else {
        valid_ranks <- ranks
    }

    df_filtered <- x |>
        dplyr::select(rank, id) |>
        dplyr::filter(rank %in% valid_ranks)

    new_df <- data.frame(x = df_filtered$id) |>
        t() |>
        as.data.frame(row.names = 1L)
    colnames(new_df) <- df_filtered$rank
    new_df
}

#' Convert taxonomic names to NCBI taxids
#'
#' \code{taxTable2Taxid} converts a table of taxonomic names to a table of
#' taxids. The format must be the same as curatedMetagenomicData (rowData) and
#' phyloseq (taxa_data).
#'
#' This function should always be run in an interactive session because it will
#' ask for clarification when two taxonomic names belong to two different
#' taxids.
#'
#' @param df A data frame.
#' @param names_from Name of the column with the names that should be used as
#' rownames.
#'
#' @return A data frame with taxIDs instead of names.
#' @export
#'
taxTable2taxid <- function(df, names_from) {

    colnames(df) <- tolower(colnames(df))
    names_from <- tolower(names_from)

    if (missing(names_from))
        stop('Missing argument "names_from".', call. = FALSE)

    if(!names_from %in% colnames(df))
        stop(names_from, ' was not found in the data.', call. = FALSE)

    col_names <- colnames(df)[which(colnames(df) %in% validRanks())]
    row_names <- df[[names_from]]

    tax_list <- df |>
        dplyr::select(tidyselect::all_of(col_names)) |>
        t() |>
        as.data.frame() |>
        as.list() |>
        magrittr::set_names(row_names) |>
        lapply(
            function(.x) {
                names(.x) <- col_names
                .x
            }
        )

    tax_list_unique <- unique(tax_list)
    names(tax_list_unique) <- vapply(
        tax_list_unique, function(x) paste0(x, collapse = '|'), character(1)
    )

    full_taxon_name <- vapply(
        tax_list, function(x) paste0(x, collapse = '|'), character(1),
        USE.NAMES = FALSE
    )

    row_data_subset <- df |>
        dplyr::select(-tidyselect::all_of(col_names)) |>
        dplyr::mutate(full_taxon_name = full_taxon_name)

    new_row_data <- tax_list_unique |>
        lapply(.getTaxonClassification) |>
        lapply(classif2Table) |>
        dplyr::bind_rows(.id = 'full_taxon_name')

    output <-
        dplyr::left_join(
            row_data_subset, new_row_data, by = 'full_taxon_name'
        ) |>
        as.data.frame()
    # new_row_data <- as.data.frame(new_row_data)

    if ('superkingdom' %in% colnames(output)) {
        n_pos <- which(colnames(output) == 'superkingdom')
        colnames(output)[n_pos] <- 'kingdom'
    }
    return(output)
}
