
## A function for creating contingency tables
.contingencyTable <- function(set, reference, sig) {
    df <- data.frame(
        reference = reference,
        set = ifelse(reference %in% set, TRUE, FALSE),
        sig = ifelse(reference %in% sig, TRUE, FALSE)
    )
    df$set <- factor(df$set, c('TRUE', 'FALSE'), labels = c('inSet', 'notInSet'))
    df$sig <- factor(df$sig, c('TRUE', 'FALSE'), labels = c('inSig', 'notInSig'))
    stats::xtabs(formula = ~ set + sig, data = df)
}

#' Microbe Set Enrichment
#'
#' \code{microbeSetEnrichment} runs ORA (Fisher's exact test) comparing
#' a set of microbes (character vector), a reference (character vector),
#' and a list of signatures (list of vectors). Furthermore,
#' \code{microbeSetEnrichment} adjusts p values (FDR), and calculate odd ratios.
#'
#' @param set A character vector.
#' @param reference A character vector.
#' @param sigs A list of character vectors.
#'
#' @return A data.frame.
#' @export
#'
microbeSetEnrichment <- function(set, reference, sigs) {

    if (is.null(names(sigs)))
        stop('List of signatures must be named.', call. = FALSE)

    cols_order <- c(
        'sig_name', 'n_sig', 'n_background', 'p_value', 'fdr', 'odds_ratio',
        'upper_ci', 'lower_ci'
    )

    vct_list <- vector('list', length(sigs))
    for (i in seq_along(sigs)) {
        ct <- .contingencyTable(set, reference, sigs[[i]])
        n_sig <- ct[1]
        n_background <- ct[1] + ct[2]
        p_value <- stats::fisher.test(ct, alternative = 'g')$p.value

        odds_ratio <- suppressWarnings(
            epitools::oddsratio.wald(ct + 0.5)$measure[2,1]
        )
        upper_ci <- exp(log(odds_ratio) + 1.96 * sqrt(sum(1 / (ct + 1) )))
        lower_ci <- exp(log(odds_ratio) - 1.96 * sqrt(sum(1 / (ct + 1) )))

        vct_list[[i]] <- data.frame(
            n_sig = n_sig, n_background = n_background, p_value = p_value,
            odds_ratio = odds_ratio, upper_ci = upper_ci, lower_ci = lower_ci
        )
    }
    df <- do.call(rbind, vct_list)
    df$fdr <- stats::p.adjust(df$p_value, method = 'fdr')
    df$sig_name <- names(sigs)
    return(df[, cols_order])
}

#' Get Sets
#'
#' \code{getSets} gets sets of signatures.
#'
#' @param bsdb Data imported from bsdb.
#' @param bsdb_body_site Body site in bsdb.
#' @param tms Tidy TypicalMicrobiomeSignatures imported with
#' \code{importTypicalTMS}
#' @param tms_body_site Body site in tms.
#' @param tax_level Taxonomic rank
#'
#' @return A list of bsdb sets and background sets (bsdb + tms).
#' @export
#'
getSets <- function(
        bsdb, bsdb_body_site,
        tms, tms_body_site,
        tax_level
) {
    bsdb_sets <- bsdb |>
        dplyr::filter(
            `Body site` == bsdb_body_site, `Host species` == 'Homo sapiens'
        ) |>
        bugsigdbr::getSignatures(
            tax.id.type = 'ncbi', tax.level = tax_level, min.size = 5
        )
    typical_sets <- tms |>
        dplyr::filter(body_site == tms_body_site, rank == tax_level) |>
        dplyr::pull(taxid) |>
        unique()
    background_sets <- bsdb_sets |>
        lapply(function(x) unique(c(x, typical_sets)))
    output <- list(bsdb_sets = bsdb_sets, background_sets = background_sets)
    return(output)
}

#' Run enrichment
#'
#' \code{runEnrichment} runs the enrichment
#'
#' @param bsdb_sets Vector of bsdb sets (named list).
#' @param background_sets Vector of background sets (named list).
#' @param bp_sigs bugphyzz signatures (named list).
#' @param bsdb BSDB data.frame.
#' @param fdr_ths Threshold for filtering.
#'
#' @return A data.frame.
#' @export
#'
runEnrichment <- function(
    bsdb_sets, background_sets, bp_sigs, bsdb, fdr_ths = NULL
) {
    res <- purrr::map2(
        .x = bsdb_sets,
        .y = background_sets,
        .f = ~ microbeSetEnrichment(.x, .y, bp_sigs)
    ) |>
        dplyr::bind_rows(.id = 'bugsigdb_sig') |>
        dplyr::mutate(`BSDB ID` = sub('_.*$', '', bugsigdb_sig)) |>
        dplyr::relocate(`BSDB ID`) |>
        dplyr::rename(bp_sig = sig_name) |>
        tibble::as_tibble() |>
        dplyr::left_join(tibble::as_tibble(bsdb), by = 'BSDB ID')
    if (!is.null(fdr_ths)) {
        res <- dplyr::filter(res, fdr < fdr_ths)
    }
    return(res)
}



