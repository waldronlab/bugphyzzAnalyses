
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
        'sig_name', 'n_set_annotated', 'n_background_annotated',
        "set_size", "background_size",
        'p_value', 'fdr', 'odds_ratio',
        'upper_ci', 'lower_ci'
    )

    vct_list <- vector('list', length(sigs))
    for (i in seq_along(sigs)) {
        ct <- .contingencyTable(set, reference, sigs[[i]])
        n_set_annotated <- ct[1]
        n_background_annotated <- ct[1] + ct[2]
        set_size <- ct[1] + ct[3]
        background_size <- ct[1] + ct[2] + ct[3] + ct[4]

        p_value <- stats::fisher.test(ct, alternative = 'g')$p.value
        res <- suppressWarnings(
            epitools::oddsratio.wald(ct + 0.5)$measure[2,]
        )
        odds_ratio <- res[["estimate"]]
        lower_ci <- res[["lower"]]
        upper_ci <- res[["upper"]]

        # upper_ci <- exp(log(odds_ratio) + 1.96 * sqrt(sum(1 / (ct + 0.5) )))
        # lower_ci <- exp(log(odds_ratio) - 1.96 * sqrt(sum(1 / (ct + 0.5) )))

        vct_list[[i]] <- data.frame(
            n_set_annotated = n_set_annotated,
            n_background_annotated = n_background_annotated,
            set_size = set_size,
            background_size = background_size,
            p_value = p_value,
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

#' Plot Enrichment Results
#'
#' \code{plotEnrichmentRes} plots the enrichment results.
#'
#' @param res Output from \code{runEnrichment}
#' @param body_site A character string for the title.
#' @param tax_level A character string for the title.
#'
#' @return A ggplot2 object.
#' @export
#'
plotEnrichmentRes <- function(res, body_site, tax_level) {
    dat <- res |>
        dplyr::filter(!is.na(`Abundance in Group 1`))

    if (!nrow(dat)) {
        warning('No results')
        return(NULL)
    }
    p <- dat |>
        dplyr::mutate(
            bugsigdb_sig = sub('^.*_(.*)_vs_.*$', '\\1', bugsigdb_sig),
            bp_sig = sub('^bugphyzz:', '', bp_sig),
            `Abundance in Group 1` = dplyr::case_when(
                `Abundance in Group 1` == 'increased' ~ 'Increased abundance in group 1',
                `Abundance in Group 1` == 'decreased' ~ 'Decreased abundance in group 1',
                TRUE ~ `Abundance in Group 1`
            ),
            `Abundance in Group 1` = forcats::fct_relevel(
                `Abundance in Group 1`, 'Increased abundance in group 1'
            )
        ) |>
        ggplot2::ggplot(ggplot2::aes(bp_sig, reorder(bugsigdb_sig, -log10(p_value)))) +
        ggplot2::geom_point(
            ggplot2::aes(
                color = -log10(p_value),
                size = n_sig / n_background
            )
        ) +
        ggplot2::scale_color_viridis_c(name = '-log10(p-value)', option = 'C') +
        ggplot2::scale_size(
            name = expression(
                frac('# annotated BugSigDB', '# annotated background')
            )
        ) +
        ggplot2::facet_wrap(~`Abundance in Group 1`, nrow = 2, scales = 'free_y') +
        ggplot2::labs(
            x = 'bugphyzz sets',
            y = 'BugSigDB sets - Group 1 name',
            title = 'Microbe set enrichment (ORA/Fisher\'s exact test)',
            subtitle = paste0(
                'Hypothesis: "greater"; FDR < 0.1;',
                body_site, '; ',  tax_level
            )
        ) +
        ggplot2::theme_bw()  +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    return(p)
}


