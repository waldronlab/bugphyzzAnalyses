
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
        'sig_name', 'n_sig', 'n_total', 'p_value', 'fdr', 'odds_ratio',
        'upper_ci', 'lower_ci'
    )

    vct_list <- vector('list', length(sigs))
    for (i in seq_along(sigs)) {
        ct <- .contingencyTable(set, reference, sigs[[i]])
        n_sig <- ct[1]
        n_total <- ct[1] + ct[2]
        p_value <- fisher.test(ct, alternative = 'g')$p.value

        odds_ratio <- suppressWarnings(
            epitools::oddsratio.wald(ct + 0.5)$measure[2,1]
        )
        upper_ci <- exp(log(odds_ratio) + 1.96 * sqrt(sum(1 / (ct + 1) )))
        lower_ci <- exp(log(odds_ratio) - 1.96 * sqrt(sum(1 / (ct + 1) )))

        vct_list[[i]] <- data.frame(
            n_sig = n_sig, n_total = n_total, p_value = p_value,
            odds_ratio = odds_ratio, upper_ci = upper_ci, lower_ci = lower_ci
        )
    }
    df <- do.call(rbind, vct_list)
    df$fdr <- stats::p.adjust(df$p_value, method = 'fdr')
    df$sig_name <- names(sigs)
    return(df[, cols_order])
}
