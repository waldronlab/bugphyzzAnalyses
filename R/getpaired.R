
#' Get paired signatures
#'
#' @param dec Signatures decreased.
#' @param inc Signatures increased.
#' @param rank Taxonomic rank.
#' @param cat if TRUE concatenate signatures. Default is FALSE.
#'
#' @return A list.
#' @export
#'
getPairedSigs <- function(dec, inc, rank, cat = FALSE) {
    sigs_inc <- getSignatures(
        inc, tax.id.type = "ncbi", tax.level = rank, min.size = 1,
        exact.tax.level = TRUE # Default
    )

    sig_names_inc <- names(sigs_inc)
    names(sig_names_inc) <- sub("^(bsdb:\\d+/\\d+)/\\d+_.*", "\\1", names(sigs_inc))
    # names(sigs_inc) <- sub("^(bsdb:\\d+/\\d+)/\\d+_.*", "\\1", names(sigs_inc))

    sigs_dec <- getSignatures(
        dec, tax.id.type = "ncbi", tax.level = rank, min.size = 1,
        exact.tax.level = TRUE  # Default
    )
    sig_names_dec <- names(sigs_dec)
    names(sig_names_dec) <- sub("^(bsdb:\\d+/\\d+)/\\d+_.*", "\\1", names(sigs_dec))
    # names(sigs_dec) <- sub("^(bsdb:\\d+/\\d+)/\\d+_.*", "\\1", names(sigs_dec))

    common_exps <- intersect(names(sig_names_dec), names(sig_names_inc))

    sig_names_dec <- sig_names_dec[common_exps]
    sig_names_inc <- sig_names_inc[common_exps]

    sigs_dec <- sigs_dec[sig_names_dec]
    sigs_inc <- sigs_inc[sig_names_inc]

    # common_sigs <- intersect(names(sigs_inc), names(sigs_dec))
    # sigs_inc <- sigs_inc[common_sigs]
    # sigs_dec <- sigs_dec[common_sigs]
    # list(dec = sigs_dec, inc = sigs_inc)

    if (cat) {
        output <- catSigs(sigs_dec, sigs_inc)
    } else {
        output <- list(dec = sigs_dec, inc = sigs_inc)
    }
    return(output)
}

catSigs <- function(dec, inc) {
    sigs_dec <- unlist(dec, use.names = FALSE)
    attr(sigs_dec, "nexp") <- length(dec)
    attr(sigs_dec, "signames") <- names(inc)
    # attr(sigs_dec, "signames") <- sub("^(.*\\d+)_.+", "\\1", names(sigs_dec))
    # attr(sigs_dec, "signames") <- sub("^(.*\\d+)_.+", "\\1", names(sigs_dec))

    sigs_inc <- unlist(inc, use.names = FALSE)
    attr(sigs_inc, "nexp") <- length(dec)
    attr(sigs_inc, "signames") <- names(inc)
    # attr(sigs_inc, "signames") <- sub("^(.*\\d+)_.+", "\\1", names(sigs_inc))

    list(dec = sigs_dec, inc = sigs_inc)
}
