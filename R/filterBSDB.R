#' Concatenate signatures of a BugSigDB data.frame
#'
#' @param bsdb A BugSigDB data.frame.
#' @param tax_level Character string. "genus" or "species".
#' @param min_size Integer.
#'
#' @return Named nested list of signatures.
#' @export
#'
getCatSignatures <- function(bsdb, tax_level, min_size = 5) {

    sigs <- bsdb |>
        {\(y) split(y, y$`Abundance in Group 1`)}() |>
        purrr::map(~ {
            bugsigdbr::getSignatures(
                df = .x, tax.id.type = "ncbi", tax.level = tax_level,
                min.size = 1
            )
        })

    cond_names <- c(names(sigs$decreased), names(sigs$increased))
    cond_names <- sub("^bsdb:\\d+/\\d+/\\d+_(.+):.*$", "\\1", cond_names)
    cond_names <- sort(unique(cond_names))

    decreased <- vector("list", length(cond_names))
    increased <- vector("list", length(cond_names))

    for (i in seq_along(cond_names)) {
        rgx <- paste0("bsdb:\\d+/\\d+/\\d+_", cond_names[i], ":")

        names(decreased)[i] <- cond_names[i]
        pos <- grep(rgx, names(sigs$decreased))
        if (!length(pos))
            next
        dec <- sigs$decreased[pos]
        dec <- unlist(dec, use.names = FALSE)
        decreased[[i]] <- dec
        attr(decreased[[i]], "nexp") <- length(pos)

        names(increased)[i] <- cond_names[i]
        pos <- grep(rgx, names(sigs$increased))
        if (!length(pos))
            next
        inc <- sigs$increased[pos]
        inc <- unlist(inc, use.names = FALSE)
        increased[[i]] <- inc
        attr(increased[[i]], "nexp") <- length(pos)
    }

    decreased <- discard(decreased, ~ length(.x) < min_size)
    increased <- discard(increased, ~ length(.x) < min_size)

    common_cond_names <- sort(intersect(names(decreased), names(increased)))

    decreased <- decreased[common_cond_names]
    increased <- increased[common_cond_names]

    if (!length(increased))
        return(NULL)
    return(list(decreased = decreased, increased = increased))

}

#' Filter even dataset from BugSigDB
#'
#' @param bsdb  A BSDB data.frame.
#' @param rank Character string "genus" or "species".
#' @param min_size Number of signatures.
#'
#' @return A data.frame.
#' @export
#'
filterEvenDat <- function(bsdb, rank, min_size = 5) {
    bsdb_ids <- getSignatures(
        df = bsdb, tax.id.type = "ncbi", tax.level = rank, min.size = min_size
    )
    if (!length(bsdb_ids)) {
        return(NULL)
    }
    bsdb_ids <- bsdb_ids |>
        names() |>
        sub("^(bsdb:\\d+/\\d+/\\d+)_.*$", "\\1", x = _)
    dats <- bsdb |>
        filter(`BSDB ID` %in% bsdb_ids) |>
        group_by(Study, Experiment) |>
        mutate(count = n()) |>
        ungroup() |>
        filter(count == 2) |>
        group_by(Study, Experiment) |>
        arrange(`Abundance in Group 1`) |>
        mutate(comb = paste0(sort(`Abundance in Group 1`), collapse = "-")) |>
        ungroup() |>
        filter(comb == "decreased-increased") |>
        arrange(`BSDB ID`) |>
        {\(y) split(y, y$`Abundance in Group 1`)}()
    if (!length(dats)) {
        return(NULL)
    }
    return(dats)
}


#' Filter even signatures from experiments
#'
#' \code{filterEvenSignatures} selects only signatures from experiments
#' in BSDB that have exactly two signatures: increased and decreased.
#' From case-control studies only.
#'
#' @param bsdb A BugSigDB dataset from importBugSigDB.
#' @param rank Character string. "genus" or "species"
#' @param min_size Minimum number of the signatures to be included.
#'
#' @return A list of signatures nested by direction and condition.
#' @export
#'
filterEvenSignatures <- function(bsdb, rank, min_size = 10) {

    dats <- filterEvenDat(bsdb, rank, min_size)
    # bsdb_ids <- getSignatures(
    #     df = bsdb, tax.id.type = "ncbi", tax.level = rank, min.size = min_size
    # )
    # if (!length(bsdb_ids)) {
    #     return(NULL)
    # }
    # bsdb_ids <- bsdb_ids |>
    #     names() |>
    #     sub("^(bsdb:\\d+/\\d+/\\d+)_.*$", "\\1", x = _)
    # dats <- bsdb |>
    #     filter(`BSDB ID` %in% bsdb_ids) |>
    #     group_by(Study, Experiment) |>
    #     mutate(count = n()) |>
    #     ungroup() |>
    #     filter(count == 2) |>
    #     group_by(Study, Experiment) |>
    #     arrange(`Abundance in Group 1`) |>
    #     mutate(comb = paste0(sort(`Abundance in Group 1`), collapse = "-")) |>
    #     ungroup() |>
    #     filter(comb == "decreased-increased") |>
    #     arrange(`BSDB ID`) |>
    #     {\(y) split(y, y$`Abundance in Group 1`)}()
    # if (!length(dats)) {
    #     return(NULL)
    # }
    sigs <- dats |>
        map(~ {
            getSignatures(
                df = .x, tax.id.type = "ncbi",
                tax.level = rank, min.size = min_size
            )
        })
    names(sigs$decreased) <- sub(
        # "^(bsdb:\\d+/\\d+)/\\d+_.*$", "\\1", names(sigs$decreased)
        "^(bsdb:\\d+/\\d+)/\\d+_(.+):.*$", "\\1_\\2", names(sigs$decreased)
    )
    names(sigs$increased) <- sub(
        # "^(bsdb:\\d+/\\d+)/\\d+_.*$", "\\1", names(sigs$increased)
        "^(bsdb:\\d+/\\d+)/\\d+_(.+):.*$", "\\1_\\2", names(sigs$increased)
    )
    return(sigs)
}

#' Concatenate even signatures
#'
#' \code{concatenateEvenSigs} concatenate signatures that are the output
#' of the \code{filterEvenSignatures}.
#'
#' @param l A list of signatures.
#'
#' @return A nested list of signatures by direction and condition (concatenated).
#' @export
#'
concatenateEvenSigs <- function(l, min_size = 5) {
    cond_names <- unique(
        sub("^bsdb:\\d+/\\d_", "", names(l$decreased)),
        sub("^bsdb:\\d+/\\d_", "", names(l$increased))
    ) |>
        unique() |>
        sort()
    decreased <- vector("list", length(cond_names))
    increased <- vector("list", length(cond_names))
    for (i in seq_along(cond_names)) {
        rgx <- paste0("/\\d+_", cond_names[i], "$")
        select_var <- grep(rgx, names(l$decreased), value = TRUE)

        names(decreased)[i] <- cond_names[i]
        n_exp_dec <- length(l$decreased[select_var])
        sig_dec <- unlist(l$decreased[select_var], use.names = FALSE)

        names(increased)[i] <- cond_names[i]
        n_exp_inc <- length(l$increased[select_var])
        sig_inc <- unlist(l$increased[select_var], use.names = FALSE)

        vct_lgl <- length(sig_dec) >= min_size && length(sig_inc) >= min_size

        if (isFALSE(vct_lgl)) {
            next
        }

        decreased[[i]] <- sig_dec
        increased[[i]] <- sig_inc
        attr(decreased[[i]], "nexp") <- n_exp_dec
        attr(increased[[i]], "nexp") <- n_exp_inc
    }

    decreased <- purrr::discard(decreased, is.null)
    increased <- purrr::discard(increased, is.null)

    cat_sigs <- list(decreased = decreased, increased = increased)
    return(cat_sigs)
}
