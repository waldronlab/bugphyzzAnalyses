getCatSignatures <- function(bsdb, tax_level, min_size = 5) {
    y <- bsdb |>
        {\(y) split(y, y$`Abundance in Group 1`)}() |>
        purrr::map(~ {
            bugsigdbr::getSignatures(
                df = .x, tax.id.type = "ncbi", tax.level = tax_level,
                min.size = 1
            )
        })
    cond_names <- c(names(y$decreased), names(y$increased))
    cond_names <- sub("^bsdb:\\d+/\\d+/\\d+_(.+):.*$", "\\1", cond_names)
    cond_names <- sort(unique(cond_names))
    decreased <- vector("list", length(cond_names))
    for (i in seq_along(decreased)) {
        names(decreased)[i] <- cond_names[i]
        rgx <- paste0("bsdb:\\d+/\\d+/\\d+_", cond_names[i], ":")
        pos <- grep(rgx, names(y$decreased))
        if (!length(pos))
            next
        dec <- y$decreased[pos]
        dec <- unlist(dec, use.names = FALSE)
        decreased[[i]] <- dec
    }
    decreased <- discard(decreased, ~ length(.x) < 5)

    increased <- vector("list", length(cond_names))
    for (i in seq_along(increased)) {
        names(increased)[i] <- cond_names[i]
        rgx <- paste0("bsdb:\\d+/\\d+/\\d+_", cond_names[i], ":")
        pos <- grep(rgx, names(y$increased))
        if (!length(pos))
            next
        inc <- y$increased[pos]
        inc <- unlist(dec, use.names = FALSE)
        increased[[i]] <- inc
    }
    increased <- discard(increased, ~ length(.x) < min_size)
    common_cond_names <- sort(intersect(names(decreased), names(increased)))
    decreased <- decreased[common_cond_names]
    increased <- increased[common_cond_names]
    if (!length(increased))
        return(NULL)
    return(list(decreased = decreased, increased = increased))
}

