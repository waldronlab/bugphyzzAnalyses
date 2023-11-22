library(bugphyzz)
library(taxizedb)
library(purrr)
library(dplyr)
bp <- importBugphyzz(force_download = TRUE)
taxids <- unique(as.character(bp$NCBI_ID))
system.time({
   taxonomy <- classification(taxids, db = 'ncbi')
})
valRanks <- c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')
system.time({
    x <- map(taxonomy, ~ {
        dat <- .x |>
            filter(rank %in% valRanks)
        ids <- dat |>
            pull(id)
        ranks <- dat |>
            pull(rank)
        mat <- matrix(ids, nrow = 1)
        colnames(mat) <- ranks
        as.data.frame(mat)
    })

})
taxonomy_ids <- bind_rows(x, .id = 'NCBI_ID')
write.table(
    x = taxonomy_ids, file = 'inst/extdata/taxonomy_ids.tsv', sep = '\t'
)
system.time({
    y <- map(taxonomy, ~ {
        dat <- .x |>
            filter(rank %in% valRanks)
        names <- dat |>
            pull(name)
        ranks <- dat |>
            pull(rank)
        mat <- matrix(names, nrow = 1)
        colnames(mat) <- ranks
        as.data.frame(mat)
    })

})
taxonomy_names <- bind_rows(y, .id = 'NCBI_ID')
write.table(
    x = taxonomy_names, file = 'inst/extdata/taxonomy_names.tsv', sep = '\t'
)
