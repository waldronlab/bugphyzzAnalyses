## Run this after running the vignette with high-throughput analysis

bp_sigs |> length()
bsdb_sigs |> length()
tms_sigs |> names() |> head()

taxids <- unique(unlist(c(bp_sigs, bsdb_sigs, tms_sigs), recursive = TRUE, use.names = FALSE))

## Get names
taxnames <- taxizedb::taxid2name(taxids, db = "ncbi")
names(taxnames) <- taxids
missing_pos <- which(is.na(taxnames))
missing_names <- taxize::classification(names(taxnames)[missing_pos], db = "ncbi") |>
    map_chr(~ as.character(pull(tail(.x, 1), name)))
taxnames[missing_pos] <- missing_names

## Get ranks
taxranks <- taxizedb::taxid2rank(taxids, db = "ncbi")
names(taxranks) <- taxids
missing_pos <- which(is.na(taxranks))
missing_ranks <- taxize::classification(names(taxranks)[missing_pos], db = "ncbi") |>
    map_chr(~ as.character(pull(tail(.x, 1), rank)))
taxranks[missing_pos] <- missing_ranks

myDat <- data.frame(
    taxid = as.character(unname(taxids)),
    taxname = unname(taxnames),
    rank = unname(taxranks)

)

list2df <- function(l) {
    output <- vector("list", length(l))
    for (i in seq_along(output)) {
        output[[i]] <- data.frame(
            taxid = as.character(l[[i]]),
            sig_name = names(l)[i]
        )
    }
    bind_rows(output)
}

x <- list2df(bsdb_sigs) |>
    rename(bsdb_sig = sig_name)
y <- list2df(bp_sigs) |>
    rename(bp_sig = sig_name)
z <- list2df(tms_sigs) |>
    rename(tms_sig = sig_name)

sigTbl <- reduce(
    .x = list(myDat, x, y, z),
    .f = ~ left_join(.x, .y,  by = "taxid", relationship = "many-to-many")
) |>
    filter(!(is.na(bsdb_sig) & is.na(tms_sig)))

# fname <- "inst/extdata/sig_table_signatures.tsv"
# fname <- "inst/extdata/sig_table_metasignatures.tsv"
fname <- "inst/extdata/sig_table_metasignatures_nointersect.tsv"

write.table(
    sigTbl, quote = FALSE, row.names = FALSE, sep = "\t", file = fname
)

