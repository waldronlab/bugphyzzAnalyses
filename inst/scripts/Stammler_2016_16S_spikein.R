
## I'll have to do this with some manual code
library(MicrobiomeBenchmarkData)
library(tidyr)

tse <- getBenchmarkData('Stammler_2016_16S_spikein', dryrun = FALSE)[[1]]
row_data <- tse |>
    rowData() |>
    as.data.frame()

taxa_list <- row_data |>
    t() |>
    as.data.frame() |>
    as.list()
names(taxa_list) <- rownames(tse)
taxa_list <- lapply(
    taxa_list, function(.x) {
        splitted <- stringr::str_split(.x, pattern = "; __")
        splitted[[1]]
    }
)
unique_taxa_list <- unique(taxa_list)



