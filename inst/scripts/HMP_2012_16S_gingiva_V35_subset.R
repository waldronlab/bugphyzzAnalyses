
## This script is for getting the taxids at the genus level
## of the OTU's of the HMP_2012_16S_gingival_V35_subset dataset
## In MicrobiomeBenchmarkData

library(MicrobiomeBenchmarkData) # github: waldronlab/MicrobiomeBenchmarkData
library(taxize)
library(dplyr)
# library(bugphyzzAnalyses) #

dat_name <- 'HMP_2012_16S_gingival_V35_subset'

tse <- getBenchmarkData(dat_name, dryrun = FALSE)[[1]]
row_data <- tse |>
    rowData() |>
    as.data.frame() |>
    tibble::rownames_to_column('OTU') |>
    tibble::as_tibble()

new_row_data <- taxNames2TaxIDs(df = row_data, names_from = 'OTU')
if (all(rownames(tse) == new_row_data$otu)) {
    message('all good')
    rownames(new_row_data) <- new_row_data$otu
}
new_row_data <- select(new_row_data, -otu)
new_row_data <- relocate(new_row_data, taxon_annotation, .after = 'genus')

## Export files

output_file <- paste0(dat_name, '_taxids.tsv')
write.table(
    x = new_row_data, file = paste0('inst/extdata/', output_file),
    row.names = TRUE, sep = '\t'
)



