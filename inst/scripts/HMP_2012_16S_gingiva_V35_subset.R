
## This script is for getting the taxids at the genus level
## of the OTU's of the HMP_2012_16S_gingival_V35_subset dataset
## In MicrobiomeBenchmarkData

library(MicrobiomeBenchmarkData) # github: waldronlab/MicrobiomeBenchmarkData
library(taxize)
library(dplyr)
library(bugphyzzAnalyses) # github: waldronlab/bugphyzzAnalyses

dat_name <- 'HMP_2012_16S_gingival_V35_subset'

tse <- getBenchmarkData(dat_name, dryrun = FALSE)[[1]]
row_data <- tse |>
    rowData() |>
    as.data.frame() |>
    tibble::rownames_to_column('OTU') |>
    tibble::as_tibble()

# taxids1 <- df |>
    # filter(is.na(genus))
# taxids1 <- taxids1[1,]
# taxids1 <- as.character(taxids1)
# taxids1 <- taxids1[length(taxids1):1]


# x <- my_fun(taxids1)

# output <- vector('list', length(taxids1))
df <- taxNames2TaxIDs(x = row_data, names_from = 'OTU')

test_out <- lapply(df, .getTaxonClassification)
lapply(test_out, classif2Table) |>
    bind_rows()


tax_names <- unique(row_data$genus)
tax_names <- tax_names[!is.na(tax_names)]

classif <- classification(tax_names, db = 'ncbi')
classif <- classif[!is.na(classif)]

taxids_df <- classif |>
    lapply(classif2Table) |>
    dplyr::bind_rows(.id = 'genus_name') |>
    dplyr::rename(kingdom = superkingdom) |>
    {\(y) magrittr::set_colnames(y, paste0(colnames(y), '_taxid'))}()


final <- left_join(row_data, taxids_df, by = c('genus' = 'genus_name_taxid'))
final <- final |>
    dplyr::select(OTU, tidyselect::ends_with("_taxid"), taxon_annotation)
colnames(final) <- sub(".y$", "", colnames(final))

new_row_data <- final |>
    relocate(genus, .after = 'family') |>
    relocate(taxon_annotation, .after = 'genus') |>
    tibble::column_to_rownames('OTU')


## Export files

output_file <- paste0(dat_name, '_taxids.tsv')
tse2 <- tse
rowData(tse2) <- S4Vectors::DataFrame(new_row_data)

readr::write_tsv(x= genera_taxids, file = paste0('inst/extdata/', output_file))

write.table(
    x = new_row_data, file = paste0('inst/extdata/', output_file),
    row.names = TRUE, sep = '\t'
)



