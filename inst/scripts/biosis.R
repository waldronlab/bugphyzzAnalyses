library(dplyr)
fname <- system.file('extdata/biosis.tsv', package = 'nychanesmicrobiome')
col_names <- c('Genus', 'Attribute')
biosis <- read.table(fname, header = TRUE, sep = '\t', col.names = col_names)
biosis <- biosis |>
    mutate(Attribute = tolower(Attribute)) |>
    mutate(
        Attribute = ifelse(
            Attribute == 'f anaerobic', 'facultatively anaerobic', Attribute
        )
    ) |>
    mutate(Attribute = gsub(' ', '_', Attribute))
## This needs to be done in an interactive session because some queries
## need manual confirmation
taxids <- taxize::get_uid(biosis$Genus, db = 'ncbi')
biosis$GenusID <- as.character(taxids)

write.table(
    biosis, file = 'inst/extdata/biosis.tsv', row.names = FALSE,
    col.names = TRUE, sep = '\t'
)
