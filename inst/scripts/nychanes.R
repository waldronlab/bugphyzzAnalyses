library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(tidySummarizedExperiment)
library(mia)
library(nychanesmicrobiome)
library(dplyr)
library(tibble)

lvls <- c('Never smoker', 'Cigarette')
nh <- loadQiimeData() |>
    annotateFactors() |>
    makeTreeSummarizedExperimentFromPhyloseq() |>
    filter(smokingstatus %in% c('Cigarette', 'Never smoker')) |>
    mutate(smokingstatus = factor(smokingstatus, levels = lvls))
rowData(nh)[c("New.CleanUp.ReferenceOTU8965", "New.ReferenceOTU179"), 'Genus'] <- 'uncultured'
rowData(nh)['New.ReferenceOTU72', 'Genus'] <- 'Family XIII'
rowData(nh)[is.na(rowData(nh)$Genus), 'Genus'] <- 'unclassified'

# taxids <- as.character(taxize::get_uid(rowData(nh)$Genus))
df <- data.frame(Genus = rowData(nh)$Genus, GenusID = taxids)
row_data <- rownames_to_column(as.data.frame(rowData(nh)), var = 'rownames')
row_data <- left_join(
    row_data, df, by = 'Genus', relationship = 'many-to-many'
) |>
    distinct()
df2 <- column_to_rownames(row_data, var = 'rownames')
df2 <- df2[rownames(nh),]
rowData(nh) <- S4Vectors::DataFrame(df2)

saveRDS(object = nh, file = 'inst/extdata/nychanes.RDS')
j <- readRDS('inst/extdata/nychanes.RDS')




