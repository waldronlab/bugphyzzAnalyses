library(bugphyzz)
library(taxPPro)
library(tidyr)
library(ggplot2)

phys <- physiologies(remove_false = TRUE)

## Only use physiologies with categorical values.
## For numeric and rages use thresholds instead of real values.
phys <- phys[vapply(phys, function(x) unique(x$Attribute_type) == 'categorical', logical(1))]
categorical <- lapply(
    phys, preSteps, tax.id.type = 'NCBI_ID', remove_false = TRUE
)

exclude_names <- c(
    'antimicrobial resistance',
    'isolation site'
)

categorical <- categorical[!names(categorical) %in% exclude_names]
propagation <- vector('list', length(categorical))
for (i in seq_along(propagation)) {
    current_name <- names(categorical)[i]
    message('>>> Propagating ',  names(categorical)[i], '. <<<')
    stat_time <- Sys.time()
    propagation[[i]] <- propagate(categorical[[i]])
    end_time <- Sys.time()
    message('>>> Done ', difftime(end_time, stat_time), '. <<<')
}
names(propagation) <- names(categorical)

df <- data.frame(
    before = vapply(categorical[names(propagation)], nrow, integer(1)),
    after = vapply(propagation, nrow, integer(1))
)
df$dataset <- rownames(df)

tidy_df <- df |>
    pivot_longer(
        names_to = 'type', values_to = 'n', cols = c('before', 'after')
    )
tidy_df$type <- factor(tidy_df$type, levels = c('before', 'after'))
tidy_df |>
    ggplot(aes(dataset, n)) +
    geom_col(aes(fill = type), position = 'dodge') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Get all numbers (including informal, unclassified, and uncultured)

## Let's get total by species with base numbers
summary_species <- base_numbers |>
    filter(Rank == 'species') |>
    count(kingdom, phylum, name = 'n_ncbi_sp') |>
    arrange(desc(kingdom), -n_ncbi_sp)

## Let's get totals by genus with base numbers
summary_genus <- base_numbers |>
    filter(Rank == 'genus') |>
    count(kingdom, phylum, name = 'n_ncbi_gn') |>
    arrange(desc(kingdom), -n_ncbi_gn)

## Let's get totals by strain with whole numbers?
summary_strain <- all_numbers |>
    filter(Rank == 'strain') |>
    count(kingdom, phylum, name = 'n_ncbi_st') |>
    arrange(desc(kingdom), -n_ncbi_st)











