library(taxPPro)
library(dplyr)
library(purrr)

## Get different sets of taxids
b <- get_ncbi_taxids(keyword = 'b')

bi <- get_ncbi_taxids(keyword = 'bi')
bl <- get_ncbi_taxids(keyword = 'bl')
bt <- get_ncbi_taxids(keyword = 'bt')

blt <- get_ncbi_taxids(keyword = 'blt')
bit <- get_ncbi_taxids(keyword = 'bit')
bli <- get_ncbi_taxids(keyword = 'bli')

blit <- get_ncbi_taxids(keyword = 'blit')


calc_numbers <- function(x) {
    ranks <- c('genus', 'species', 'strain')
    output <- purrr::map_int(ranks, ~ {
        df <- x[x$Rank == .x,]
        length(unique(df$NCBI_ID))
    })
    names(output) <- ranks
    output
}

data <- list(
    base = b,
    `base + informal names` = bi,
    `base + unclassified` = bl,
    `base + uncultured` = bt,
    `base + unclassified + uncultured` = blt,
    `base + informal names + uncultured` = bit,
    `base + informal names + unclassified` = bli,
    `base + informal names + unclassified + uncultured` = blit

)

summary <- bind_rows(map(data, calc_numbers), .id = 'ncbi_sets')
summary <- summary |>
    mutate(
        ncbi_sets = case_when(
            ncbi_sets == 'base' ~ paste0('[gn/sp] base'),
            ncbi_sets == 'base + informal names + unclassified + uncultured' ~ paste0('[st] base + informal name + unclassified + uncultured'),
            TRUE ~ ncbi_sets
        )
    )
summary






g <- b |>
    filter(
        NCBI_ID %in% unique(ncbi_taxonomy$NCBI_ID),
        kingdom == 'Bacteria', phylum == 'Thermotogae',
        Rank == 'genus'
    )

s <- b |>
    filter(
        NCBI_ID %in% unique(ncbi_taxonomy$NCBI_ID),
        kingdom == 'Bacteria', phylum == 'Thermotogae',
        Rank == 'species'
    )


genus <- propagation$aerophilicity |>
    filter(
        NCBI_ID %in% ncbi_taxonomy$NCBI_ID,
        Rank == 'genus'
    ) |>
    left_join(ncbi_taxonomy, by = 'NCBI_ID') |>
    # select(kingdom, phylum, NCBI_ID) |>
    distinct() |>
    filter(phylum == 'Thermotogae')

species <- propagation$aerophilicity |>
    filter(
        NCBI_ID %in% ncbi_taxonomy$NCBI_ID,
        Rank == 'species'
    ) |>
    left_join(ncbi_taxonomy, by = 'NCBI_ID') |>
    # select(kingdom, phylum, NCBI_ID) |>
    distinct() |>
    filter(phylum == 'Thermotogae')

View(species)









ncbi_summary_sp |>
    filter(phylum == 'Thaumarchaeota')

prop_summary_sp |>
    filter(phylum == 'Thaumarchaeota')



