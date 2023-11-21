lenGn <- bp |>
    filter(Attribute_group == 'length') |>
    filter(Rank == 'genus')
gn <- node_data |>
    filter(Rank == 'genus') |>
    pull(taxid) |>
    unique()
lenGn1 <- lenGn |>
    filter(Evidence %in% c('exp', 'nas', 'tas', 'igc')) |>
    pull(NCBI_ID)
mean(lenGn1 %in% gn)


lenSp <- bp |>
    filter(Attribute_group == 'length') |>
    filter(Rank == 'species')


lenSp1 <- lenSp |>
    filter(Evidence == 'inh')


len_sp <- lenSp1 |>
    pull(NCBI_ID) |>
    unique() |>
    taxizedb::classification(db = 'ncbi')

a <- map(len_sp, ~ {
    filter(.x, rank == 'genus') |>
        pull(id)
}) |>
    unlist() |>
    unique()

mean(a %in% gn)

mean(gn %in% a)
