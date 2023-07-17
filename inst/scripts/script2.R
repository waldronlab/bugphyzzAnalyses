library(dplyr)
library(tibble)

tbl1 <- tribble(
    ~NCBI_ID, ~ Attribute, ~ Score, ~ Evidence,
    '1', 'A', 1, 'exp',
    '1', 'B', 0, '',
    '1', 'C', 0, '',
    '1', 'D', 0, ''
)

tbl2 <- tribble(
    ~ NCBI_ID, ~ Attribute, ~ Score, ~ Evidence,
    '2', 'A', 0.5, 'igc',
    '2', 'B', 0.5, 'exp',
    '2', 'C', 0, '',
    '2', 'D', 0, ''
)

tbl3 <- tribble(
    ~ NCBI_ID, ~ Attribute, ~ Score, ~ Evidence,
    '3', 'A', 0.8, 'igc',
    '3', 'B', 0, '',
    '3', 'C', 0, '',
    '3', 'D', 0, ''
)

tbl4 <- tribble(
    ~ NCBI_ID, ~ Attribute, ~ Score, ~ Evidence,
    '4', 'A', 0, 'igc',
    '4', 'B', 0, 'exp',
    '4', 'C', 0, '',
    '4', 'D', 0, ''
)

x <- bind_rows(tbl1, tbl2, tbl3, tbl4)

x |>
    group_by(NCBI_ID) |>
    mutate(Score =  Score / sum(Score)) |>
    mutate(Score = ifelse(is.na(Score), 0, Score)) |>
    group_by(Attribute) |>
    reframe(
        Score = sum(Score),
        Evidence = paste0(Evidence, collapse = '|'),
    ) |>
    mutate(Score = Score / sum(Score)) |>
    mutate(Score = ifelse(is.na(Score), 0, Score)) |>
    mutate(Evidence = sub('^(\\|*)', '', Evidence)) |>
    mutate(Evidence = sub('\\|\\|+', '|', Evidence))
