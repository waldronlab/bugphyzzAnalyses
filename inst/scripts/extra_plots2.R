
library(ComplexHeatmap)


url <- 'https://raw.githubusercontent.com/waldronlab/BugSigDBPaper/main/inst/extdata/condition2category.txt'
cond <- read.csv(url, header = FALSE, col.names = c('Condition', 'cat'))
df <- bind_rows(enrichment_res, .id = 'comb') |>
    separate(col = 'comb', into = c('body_site', 'rank'), sep = '__')
df <- left_join(df, cond, by = 'Condition') |>
    mutate(Attribute_group = sub('^bp:', '', sub('\\|.*$', '', bp_sig)) )

## Frequencies per attribute group
myDF1 <- df |>
    count(Attribute_group) |>
    arrange(-n)
myPlot1 <- myDF1 |>
    mutate(Attribute_group = forcats::fct_inorder(Attribute_group)) |>
    slice_max(n, n = 10) |>
    ggplot(aes(x = reorder(Attribute_group, desc(n)), y = n)) +
    geom_col(fill = 'dodgerblue3') +
    labs(
        x = 'Attribute group', y = '# BSDB signatures'
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust  = 1)
    )

## Attribue x disease
myDat1 <- df |>
    filter(!is.na(cat)) |>
    filter(cat != 'environment, sociodemographics, lifestyle') |>
    count(Attribute_group, cat) |>
    rename(Disease_type = cat) |>
    pivot_wider(
        names_from = 'Disease_type', values_from = n, values_fill = 0
    ) |>
    tibble::column_to_rownames(var = 'Attribute_group') |>
    as.matrix()
myDat2 <- df |>
    filter(!is.na(cat)) |>
    count(Attribute_group, body_site) |>
    # rename(Disease_type = body_site) |>
    pivot_wider(
        names_from = 'body_site', values_from = n, values_fill = 0
    ) |>
    tibble::column_to_rownames(var = 'Attribute_group') |>
    as.matrix()

max_n <- max(max(myDat1), max(myDat2))
color_fun <- circlize::colorRamp2(
    breaks = c(0, max_n), colors = c('white', 'dodgerblue4')
)
lgd <- Legend(col_fun = color_fun, title = "# sigs")

ht1 <- Heatmap(
    matrix = myDat1,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_side = 'left',
    show_heatmap_legend = FALSE,
    col = color_fun,
    column_names_rot = 45,
    column_title = "Disease/condition type",
    row_names_max_width = max_text_width(
        rownames(myDat),
        gp = gpar(fontsize = 12)
    )
    # cell_fun = function(j, i, x, y, width, height, fill) {
    #       grid.rect(x = x, y = y, width = width, height = height,
    #           gp = gpar(col = "black", fill = NA))
    # }
)

ht2 <- Heatmap(
    matrix = myDat2,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_side = 'left',
    show_heatmap_legend = TRUE,
    name = '# sigs',
    col = color_fun,
    column_title = "Body site",
    column_names_rot = 45,
    row_names_max_width = max_text_width(
        rownames(myDat),
        gp = gpar(fontsize = 12)
    )
    # cell_fun = function(j, i, x, y, width, height, fill) {
    #       grid.rect(x = x, y = y, width = width, height = height,
    #           gp = gpar(col = "black", fill = NA))
    # }
)

draw(ht1 + ht2, ht_gap = unit(0.4, "cm"))


png(
    filename = 'vignettes/articles/attr_bs_dt.png', width = 8.4, height = 9.5,
    units = 'in', res = 300
)
draw(ht1 + ht2, ht_gap = unit(0.4, "cm"))
dev.off()


## Attribute x body site

