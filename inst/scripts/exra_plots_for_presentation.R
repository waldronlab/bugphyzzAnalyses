## Temporary file for exporting plots for my presentation
## Run this code after unning the vignette for stats (and possiblty others)




png(
    filename = 'vignettes/articles/completeness_heatmap.png', width = 8.3, height = 8,
    units = 'in', res = 150
)
draw(ht_b + ht_a, ht_gap = unit(0.4, "cm"))
dev.off()


summary_source_p <- summary_source |>
    mutate(origin = ifelse(origin == 'Original source', 'Source', origin)) |>
    mutate(origin = fct_reorder(origin, desc(n))) |>
    ggplot(aes(reorder(Attribute_source, -n), n)) +
    geom_col(aes(fill = origin)) +
    labs(
        x = 'Origin', y = 'Number of annotations'
    ) +
    scale_y_continuous(limits = c(0, 1500000), labels = scales::comma) +
    scale_y_cut(breaks = c(100, 20000)) +
    scale_fill_discrete(name = '') +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
    )
ggsave(
    filename = 'vignettes/articles/summary_source_p.png',
    plot = summary_source_p, width = 7.8, height = 6.25
)


atr_gp_ev_p

ggsave(
    filename = 'vignettes/articles/total_counts.png', plot = attr_gp_rk_p,
    width = 11.3, height = 6.2
)




## Extra code

bp |>
    filter(Evidence %in% c('asr', 'inh')) |>
    pull(NCBI_ID) |>
    unique() |>
    length()

bp |>
    filter(!Evidence %in% c('asr', 'inh')) |>
    pull(NCBI_ID) |>
    unique() |>
    length()

selectCols <- c(
    "Bergey's Manual", "BacDive", "Browne_2021"
)
selectRows <- c(
    'aerophilicity', 'growth temperature', 'shape', 'spore shape'
)


data <- original |>
    filter(Attribute_source %in% selectCols) |>
    filter(Attribute_group %in% selectRows) |>
    count(Attribute_group, Attribute_source, Rank) |>
    mutate(new_column = paste0(Attribute_group, '|', Rank)) |>
    select(new_column, Attribute_source, n) |>
    pivot_wider(
        names_from = 'Attribute_source', values_from = 'n', values_fill = 0
    ) |>
    arrange(new_column) |>
    tibble::column_to_rownames(var = 'new_column') |>
    as.matrix()


bar_vct <- sub('^.*\\|', '', rownames(data))
split_vct <- sub('\\|.*$', '', rownames(data))
bar_color <- viridis::viridis(n = 3, option = 'C')
names(bar_color) <- c('genus', 'species', 'strain')
row_ha <- rowAnnotation(
    'Rank' = bar_vct, col = list(Rank = bar_color),
    show_annotation_name = FALSE
)
# log2_data <- log2(data + 1)
data[data > 0] <- 1

color_fun2 <- circlize::colorRamp2(
    breaks = c(0, 1), colors = c('white', 'black')
)
ht_sources <- Heatmap(
    matrix = data,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE,
    row_names_side = 'left',
    name = 'Presence/ausence',
    show_heatmap_legend = TRUE,
    col = color_fun2(c(0, 1)),
    left_annotation = row_ha,
    row_split = split_vct,
    column_names_rot = 45,
    row_title_rot = 0, gap = unit(0.1, 'inches'),
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height,
                  gp = gpar(col = "black", fill = NA))
    }
)

# png(filename = 'sup3.png', width = 10, height = 11, units = 'in', res = 300)
# draw(ht_sources)
# dev.off()

ht_sources




# png(filename = 'sup3.png', width = 10, height = 11, units = 'in', res = 300)
# draw(ht_sources)
# dev.off()

ht_sources






# png(filename = 'sup3.png', width = 10, height = 11, units = 'in', res = 300)
# draw(ht_sources)
# dev.off()

ht_sources
