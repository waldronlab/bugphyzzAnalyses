
library(ComplexHeatmap)
library(bugphyzzAnalyses)
library(stringr)
library(dplyr)
library(bugsigdbr)

bsdb <- importBugSigDB()
bsdb$Condition <- tolower(bsdb$Condition)

url <- 'https://raw.githubusercontent.com/waldronlab/BugSigDBPaper/main/inst/extdata/condition2category.txt'
cond_cat_paper <- read.csv(url, header = FALSE, col.names = c('Condition', 'Category'))
cond_cat_paper$Condition <- tolower(str_squish(cond_cat_paper$Condition))

url2 <- system.file('extdata', 'condition2category.tsv',  package = 'bugphyzzAnalyses', mustWork = TRUE)
cond_cat_ba <- read.table(url2, sep = '\t', header = TRUE)
cond_cat_ba$Condition <- tolower(str_squish(cond_cat_ba$Condition))

conditions_bsdb <- unique(tolower(str_squish(bsdb$Condition)))


mean(conditions_bsdb %in% unique(cond_cat_paper$Condition))
conditions_bsdb[!conditions_bsdb %in% unique(cond_cat_ba$Condition)]




paper_conditions <- unique(cond_cat_paper$Condition)
bsdb_current_new_conditions <- bsdb |>
    mutate(Condition = tolower(str_squish(Condition))) |>
    filter(!is.na(Condition)) |>
    filter(!Condition %in% paper_conditions) |>
    pull(Condition) |>
    unique()









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
)

draw(ht1 + ht2, ht_gap = unit(0.4, "cm"))


png(
    filename = 'vignettes/articles/attr_bs_dt.png', width = 8.4, height = 9.5,
    units = 'in', res = 300
)
draw(ht1 + ht2, ht_gap = unit(0.4, "cm"))
dev.off()


## Attribute x body site

x <- summary1 |>
    select(-total_n) |>
    pivot_wider(
        names_from = 'rank', values_from = 'n', values_fill = 0
    ) |>
    tibble::column_to_rownames(var = 'Attribute_group') |>
    as.matrix()

anno_width <- unit(3, 'cm')
rank_col <- c(species = 'blue', genus = 'cyan')
row_anno <- rowAnnotation("Frequency" = anno_barplot(x[attr_grp_nms,], bar_width = 1, gp = gpar(fill = rank_col),
                                        width = anno_width), show_annotation_name = FALSE)
ht1 + row_anno




cond_names <- colnames(summary_condition_mat)



col_anno <- columnAnnotation(
    "x" = anno_barplot(bsdb_summary_cat_matrix[cond_names,], bar_height = 1, gp = gpar(fill = 'green'),
                                        height = anno_width), show_annotation_name = FALSE)

columnAnnotation(
    conditions = anno_barplot()
)


ht1 + row_anno



ht1 + col_anno

HeatmapAnnotation(df = bsdb_summary_cat_matrix[,cond_names,]




col_anno <- HeatmapAnnotation(
    df = bsdb_summary_cat_matrix[cond_names,], which = 'column'
)

ht1 + col_anno




