
#' Create matrix for dbBact enrichment approach
#'
#' \code{dbMat} creates a count matrix for control (decreased) and
#' case (incraesed) sets and signature terms (e.g., from bugphyzz).
#'
#' @param control Character vector. Down-regulated set (control/decreased).
#' Number of experiments must have been added as the "nexp" attribute to the
#' signature.
#' @param case Character vector. Up-regulated set (case/increased).
#' Number of experiments must have been added as the "nexp" attribute to the
#' signature.
#' @param term_list A list of signature terms (e.g., from bugphyzz).
#' @param freq Integer. Minimum frequency. Default is 1.
#'
#' @return A SummarizedExperiment.
#' @export
#'
dbMat <- function(control, case, term_list, freq = 1) {

    mat_list <- vector("list", length(term_list))

    for (i in seq_along(mat_list)) {
        term <- term_list[[i]]

        co <- as.integer(c(control %in% term))
        names(co) <- control
        co <- tapply(co, names(co), sum)
        co <- co[sort(names(co))]

        ca <- as.integer(c(case %in% term))
        names(ca) <- case
        ca <- tapply(ca, names(ca), sum)
        ca <- ca[sort(names(ca))]

        mat <- matrix(data = c(co, ca), nrow = 1)
        colnames(mat) <- c(paste0("co_", names(co)), paste0("ca_", names(ca)))
        rownames(mat) <- names(term_list)[i]
        mat_list[[i]] <- mat
    }

    mat <- do.call("rbind", mat_list)
    mat <-  mat[rowSums(mat) > 0, , drop = FALSE]
    if (nrow(mat) == 0) {
        return(NULL)
    }
    rownames(mat) <- sub("^bugphyzz:", "", rownames(mat))

    ## nexp is assigned by the concatenateEvenSigs
    if (attr(control, "nexp") != attr(case, "nexp")) {
        warning(
            "Number of experiments don't match in control and case.",
            call. = FALSE
        )
    }

    n_exp_ctrl <- attr(control, "nexp")
    n_exp_case <- attr(case, "nexp")

    control_freq <- table(control)
    case_freq <- table(case)

    cond <- c(rep("Control", length(co)), rep("Case", length(ca)))
    cond <- factor(cond, levels = c("Control", "Case"))

    Taxon <- c(names(co), names(ca))
    Taxon2 <- taxizedb::taxid2name(Taxon, db = "ncbi", verbose = FALSE)
    names(Taxon2) <- Taxon

    for (i in seq_along(Taxon2)) {
        if (is.na(Taxon2[i])) {
            # classif <- taxize::classification(names(Taxon2)[i], db = "ncbi")[[1]]
            # Taxon2[i] <- tail(classif$name, 1)
            Taxon2[i] <- names(Taxon2)[i]
        }
    }

    Frequency = c(control_freq[names(co)], case_freq[names(ca)])

    dat <- data.frame(
        Condition = cond,
        Taxon = Taxon,
        Taxon2 = Taxon2,
        Frequency = Frequency
    )

    rownames(dat) <- c(paste0("co_", names(co)), paste0("ca_", names(ca)))

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(Scores = mat),
        colData = S4Vectors::DataFrame(dat),
        metadata = list(
            nexp_ctrl = n_exp_ctrl,
            nexp_case = n_exp_case,
            control = control,
            case = case
        )
    )

    pass_control <- sum(control_freq >= freq) >= 5
    pass_case <- sum(case_freq >= freq) >= 5

    if (isFALSE(pass_control && pass_case)) {
        return(NULL)
    }

    se <- se[, which(Frequency >= freq)]
    return(se)
}

calcEffectSize <- function(se, perm = 1000) {
    mat <- SummarizedExperiment::assay(se, "Scores")
    col_data <- SummarizedExperiment::colData(se)
    Control <- col_data$Taxon[col_data$Condition == "Control"]
    Case <- col_data$Taxon[col_data$Condition == "Case"]
    es <- apply(mat, 1, function(x) {
        co <- x[1:length(Control)]
        ca <- x[(length(Control)+1):ncol(mat)] # increased
        # round((sum(ca) / length(Case)) - (sum(co) / length(Control)), 3)
        mean(ca) - mean(co)
    })

    pVals <- apply(mat, 1, function(x) {
        co <- x[1:length(Control)]
        ca <- x[(length(Control)+1):ncol(mat)] # increased
        observed_es <- mean(ca) - mean(co)
        combined <- c(co, ca)
        permuted_es <- vector("numeric", perm)
        for (i in seq_along(permuted_es)) {
            shuffled_data <- sample(combined)
            shuffled_group1 <- shuffled_data[1:length(co)]
            shuffled_group2 <- shuffled_data[(length(co) + 1):length(shuffled_data)]
            permuted_es[i] <- mean(shuffled_group2) - mean(shuffled_group1)
        }
        mean(abs(permuted_es) >= abs(observed_es))
    })
    dir <- dplyr::case_when(
        es < 0 ~ "Control",
        es > 0 ~ "Case",
        es == 0 ~ "Unchanged",
        TRUE ~ NA
    )
    dir <- factor(dir, levels = c("Case", "Unchanged", "Control"))
    dat <- data.frame(
        Effect_size = es,
        Direction = dir,
        PermP = pVals
    )
    rownames(dat) <- names(es)
    SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(dat)
    return(se)
}

calcPvalue <- function(se) {
    mat <- SummarizedExperiment::assay(se, "Scores")
    col_data <- SummarizedExperiment::colData(se)
    Control <- col_data$Taxon[col_data$Condition == "Control"]
    Case <- col_data$Taxon[col_data$Condition == "Case"]
    p_val <- apply(mat, 1, function(x) {
        co <- x[1:length(Control)]
        ca <- x[(length(Control)+1):ncol(mat)] # increased
        res <- wilcox.test(co, ca, exact = FALSE)
        # res <- do.call(what = f, args = list(ca, co))
        round(res$p.value, 3)
    })
    if (all(rownames(se) == names(p_val))) {
        SummarizedExperiment::rowData(se)$P_value <- p_val
    }
    return(se)
}

#' Run enrichment using the dbBact approach
#'
#' @param control Character vector. Down-regulated set (control/decreased).
#' Number of experiments must have been added as the "nexp" attribute to the
#' signature.
#' @param case Character vector. Up-regulated set (case/increased).
#' Number of experiments must have been added as the "nexp" attribute to the
#' signature.
#' @param term_list A list of signature terms (e.g., from bugphyzz).
#' @param freq Integer. Minimum frequency. Default is 1.
#' @param perm Number of permutations for calculating P-value with mean
#' difference.
#'
#' @export
#'
dbEn2 <- function(control, case, term_list, freq = 1, perm = 1000) {
    se <- dbMat(
        control = control, case = case, term_list = term_list, freq = freq
    )
    if (is.null(se)) {
        return(NULL)
    }
    if (!ncol(se)) {
        warning(
            "Not enough prevalence", call. = FALSE
        )
        return(NULL)
    }
    se <- calcEffectSize(se, perm = perm)
    se <- calcPvalue(se) # Wilcox.test
    ef <- SummarizedExperiment::rowData(se)$Effect_size
    names(ef) <- rownames(se)
    se <- se[names(sort(ef, decreasing = TRUE)),]
    return(se)
}

#' Plot heatmap of output of dbEn2
#'
#' @param se A SummarizedExperiment
#' @param col_pad Padding at the bottom of the heatmap. Default is 2.
#' @param pCol Column to be used for getting the P-values. "P_value" or
#' "PermP". Default is P_value.
#'
#' @return A Heatmap.
#' @export
#'
dbHt <- function(se, col_pad = 2, pCol = "P_value") {

    if (is.null(se)) {
        warning("Input se was NULL. Returning NULL.", call. = FALSE)
        return(NULL)
    }

    se <- se[which(rowData(se)$Effect_size != 0),]

    mat <- SummarizedExperiment::assay(se, "Scores")
    n_exp_ctrl <- S4Vectors::metadata(se)$nexp_ctrl
    n_exp_case <- S4Vectors::metadata(se)$nexp_case

    score_name <- paste0("Score (nexp:", n_exp_ctrl, "|", n_exp_case, ")")

    col_data <- SummarizedExperiment::colData(se)
    row_data <- SummarizedExperiment::rowData(se)

    colnames(mat) <- col_data$Taxon2

    ## Color scale
    htColor <- function(mat) {
        circlize::colorRamp2(
            breaks = c(0, max(mat, na.rm = TRUE)),
            colors = c("white", "gray10")
        )
    }

    ## Top annotation
    top_ha <-  ComplexHeatmap::HeatmapAnnotation(
        "Frequency" = ComplexHeatmap::anno_barplot(
            col_data$Frequency,
            height = unit(5, "cm"),
            bar_width = 1,
            gp = grid::gpar(fill = "gold")
        ),
        foo = ComplexHeatmap::anno_block(
            gp = grid::gpar(fill = c("dodgerblue", "firebrick")),
            height = unit(1, "cm"),
            labels = levels(col_data$Condition),
            labels_gp = grid::gpar(col = "white", fontsize = 20, fontface = "bold")
        )
    )

    ## Left annotation
    vct_chr <- levels(droplevels(row_data$Direction))
    if (length(vct_chr) == 2) {
        vct_int <- c("firebrick", "dodgerblue")
    } else if (length(vct_chr) == 3) {
        vct_int <- c("firebrick", "white", "dodgerblue")
        vct_chr <- c("Case", "", "Control")
    }
    left_ha <- ComplexHeatmap::HeatmapAnnotation(
        foo2 = ComplexHeatmap::anno_block(
            gp = grid::gpar(fill = vct_int),
            width = unit(1, "cm"),
            labels = vct_chr,
            labels_gp = grid::gpar(col = "white", fontsize = 20, fontface = "bold")
        ),
        which = "row"
    )

    ## Right annotation
    pValCol <- function(values) {
        circlize::colorRamp2(
            breaks = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)),
            colors = c("white", "darkcyan")
        )
    }
    pValColWhite <- function(values) {
        circlize::colorRamp2(
            breaks = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)),
            colors = c("white", "white")
        )
    }
    log10_pval <- -log10(row_data[[pCol]] + 1)
    pch_var <- case_when(
        row_data[[pCol]] < 0.1 ~ 8,
        TRUE ~ NA
    )
    right_ha <- ComplexHeatmap::HeatmapAnnotation(
        "-log10(pval+1)" = ComplexHeatmap::anno_simple(
            x = -log10(row_data[[pCol]] + 1),
            col = pValCol(c(-0.30103, 0)),
            border = TRUE
        ),
        "pval < 0.1" = ComplexHeatmap::anno_simple(
            x = -log10(row_data[[pCol]] + 1),
            col = pValColWhite(-log10(row_data[[pCol]] + 1)),
            pch = pch_var
        ),

        "Mean difference" = ComplexHeatmap::anno_barplot(
            row_data$Effect_size, width = unit(5, "cm"),
            bar_width = 1,
            gp = grid::gpar(fill = "gold")
        ),
        which = "row"
    )

    ## Legend for draw
    pval_lgd <- ComplexHeatmap::Legend(
        col_fun = pValCol(c(-0.30103, 0)),
        title = "-log10(pval+1)", direction = "horizontal",
        title_position = "topcenter"
    )

    ## Create Heatmap
    ht <- ComplexHeatmap::Heatmap(
        matrix = mat,

        ## Heatmap color and legend
        name = score_name,
        col = htColor(mat),
        heatmap_legend_param = list(
            direction = "horizontal", title_position = "topcenter"
        ),

        ## Clustering
        cluster_columns = FALSE,
        cluster_rows = FALSE,

        ## Splitting
        column_split = col_data$Condition,
        column_title = NULL,
        column_gap = unit(3, "mm"),

        row_split = row_data$Direction,
        row_title = NULL,
        row_gap = unit(3, "mm"),

        ## Row names
        row_names_side = "right",
        row_names_max_width = ComplexHeatmap::max_text_width(
            rownames(mat),
            gp = grid::gpar(fontsize = 12)
        ),

        ## Column names
        column_names_gp = grid::gpar(fontface = "italic"),

        ## Annotations
        top_annotation = top_ha,
        left_annotation = left_ha,
        right_annotation = right_ha,

        ## Borders
        border = TRUE,
        rect_gp = grid::gpar(col = "gray50", lwd = 0.25)

    )
    ComplexHeatmap::draw(
        ht, heatmap_legend_side = "bottom",
        annotation_legend_list = list(pval_lgd),
        annotation_legend_side = "bottom",
        merge_legends = TRUE,
        padding = unit(c(col_pad, 2, 2, 2), "mm")
    )
}

## Functions
## Function for calculating scores
## Scores need to be added to the signatures (or not)
dbEn <- function(s1, s2, t) {
    size1 <- length(s1)
    size2 <- length(s2)
    scores1 <- attr(t, "Scores")[s1] |>
        {\(y) y[!is.na(y)]}()
    scores2 <- attr(t, "Scores")[s2] |>
        {\(y) y[!is.na(y)]}()
    (sum(scores1) / size1) - (sum(scores2) / size2)
}

## Function for
dbEnBoot <- function(s1, s2, t, R = 1000) {
    observed_score <- dbEn(s1, s2, t)
    perm_scores <- vector("numeric", R)
    for (i in seq_along(perm_scores)) {
        cat_sets <- c(s1, s2)
        permuted_sets <- sample(cat_sets)
        perm_s1 <- permuted_sets[1:length(s1)]
        perm_s2 <- permuted_sets[(length(s1) + 1):length(permuted_sets)]
        perm_scores[i] <- dbEn(perm_s1, perm_s2, t)
    }
    p_value <- mean(abs(perm_scores) >= abs(observed_score))
    list(score = observed_score, p_value = p_value)
}
