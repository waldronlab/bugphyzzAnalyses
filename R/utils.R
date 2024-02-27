#' Filter taxa based on abundance values per sample
#'
#' \code{filterTaxa} filters the number of taxa per sample based on a minimum
#' value of abundance. This functions works with both phyloseq and
#' TreeSummarizedExperiment objects.
#'
#' @param x A phyloseq or TreeSummarizedExperiment object with otu_table/assay
#' and sample_data/colData
#' @param min_ab The minimum value of abundance for taxon to be
#' considered as present in a sample. Default is 1. The default value of 1
#' could be good for counts. Relative abundance or other data
#' transformations might require another threshold value.
#' @param min_per minimum percentage of samples in which each taxon must be
#' present in order to be kept in the data. Default is 0.2. Taxon presence is
#' dtermined by the `min_ab` argument (see above).
#'
#' @return The filtered phyloseq/TreeSummarized object
#' @export
#'
filterTaxa <- function(x, min_ab = 1, min_per = 0.2) {

    if (is(x, 'TreeSummarizedExperiment') || is(x, 'SummarizedExperiment')) {
        m <- SummarizedExperiment::assay(x)
        n_samples <- ncol(m)
        min_n_samples <- round(n_samples * min_per)
        return(x[rowSums(m >= min_ab) >= min_n_samples,])

    } else if (is(x, 'phyloseq')) {
        m <- phyloseq::otu_table(x)
        n_samples <- ncol(m)
        min_n_samples <- round(n_samples * min_per)
        return(phyloseq::prune_taxa(rowSums(m >= min_ab) >= min_n_samples, x))
        # m <- phyloseq::otu_table(ps)
        # return(phyloseq::prune_samples(colSums(m > 1) >= 2, ps)) ## apply higher filtering
    }
}

#' Import taxids
#'
#' @param x Name of the dataset from Microbiome BenchmarkData
#'
#' @return A data frame.
#' @export
#'
import_taxids <- function(x) {
    fname <- paste0('extdata/', x, '_taxids.tsv')
    file_path <- system.file(fname, package = 'bugphyzzAnalyses')
    read.table(file_path, header = TRUE, sep = '\t', row.names = 1)
}

#' Create DataTable
#'
#' \code{craeteDT} creates a data table.
#'
#' Reference for this code: https://www.r-bloggers.com/2019/12/vignette-downloadable-tables-in-rmarkdown-with-the-dt-package/
#' The link was valid on Dec 2, 2022.
#'
#' @param x A data frame.
#'
#' @return A DataTable output.
#'
#' @export
#'
createDT <- function(x){
    DT::datatable(
        x, extensions = 'Buttons', rownames = FALSE,
         options = list(
            dom = 'Brti',
            buttons = list(
                list(extend = 'copy', title = NULL),
                list(extend = 'csv', title = NULL),
                list(extend = 'excel', title = NULL),
                list(extend = 'pdf', title = NULL),
                list(extend = 'print', title = NULL)
            )
        )
    )
}

#' Import the NYCHANES file
#'
#' \code{imporNYCHANES} imports the annotated NYCHANES file from the
#' original paper and repo. The output is already as a TreeSummarizedExperiment.
#' The imported TSE only contains data for 'Cigarette'  and 'Never Smokers'.
#' Some filtering might be necessary since many samples were removed from the
#' full study.
#'
#' @return A TreeSummarizedExperiment
#' @export
#'
importNYCHANES <- function() {
    fname <- system.file('extdata/nychanes.RDS', package = 'bugphyzzAnalyses')
    readRDS(file = fname)
}

#' Import biosis data
#'
#' \code{importBiosis} imports manually curated annotations used in the
#' paper of Calgaro and NYCHANES. It could be useful for comparison of
#' annotations
#'
#' @return A data.frame
#' @export
#'
importBiosis <- function() {
    fname <- system.file(
        'extdata/biosis.tsv', package = 'bugphyzzAnalyses'
    )
    read.table(fname, header = TRUE, sep = '\t')
}

#' Import taxids for MicrobiomeBenchmarkData
#'
#' \code{importTaxids} imports a data.frame for MicrobiomeBenchmarkData
#'
#' @param x Character string
#'
#' @return A data.frame.
#' @export
#'
importTaxids <- function(x = 'HMP_2012_16S_gingival_V35_taxids') {
    x <- paste0('extdata/', x, '.tsv')
    fname <- system.file(x , package = 'bugphyzzAnalyses')
    df <- read.table(fname, header = TRUE, row.names = 1, sep = '\t')
    dplyr::select(df, -.data$full_taxon_name, -.data$taxon_annotation)
}

#' My DataTable
#'
#' \code{myDT} prints a datatable
#'
#' @param df A data.frame.
#' @param cap A caption.
#' @param ap = ADJ.PVAL threshold.
#'
#' @return A data.table.
#' @export
#'
myDT <- function(df, cap = 'Table. Caption...', ap = 0.1) {
    colnames(df) <- sub('GENE', 'BUG', colnames(df))
    df$BUG.SET <- ifelse(df$ADJ.PVAL < ap, paste0(df$BUG.SET, '*'), df$BUG.SET)
    df <- dplyr::arrange(df, .data$PVAL)
    DT::datatable(
        data = df,
        filter = 'top',
        rownames = FALSE,
        extensions = 'Buttons',
        options = list(
            dom = 'Bft',
            buttons = list('copy', 'print'),
            iDisplayLength = 10,
            keys = TRUE,
            autoWidth = TRUE
        ),
        caption = htmltools::tags$caption(
            style = 'caption-side: top; text-align: left;',
           cap
        )
    )
}

#' Normalize assay with limma::voom
#'
#' \code{limmaVoom} normalizes an assay (the first one) from a
#' (Tree)SummarizedExperiment with lima::voom.
#'
#' @param tse A (Tree)SummarizedExperiment.
#'
#' @return A (Tree)SummarizedExperiment.
#' @export
#'
limmaVoom <- function(tse) {
    df <- data.frame(SummarizedExperiment::colData(tse))
    design <- stats::model.matrix(~ GROUP, data = df)
    assay_voom1 <- limma::voom(
        SummarizedExperiment::assay(tse), design = design, plot = FALSE
    )
    SummarizedExperiment::assay(tse) <- assay_voom1$E
    class(SummarizedExperiment::assay(tse)) <- "matrix"
    return(tse)
}

#' Cacl prediciton stats
#'
#' \code{calcPredStats} calculates prediction stats related to AUC ROC
#'
#' @param df A data.frame
#'
#' @return A data.frame
#' @export
#'
calcPredStats <- function(df) {
    df <- df |>
        dplyr::count(label)

    TP <- df |>
        dplyr::filter(label == 'TP') |>
        dplyr::pull(n)

    if (!length(TP)) {
        TP <- 0
    }

    TN <- df |>
        dplyr::filter(label == 'TN') |>
        dplyr::pull(n)

    if (!length(TN)) {
        TN <- 0
    }

    FP <- df |>
        dplyr::filter(label == 'FP') |>
        dplyr::pull(n)

    if (!length(FP)) {
        FP <- 0
    }

    FN <- df |>
        dplyr::filter(label == 'FN') |>
        dplyr::pull(n)

    if (!length(FN)) {
        FN <- 0
    }

    tibble::tribble(
        ~ stat, ~ value,
        "sensitivity (recall)", TP / (TP + FN),
        "specificity", TN / (TN + FP),
        "precision", TP / (TP + FP),
        "FDR", FP / (FP + TP),
        "FPR", FP / (FP + TN)
    ) |>
        dplyr::mutate(value = ifelse(is.na(value), 0, value))

}





#' Print datatable with my default options
#'
#' \code{myDataTable} prints a DT::datatable with my default options.
#'
#' @param dat A data.frame.
#' @param page_len An integer. Number of rows to print.
#'
#' @return A datatable object
#' @export
#'
myDataTable <- function(dat, page_len = NULL) {
    if (!is.null(page_len)) {
        output <- dat |>
            DT::datatable(
                rownames = FALSE,
                extensions = "Buttons",
                filter = "top",
                options = list(
                    dom = 'lBfrtip',
                    paging = TRUE,
                    pageLength = page_len,
                    buttons = list(
                        list(extend = "copy", title = NULL),
                        list(extend = "print", title = NULL),
                        list(extend = "csv", title = NULL)
                    ),
                    scrollX = TRUE
                )
            )
    } else {
        output <- dat |>
            DT::datatable(
                rownames = FALSE,
                extensions = "Buttons",
                filter = "top",
                options = list(
                    dom = 'lBfrtip',
                    # lengthMenu = c(10, 20, 50, 100),
                    buttons = list(
                        list(extend = "copy", title = NULL),
                        list(extend = "print", title = NULL),
                        list(extend = "csv", title = NULL)
                    ),
                    scrollX = TRUE
                )
            )
    }
    return(output)
}

#' Elbows
#'
#' \code{elbows} returns the prevalence threshold for TypicalMicrobiomeSignatures
#' based on the elbow of the curve method. See vignette.
#'
#' @param ag Age group. A character string. "adult" or "child"
#'
#' @return Vector of type double with thresholds.
#' @export
#'
#' @examples
#'
#' elbows()
#' elbows("child")
#'
elbows <- function(ag = "adult") {
    c(
        # feces_genus = 0.01,
        # feces_species = 0.01,
        # mouth_genus = 0.01,
        # mouth_species = 0.01,
        # skin_genus = 0.01,
        # skin_species = 0.01,
        # vagina_genus = 0.01,
        # vagina_species = 0.01
    )
    adult <- c(
        feces_genus = 0.04097764,
        feces_species = 0.04430577,
        mouth_genus = 0.01363636,
        mouth_species = 0.01363636,
        skin_genus = 0.07048458,
        skin_species = 0.07488987,
        vagina_genus = 0.02105263,
        vagina_species = 0.01052632
    )
    child <- c(
        feces_genus = 0.06133333,
        feces_speies = 0.04977778,
        mouth_genus = 0.04938272,
        mouth_species = 0.04938272,
        skin_genus = 0.01,
        skin_species = 0.01
    )
    if (ag == "adult") {
        return(adult)
    } else if (ag == "child") {
        return(child)
    }
}


## Non exported functions ----------------------------------------------------

#' Import typical microbiome signatures
#'
#' \code{importTMS} Import all of the typical microbiome signatures in
#' a tidy format.
#'
#' @return A data.frame.
#' @export
#'
importTMS <- function() {
    bfc <- .getCache()
    rid <- BiocFileCache::bfcquery(bfc, query = "tms", field = "rname")$rid
    if (!length(rid)) {
        tms <- .downloadTMS()
        fpath <- BiocFileCache::bfcnew(bfc, rname = "tms", ext = ".tsv")
        readr::write_tsv(tms, file = fpath)
    } else {
        message("Using cached file...")
        fpath <- BiocFileCache::bfcpath(x = bfc, rids = rid)
        tms <- readr::read_tsv(fpath, show_col_types = FALSE)
    }
    return(tms)
}

.getCache <- function() {
    cache <- tools::R_user_dir("bugphyzzAnalyses", which = "cache")
    BiocFileCache::BiocFileCache(cache)
}

.downloadTMS <- function() {
    message("Downloading typicial microbiome signatures...")
    url <- "https://zenodo.org/records/7622129/files/waldronlab/TypicalMicrobiomeSignaturesExports-v1.0.1.zip?download=1"
    temp_dir <- tempdir()
    temp_file <- file.path(temp_dir, "tms.zip")
    download.file(url = url, destfile = temp_file)
    unzip(temp_file, exdir = temp_dir, junkpaths = TRUE)
    csv_files <- list.files(temp_dir, pattern = "csv", full.names = TRUE)
    l <- purrr::map(csv_files,  ~ {
        .x |>
            utils::read.csv() |>
            tidyr::pivot_longer(
                names_to = "Body site", values_to = "Prevalence",
                cols = 3:tidyselect::last_col()
            ) |>
            dplyr::mutate(
                `Body site` = sub("_(species|genus)_prevalence$", "", `Body site`)
            )
    })
    names(l) <- sub("^.*matrix_(.*).csv$", "\\1", csv_files)
    tms <- dplyr::bind_rows(l, .id = "rank_agegroup") |>
        tidyr::separate(
            col = "rank_agegroup", into = c("Rank", "Age group"),
            sep = "_", remove = TRUE
        ) |>
        dplyr::relocate(
            `Age group`, `Rank`, `NCBI ID` = NCBI, `Taxon name` = name,
            `Body site`, Prevalence
        ) |>
        dplyr::mutate(
            `Body site` = dplyr::case_when(
                `Body site` == "stool" ~ "feces",
                `Body site` == "oralcavity" ~ "mouth",
                TRUE ~ `Body site`
            )
        )
}
