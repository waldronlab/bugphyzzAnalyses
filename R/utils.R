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

#' Import tidy Typical Microbiome Signatures
#'
#' \code{importTidyTMS} imports prevalence data from the
#' TypicalMicrobiomeSignatures project.
#'
#' @param prevalence_threshold A numeric value: 0.01.
#'
#' @return A data.frame.
#' @export
#'
importTidyTMS <- function(prevalence_threshold = 0.01) {
    base_dir <- tools::R_user_dir('bugphyzzAnalyses', which = 'cache')
    zip_file <- .downloadTMSZip()
    unzip(zip_file, exdir = base_dir)
    extracted_dir <- paste0(base_dir, '/', list.files(base_dir))
    fnames <- grep(
        '\\.csv', list.files(extracted_dir, full.names = TRUE), value = TRUE
    )
    tms <- purrr::map(fnames, read.csv)
    names(tms) <- sub('^.*/matrix_(.*)\\.csv$', '\\1', fnames)
    output <- vector('list', length(tms))
    for (i in seq_along(output)) {
        output[[i]] <- tms[[i]] |>
            tidyr::pivot_longer(
                cols = ends_with('_prevalence'), names_to = 'body_site',
                values_to = 'prevalence'
            ) |>
            dplyr::mutate(body_site = sub('_prevalence$', '', body_site)) |>
            tidyr::separate(
                col = 'body_site', into = c('body_site', 'rank'), sep = '_'
            ) |>
            dplyr::rename(
                taxid = NCBI, taxon_name = name
            ) |>
            dplyr::mutate(age_range = sub('^.*_', '', names(tms)[i]))
    }
    typical <- do.call(rbind, output)
    typical$prevalence <- round(typical$prevalence, 2)
    typical <- typical[typical$prevalence >= 0.01,]
    return(typical)
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

.get_cache <- function() {
    cache <- tools::R_user_dir("bugphyzzAnalyses", which="cache")
    BiocFileCache::BiocFileCache(cache)
}

.downloadTMSZip <- function(verbose = FALSE) {
    fileURL <- 'https://zenodo.org/records/7544550/files/waldronlab/TypicalMicrobiomeSignaturesExports-v1.0.0.zip?download=1'
    bfc <- .get_cache()
    rid <- BiocFileCache::bfcquery(bfc, "TypicalMicrobiomeSignatures.zip", "rname")$rid
    if (!length(rid)) {
        if( verbose )
            message( "Downloading GENE file" )
            rid <- names(BiocFileCache::bfcadd(bfc, "TypicalMicrobiomeSignatures.zip", fileURL))
    }
    BiocFileCache::bfcrpath(bfc, rids = rid)
}
