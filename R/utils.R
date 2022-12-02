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
                 extend = c('copy', 'csv', 'excel', 'pdf', 'print'),
                 title = NULL
            )
        )
    )
}
