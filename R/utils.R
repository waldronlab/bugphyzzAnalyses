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
