
## Create or get the BiocFileCache object for the package
.get_cache <- function() {
    cache <- tools::R_user_dir('bugphyzzAnalyses', which = 'cache')
    BiocFileCache::BiocFileCache(cache)
}

#' Convert con
