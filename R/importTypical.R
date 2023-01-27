importTypical <- function() {
    temp_dir <- paste0(tempdir(), '/TypicalMicrobiomeSignatures')
    dir.create(temp_dir, showWarnings = FALSE)
    doi <- "10.5281/zenodo.7544550"
    zen4R::download_zenodo(doi = doi, path = temp_dir, quiet = TRUE)
    zip_file <- list.files(temp_dir, full.names = TRUE)
    base_dir <- sub('\\.zip', '', zip_file)
    utils::unzip(zip_file, exdir = base_dir, )
    extracted_dir <- paste0(base_dir, '/', list.files(base_dir))
    fnames <- grep(
        '\\.csv', list.files(extracted_dir, full.names = TRUE), value = TRUE
    )
    tms <- map(fnames, ~ {
        df <- read.csv(.x)
    })
    names(tms) <- sub('^.*/matrix_(.*)\\.csv$', '\\1', fnames)
    output <- vector('list', length(tms))
    for (i in seq_along(output)) {
        output[[i]] <- tms[[i]] |>
            tidyr::pivot_longer(
                cols = tidyr::ends_with('_prevalence'), names_to = 'body_site',
                values_to = 'pravalence'
            ) |>
            dplyr::mutate(body_site = sub('_prevalence$', '', body_site)) |>
            tidyr::separate(
                col = 'body_site', into = c('body_site', 'rank'), sep = '_'
            ) |>
            dplyr::rename(
                taxid = .data$NCBI, taxon_name = .data$name
            ) |>
            dplyr::mutate(age_range = sub('^.*_', '', names(tms)[i]))
    }
    data <- do.call(rbind, output)
    data <- data[data$pravalence > 0.01,]
    return(data)
}
