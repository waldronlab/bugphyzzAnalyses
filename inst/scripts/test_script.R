library(bugphyzz)
library(taxPPro)

phys <- physiologies(remove_false = TRUE)

## Only use physiologies with categorical values.
## For numeric and rages use thresholds instead of real values.
phys <- phys[vapply(phys, function(x) unique(x$Attribute_type) == 'categorical', logical(1))]
categorical <- lapply(
    phys, preSteps, tax.id.type = 'NCBI_ID', remove_false = TRUE
)

exclude_names <- c(
    'antimicrobial resistance'
)

cat_names <- names(categorical)
propagation <- vector('list', length(categorical))
names(propagation) <- cat_names
for (i in seq_along(propagation)) {
    current_name <- cat_names[i]
    if (current_name %in% exclude_names)
        next
    message('>>> Propagating ',  current_name, '. <<<')
    stat_time <- Sys.time()
    propagation[[i]] <- propagate(categorical[[i]])
    end_time <- Sys.time()
    message('>>> Done ', difftime(end_time, stat_time), '. <<<')
}
