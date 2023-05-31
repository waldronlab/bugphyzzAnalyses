library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(tidySummarizedExperiment)

NYC_HANES <- annotateFactors(loadQiimeData())
nh <- mia::makeTreeSEFromPhyloseq(NYC_HANES)

x <- nh |>
    filter(
        smokingstatus == 'Alternative smoker',
        CIGARETTES == 'Yes'
    )
colData(x)$Burklab_ID
