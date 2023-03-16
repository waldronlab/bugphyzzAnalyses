#' Prepare data for propagation
#'
#' \code{preSteps} prepares data for propagation.
#'
#' @param df A data.frame.
#' @param tax.id.type A character string. Default: 'NCBI_ID'.
#' @param remove_false TRUE or FALSE. Default TRUE.
#'
#' @return A data.frame.
#' @export
#'
preSteps <- function(df, tax.id.type = 'NCBI_ID', remove_false = TRUE) {
    x <- taxPPro::prepareData(
        df, tax.id.type = tax.id.type, remove_false = remove_false
    )
    if (is.null(x))
        return(NULL)
    output <- taxPPro::prepareData2(x)
    return(output)
}
