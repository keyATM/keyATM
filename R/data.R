#' Bills data
#' 
#' @format A list with following objects:
#' \describe{
#'   \item{doc_dfm}{A \code{quanteda} dfm object of 140 documents. The text data is a part of the Congressional Bills scraped from \url{https://www.congress.gov}.}
#'   \item{cov}{An integer vector which takes one if the Republican proposed the bill.}
#'   \item{keywords}{A list of length 4 which contains keywords for four selected topics.}
#'   \item{time_index}{An integer vector indicating the session number of each bill.}
#'   \item{labels}{An integer vector indicating 40 labels.}
#'   \item{labels_all}{An integer vector indicating all labels.}
#' }
#' @source \url{https://www.congress.gov}
"keyATM_data_bills"
