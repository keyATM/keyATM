

#' Create a list of word indexes W
#'
#' @param files Names of files
#' @param encoding File encoding (defaults to \code{readLines}'s default)
#' @param ... Additional arguments to \code{quanteda::tokens}
#'
#' @return A list containing W: a list of word indexes, word_to_id: a hashmap of words
#'         to 0-starting integers, doc_to_id:a hashmap of elements of files to
#'         0-starting integers
#' @export
#'
create_W <- function(files, encoding = "unknown", ...){
  tokenize_text <- function(x){
    x <- paste(readLines(x, encoding = encoding), collapse = "\n")
    as.character(quanteda::tokens(x, ...))
  }
  toklist <- lapply(files, tokenize_text)
  vocab <- unique(unlist(toklist))
  wd_id <- as.numeric(factor(vocab))
  wd_map <- hashmap::hashmap(vocab, wd_id - 1) # because we're going into C++ code
  doc_map <- hashmap(files, 1:length(files) - 1)
  list(W = lapply(toklist, function(x){ map.find(x) }),
       word_to_id = wd_map,
       doc_to_id = doc_map)
}
