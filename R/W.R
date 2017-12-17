## need wd_id -> wd   = vocab[wd_id+1])
## dictcatid_id -> dictcatname  = names(dict)[[dictcat_id]]
## docid -> docname  = files[docid]

#' Create a list of word indexes W
#'
#' @param files Names of files
#' @param dict a quanteda dictionary or named list of character vectors
#' @param k number of topics
#' @param encoding File encoding (defaults to \code{readLines}'s default)
#' @param ... Additional arguments to \code{quanteda::tokens}. (Lowercasing happens by default).
#'
#' @return A list containing W: a list of vectors of word indexes,
#'         Z: a list of vectors of topic indicators,
#'         X: a list of vectors of seed indicators
#'         vocab: a vector of vocabulary,
#'         dict: tokenized version of the dictionary passed to the function
#' @export
#'
init_XWZ <- function(files, dict, k, encoding = "unknown", ...){
  tokenize_text <- function(x){
    x <- paste(readLines(x, encoding = encoding), collapse = "\n")
    as.character(quanteda::tokens(char_tolower(x), ...))
  }
  # tokenize the dictionary entries and the documents
  dict <- lapply(dict, function(x){ as.character(tokens(x)) })
  toklist <- lapply(files, tokenize_text)

  ## construct W and a vocab list (W entries are 0 based)
  vocab <- unique(unlist(toklist))
  wd_id <- as.numeric(factor(vocab))
  wd_map <- hashmap::hashmap(vocab, as.integer(wd_id - 1)) # because we're going into C++ code
  W <- lapply(toklist, function(x){ wd_map[[x]] })

  # zx_assigner maps seed words to category ids
  seed_wdids <- unlist(lapply(dict, function(x){ wd_map$find(x) }))
  cat_ids <- rep(1:length(dict) - 1, unlist(lapply(dict, length)))
  zx_assigner <- hashmap::hashmap(as.integer(seed_wdids), as.integer(cat_ids))

  # if the word is a seed, 1, else 0
  make_x <- function(x){
    z <- zx_assigner$find(x)
    as.numeric(!is.na(z))
  }
  X <- lapply(W, make_x)

  # if the word is a seed, assign the appropriate (0 starting) Z, else a random Z
  make_z <- function(x){
    zz <- zx_assigner$find(x)
    unseeded <- is.na(zz)
    zz[unseeded] <- sample(1:k - 1, sum(unseeded), replace = TRUE)
    zz
  }
  Z <- lapply(W, make_z)

  list(W = W, Z = Z, X = X, vocab = vocab, files = files, dict = dict)
}


## init_XWZ(list.files("inst/extdata", pattern = "macavity", full.names = TRUE),
##     dict = dictionary(list(mac = c("macavity", "mystery", "cat"),
##                            victims = c("milk", "admiralty", "trellis", "drawings"))),
##     remove_numbers = TRUE, remove_punct = TRUE, remove_symbols = TRUE,
##     remove_separators = TRUE)
