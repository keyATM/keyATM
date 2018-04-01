
#' Preprocess a dictionary to make it usable as seed words
#'
#' This function take a dictionary, collapses to the top layer of
#' categories, adds all the corpus matches for any wildcarded pattern,
#' and removes anything that does not occur more than \code{min.freq}
#' times. The result is a dictionary suitable for use as seeds in
#' \code{topicdict_train}.
#'
#' @param dict A quanteda dictionary, or the filename of something quanteda can
#'             read in as a dictionary, or a named list which can be treated as
#'             a dictionary
#' @param corpus A quanteda corpus
#' @param min.freq The fewest times a word can appear and remain in the dictionary
#' @param ... Extra arguments to give to \code{dfm}
#'
#' @return A dictionary with entries filtered by frequency
#' @export
#' @importFrom stats setNames
#'
#' @examples
#' dictfile <- system.file("extdata/laver-garry-ajps.ykd", package="topicdict")
#' # Use just the economics subtree of this dictionary
#' lavergarry <- quanteda::dictionary(file = dictfile)
#' dict <- lavergarry[['Laver and Garry']][["State in Economy"]]
#' data("corpus_uk_platforms", package = "quanteda")
#' new_dict <- preprocess_dictionary(dict, corpus_uk_platforms, min.freq = 5)
preprocess_dictionary <- function(dict, corpus, min.freq = 5, ...){
  if (is.character(dict))
    dict <- dictionary(file = dict)
  else if (is.list(dict))
    dict <- dictionary(x = dict)

  dcats <- setNames(lapply(names(dict),
                           function(n){ as.character(unlist(dict[[n]])) }),
                    names(dict))
  dtm_orig <- dfm(corpus, ...)
  # words that are frequent enough
  tstats <- quanteda::textstat_frequency(dtm_orig)
  target_vocab <- tstats$feature[tstats$frequency >= min.freq]
  # process patterns, expanding globs and filtering for frequency
  expand_patterns <- function(x){
    starred <- grep("*", x, fixed = TRUE, value = TRUE)
    globs <- gsub("*", "", starred, fixed = TRUE)
    expanded_targets <- unlist(lapply(globs, grep, x = target_vocab, value = TRUE))
    non_starred_targets <- intersect(setdiff(x, starred), target_vocab)
    union(expanded_targets, non_starred_targets)
  }
  dcats <- lapply(dcats, expand_patterns)
  dictionary(dcats)
}
#
#
# stops <- estops[!(estops %in% c("her", "his", "he", "she", "their", "our"))]
# # Make a document term matrix and take a copy that's tf-idf transformed.
# # we'll use that to choose informative words.
# dtm_orig <- dfm(corp, remove = stops,
#                 remove_numbers = TRUE,
#                 remove_punct = TRUE,
#                 remove_symbols = TRUE,
#                 remove_separators = TRUE,
#                 remove_hyphens = TRUE)
# tfidf_dtm_orig <- tfidf(dtm_orig)
# # The 50 most informative words in thus corpus according to tf-idf score
# # sort(colMeans(tfidf_dtm_orig), decreasing = TRUE)[1:50]
#
# informative_wds <- data.frame(tfidf = colMeans(tfidf_dtm_orig))
# informative_wds <- informative_wds[order(informative_wds$tfidf,
#                                          decreasing = TRUE), ,drop = FALSE]
# # Construct the paper's topic counts using their dictionary.
# # This is what Bara et al think the correct, person summed, topic counts ought to be
# topics_orig_dict <- dfm(corp, remove = stops,
#                         remove_numbers = TRUE,
#                         remove_punct = TRUE,
#                         remove_symbols = TRUE,
#                         remove_separators = TRUE,
#                         remove_hyphens = TRUE,
#                         dictionary = ddict)
# tt <- data.frame(topics_orig_dict,
#                  speaker = docvars(corp, "speaker"))
# colnames(tt) <- unlist(gsub("incoming.txt.", "", colnames(tt)))
# tts <- summarise_at(group_by(tt, speaker),
#                     vars(advocacy, legal, medical,
#                          moral, procedural, social),
#                     sum)
# tts <- as.data.frame(tts)
# tt_prop <- data.frame(tts[,2:7] / rowSums(tts[,2:7]),
#                       row.names = tts$speaker)
# # tt_prop
#
# # Now construct the top 10 seed words using a tf-idf heuristic
# # and write to a file
# seeds <- lapply(ddict$incoming.txt, function(x){
#   candidates <- intersect(rownames(informative_wds), x)
#   inf <- informative_wds[candidates, ,drop = FALSE]
#   rownames(inf)[order(inf$tfidf, decreasing = TRUE)[1:10]]
# })
# lines <- as.vector(apply(t(as.data.frame(seeds)), 1,
#                          paste, collapse = " "))
# writeLines(c("##", lines), con = file("./bara_seeds.txt"))
# ## Write out the pre-processed texts to a file
# tt <- tokens(corp, remove_numbers = TRUE,
#              remove_punct = TRUE,
#              remove_symbols = TRUE,
#              remove_separators = TRUE,
#              remove_hyphens = TRUE)
# tt <- lapply(tt, function(x) {
#   x <- char_tolower(x)
#   paste(x[!(x %in% stops)], collapse = " ")
# })
# names(tt) <- gsub(" ", "_",
#                   paste0(make.unique(docvars(corp, "speaker")), ".txt"))
# if (!dir.exists("bara_paras"))
#   dir.create("bara_paras")
# devnull <- lapply(names(tt), function(x){
#   writeLines(tt[[x]], con = file.path("bara_paras", x))
# })
# ## assign these document names later
# docnames <- names(tt)
