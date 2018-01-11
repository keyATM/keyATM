library(quanteda)
data(corpus_bara_para)
corp <- corpus_subset(corpus_bara_para,
                      speaker != "Mrs Anne Kerr")

file <- system.file("extdata", "bara-et-al.ykd", package = "topicdict")
ddict <- dictionary(file = file)
ddict_words <- as.vector(unlist(ddict))
# It would be a bad idea to remove gendered pronouns as stop words because
# of the debate topic, so we don't do that
estops <- stopwords("english") # re-gendered stopwords
# stops minus some gender words
stops <- estops[!(estops %in% c("her", "his", "he", "she", "their", "our"))]
# Make a document term matrix and take a copy that's tf-idf transformed.
# we'll use that to choose informative words.
dtm_orig <- dfm(corp, remove = stops,
                remove_numbers = TRUE,
                remove_punct = TRUE,
                remove_symbols = TRUE,
                remove_separators = TRUE,
                remove_hyphens = TRUE)
tfidf_dtm_orig <- tfidf(dtm_orig)
# The 50 most informative words in thus corpus according to tf-idf score
# sort(colMeans(tfidf_dtm_orig), decreasing = TRUE)[1:50]

informative_wds <- data.frame(tfidf = colMeans(tfidf_dtm_orig))
informative_wds <- informative_wds[order(informative_wds$tfidf,
                                         decreasing = TRUE), ,drop = FALSE]
# Construct the paper's topic counts using their dictionary.
# This is what Bara et al think the correct, person summed, topic counts ought to be
topics_orig_dict <- dfm(corp, remove = stops,
                        remove_numbers = TRUE,
                        remove_punct = TRUE,
                        remove_symbols = TRUE,
                        remove_separators = TRUE,
                        remove_hyphens = TRUE,
                        dictionary = ddict)
tt <- data.frame(topics_orig_dict,
                 speaker = docvars(corp, "speaker"))
colnames(tt) <- unlist(gsub("incoming.txt.", "", colnames(tt)))
tts <- summarise_at(group_by(tt, speaker),
                    vars(advocacy, legal, medical,
                         moral, procedural, social),
                    sum)
tts <- as.data.frame(tts)
tt_prop <- data.frame(tts[,2:7] / rowSums(tts[,2:7]),
                      row.names = tts$speaker)
# tt_prop

# Now construct the top 10 seed words using a tf-idf heuristic
# and write to a file
seeds <- lapply(ddict$incoming.txt, function(x){
  candidates <- intersect(rownames(informative_wds), x)
  inf <- informative_wds[candidates, ,drop = FALSE]
  rownames(inf)[order(inf$tfidf, decreasing = TRUE)[1:10]]
})
lines <- as.vector(apply(t(as.data.frame(seeds)), 1,
                         paste, collapse = " "))
writeLines(c("##", lines), con = file("./bara_seeds.txt"))
## Write out the pre-processed texts to a file
tt <- tokens(corp, remove_numbers = TRUE,
             remove_punct = TRUE,
             remove_symbols = TRUE,
             remove_separators = TRUE,
             remove_hyphens = TRUE)
tt <- lapply(tt, function(x) {
  x <- char_tolower(x)
  paste(x[!(x %in% stops)], collapse = " ")
})
names(tt) <- gsub(" ", "_",
                  paste0(make.unique(docvars(corp, "speaker")), ".txt"))
if (!dir.exists("bara_paras"))
  dir.create("bara_paras")
devnull <- lapply(names(tt), function(x){
  writeLines(tt[[x]], con = file.path("bara_paras", x))
})
## assign these document names later
docnames <- names(tt)
