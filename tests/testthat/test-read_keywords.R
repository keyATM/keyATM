require(quanteda)
data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
keyATM_docs <- keyATM_read(bills_dfm)
bills_keywords <- keyATM_data_bills$keywords

test_that("normal usage", {
    bill_dictionary <- quanteda::dictionary(bills_keywords)
    x <- read_keywords(dictionary = bill_dictionary, docs = bills_keywords)
    expect_equal(x$Drug, "drug")
    ## defensive programming
    expect_error(read_keywords())
    expect_warning(read_keywords(dictionary = bill_dictionary))
})

test_that("resolving and split", {
    ## a simplfied example from koheiw/quanteda.seededlda
    corp <- corpus(c("air force 1", "soldier soldiers navy army politician", "party members", "parliament is politicians' party politics leaders leadership"))
    toks <- tokens(corp, remove_punct = TRUE)
    dfmt <- dfm(toks) %>% 
        dfm_select("^[A-Za-z]+$", valuetype = "regex") %>% 
        dfm_remove(stopwords('en'))
    docs <- keyATM_read(texts = dfmt)
    smalldictionary <- dictionary(list(politics = c("parliament*", "congress*", "party leader*", "party member*", "voter*", "lawmaker*", "politician*"), military = c("military", "soldier*", "air force", "marine", "navy", "army")))
    ## default: split = TRUE
    x <- read_keywords(dictionary = smalldictionary, docs = docs)
    expect_true("soldiers" %in% x$military)
    expect_true("air" %in% x$military)
    y <- read_keywords(dictionary = smalldictionary, docs = docs, split = FALSE)
    expect_true(!"air" %in% y$military)
})

test_that("file I/O", {
    ## the example of the function
    skip_on_cran() ; skip_on_os("linux")
    dictfile <- tempfile()
    download.file("http://bit.ly/37cV95h", dictfile)
    ## first position is file.
    x <- read_keywords(dictfile, docs = keyATM_docs, format = "LIWC")
    expect_true("terrorism" %in% x$IngroupVice)
})

test_that("Integration",{
    skip_on_cran() ; skip_on_os("linux")
    ## exact example: koheiw/quanteda.seededlda
    rdsfile <- tempfile()
    download.file("https://github.com/koheiw/quanteda.seededlda/raw/master/tests/data/data_corpus_sputnik.RDS", rdsfile)

    corp <- readRDS(rdsfile)
    toks <- tokens(corp, remove_punct = TRUE)
    dfmt <- dfm(toks) %>% 
    dfm_select("^[A-Za-z]+$", valuetype = "regex") %>% 
    dfm_remove(stopwords('en')) %>% 
    dfm_trim(min_termfreq = 0.90, termfreq_type = "quantile", 
             max_docfreq = 0.1, docfreq_type = "prop")
    docs <- keyATM_read(texts = dfmt)
    dictfile <- tempfile(fileext = ".yml")
    download.file("https://raw.githubusercontent.com/koheiw/quanteda.seededlda/master/tests/data/topics.yml", dictfile)
    kw <- read_keywords(dictfile, docs = docs)
    res <- keyATM(docs, no_keyword_topics = 0, keywords = kw, model = "base", options = list(seed = 1234, iterations = 20))
    ## syria is related to military
    expect_true("syria" %in% top_words(res, n = 20, show_keyword = FALSE)[,"5_military"])
    ## German is not related to military
    expect_true(!"german" %in% top_words(res, n = 20, show_keyword = FALSE)[,"5_military"])
})
