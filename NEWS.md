# keyATM 0.5.1

### Minor changes
* Adding options to set the hyperparameters `eta`.
* Documentation updates to pass CRAN checks.


# keyATM 0.5.0
### Major changes
* Migrating to C++17 to follow [the new CRAN check](https://developer.r-project.org/blosxom.cgi/R-devel/2023/01/27#n2023-01-27). We edited `shuffled_indexes()` that internally used `std::random_shuffle()`. This change does not guarantee backward compatibility across all platforms.
* Using the package `cli` instead of the base R `message()` and `warning()` functions.
* A new feature to resume the iteration.
* Stopped support for the `label` model in `keyATM()` (it was an experimental feature).

### Minor changes
* Supporting [tidyselect 1.2.0](https://www.tidyverse.org/blog/2022/10/tidyselect-1-2-0/#using-data-inside-selections). Updating some internal functions.
* Including the state information of HMM in the `plot_timetrend()` output (thank you [@WenHanGao](https://github.com/WenHanGao) for the suggestion in [#188](https://github.com/keyATM/keyATM/issues/188)).

### Bug fix
* A bug fix in `summary.keyATM_docs()`.

# keyATM 0.4.2

### Major changes
* Adding the `plot_topicprop()` function.
* Implementation of semantic coherence diagnostics from Mimno et al. (2011) via the `semantic_coherence()` function (thanks to [Seo-young Silvia Kim](https://sysilviakim.com/) for your [suggestion](https://github.com/keyATM/keyATM/pull/191)).

### Minor changes
* Pure R text loading to address issues related to UTF-8 encoding in Windows ([#189](https://github.com/keyATM/keyATM/issues/189))

# keyATM 0.4.1

### Minor changes
* Remove an unused argument (`width`) in the `plot_timetrend()` function.
* Use `Rcpp::message()` if `verbose = TRUE`.
* Completely remove `parallel::mclapply`.

### Bug fix
* `by_strata_DocTopic()` takes the correct arguments ([#180](https://github.com/keyATM/keyATM/issues/180), thank you [@pmeiners](https://github.com/pmeiners) for reporting this!).

# keyATM 0.4.0
[Roadmap](https://github.com/keyATM/keyATM/projects/3) for this version.

### Major changes
* Implementation of Polya-Gamma covariate keyATM.
* Use `future.apply` instead of `parallel` (no backward compatibility if you use the `init_parallel` option).
* The `keyATM_read()` function returns a list of objects (e.g., text and document index).
* An option to store document names in a quanteda dfm object. The `keep_docnames` option in the `keyATM_read()` function (thank you [Morgan 'Les' DeBusk-Lane](https://github.com/debusklaneml) for the suggestion!).
* An option to split a dfm to choose keywords with an unsupervised topic model.

### Bug fix
* Using just first 58 speeches of inaugural corpus in test (thank you [Ken Benoit](https://github.com/kbenoit) for catching this!).

# keyATM 0.3.1
### Major changes
* Changes related to release of dplyr 1.0.1.


# keyATM 0.3.0
[Roadmap](https://github.com/keyATM/keyATM/projects/2) for this version.

### Major changes
* Use the Highest Density Interval as a default (`method = "hdi"`) in `plot.strata_doctopic()`, `plot_timetrend()`, and `plot_pi()`. The previous version uses the Equal-tailed Interval (`method = "eti"`).
* Add `read_keywords` for reading dictionary files (e.g. YAML, LIWC).
* Add the `predict()` function for the covariate keyATM (thank you [Sanja Hajdinjak](https://sanjahajdinjak.wordpress.com/) for the suggestion!).
* Detailed options for standardization in the covariate keyATM:
    * The `standardize` option in `model_settings` argument of the `keyATM()` function now takes one of `"all"`, `"none"`, or `"non-factor"` (default).
    * `"all"` standardizes all covariates (except the intercept), `"none"` does not standardize any covariates, and `"non-factor"` standardizes non-factor covariates.
    * In previous versions, this option takes either `TRUE` (default, standardizing all covariates) or `FALSE`.
* A bug fix in the `by_strata_DocTopic()` function.
* The output of the `keyATM()` includes the index of documents used for fitting (this will be useful if the input includes documents with zero length).
* Add a `progress_bar` option in the `keyATM_read()` function (thank you [Jae Yeon Kim](https://jaeyk.github.io/) for the suggestion!).

### Bug fix
* Fix checking time index input (thank you [Jae Yeon Kim](https://jaeyk.github.io/) for pointing out this issue!).

# keyATM 0.2.2
### Major changes
* Updates for dplyr 1.0.0.
* Update tests.

# keyATM 0.2.1
### Major changes
* Temporary update `test-Initialization.R` to deal with some errors.

# keyATM 0.2.0
[Roadmap](https://github.com/keyATM/keyATM/projects/1) for this version.

### Major changes
* Update the `by_strata_DocTopic()` function.
* Make examples runnable (thank you [Chung-hong Chan](https://github.com/chainsawriot) for the suggestion!).
* Speed up (about 15% faster).
* `save_fig()` function.
* Automatically drops documents with length 0, raising a warning (thank you [Francesco Grossetti](https://github.com/contefranz) for the suggestion!).
* Update `plot.strata_doctopic()`: showing by topic by default (thank you [Soichiro Yamauchi](https://soichiroy.github.io/) for the suggestion!).

### Bug fix
* `weightedLDA()` without specifying the number of iterations ([Chung-hong Chan](https://github.com/chainsawriot) independently reported this bug, thank you!).
* Log-likelihood of dynamic models.
* Saving figures
* Topic labels when there is no keyword topic.
* `summary.strata_doctopic()`: the last topic is removed when the number of no-keyword topic is 0 (thank you [Emma Ebowe](https://www.gov.harvard.edu/directory/emma-ebowe/) for pointing out this issue!).

# keyATM 0.1.0
### Major changes
* The first CRAN version.
* Organize functions into a package.
* Add keyATM Label.
* Replace `hashmap` with `fastmap`.
* Thank you [Santiago Olivella](https://github.com/solivella) for finding several bugs!

# keyATM 0.0.7
### Major changes
* We have a new syntax (this version does not support objects made in older keyATM).
* Faster read functions.
* Memory efficiency.

# keyATM 0.0.6
### Major changes
* Add keyATM Dynamic.

# keyATM 0.0.5
### Major changes
* Add keyATM Covariate.

# keyATM 0.0.4
### Major changes
* Faster estimation.

# keyATM 0.0.3
### Major changes
* This is the first stable version.

# keyATM 0.0.2
### Major changes
* This version implements weighted model.

# keyATM 0.0.1
### Major changes
* This is the first release of keyATM.
* It includes the Base model and the first version of the Covariate model.
