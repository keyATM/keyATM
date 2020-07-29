# keyATM 0.3.1
### Major changes
* Changes related to release of dplyr 1.0.1


# keyATM 0.3.0
[Roadmap](https://github.com/keyATM/keyATM/projects/2) for this version.

### Major changes
* Use the Highest Density Interval as a default (`method = "hdi"`) in `plot.strata_doctopic()`, `plot_timetrend()`, and `plot_pi()`. The previous version uses the Equal-tailed Interval (`method = "eti"`).
* Add `read_keywords` for reading dictionary files (e.g. YAML, LIWC).
* Add the `predict()` function for the covariate keyATM (thank you [Sanja Hajdinjak](https://sanjahajdinjak.wordpress.com/) for the suggestion!).
* Detailed options for standardization in the covariate keyATM
    * The `standardize` option in `model_settings` argument of the `keyATM()` function now takes one of `"all"`, `"none"`, or `"non-factor"` (default).
    * `"all"` standardizes all covariates (except the intercept), `"none"` does not standardize any covariates, and `"non-factor"` standardizes non-factor covariates.
    * In previous versions, this option takes either `TRUE` (default, standardizing all covariates) or `FALSE`.
* A bug fix in the `by_strata_DocTopic()` function.
* The output of the `keyATM()` includes the index of documents used for fitting (this will be useful if the input includes documents with zero length).
* Add a `progress_bar` option in the `keyATM_read()` function (thank you [Jae Yeon Kim](https://jaeyk.github.io/) for the suggestion!)

### Bug fix
* Fix checking time index input (thank you [Jae Yeon Kim](https://jaeyk.github.io/) for pointing out this issue!)

# keyATM 0.2.2
### Major changes
* Updates for dplyr 1.0.0
* Update tests

# keyATM 0.2.1
### Major changes
* Temporary update `test-Initialization.R` to deal with some errors.

# keyATM 0.2.0
[Roadmap](https://github.com/keyATM/keyATM/projects/1) for this version.

### Major changes
* Update the `by_strata_DocTopic()` function.
* Make examples runnable (thank you [Chung-hong Chan](https://github.com/chainsawriot) for the suggestion!)
* Speed up (about 15% faster)
* `save_fig()` function
* Automatically drops documents with length 0, raising a warning (thank you [Francesco Grossetti](https://github.com/contefranz) for the suggestion!)
* Update `plot.strata_doctopic()`: showing by topic by default (thank you [Soichiro Yamauchi](https://soichiroy.github.io/) for the suggestion!)

### Bug fix
* `weightedLDA()` without specifying the number of iterations ([Chung-hong Chan](https://github.com/chainsawriot) independently reported this bug, thank you!)
* log-likelihood of dynamic models 
* saving figures
* topic labels when there is no keyword topic
* `summary.strata_doctopic()`: the last topic is removed when the number of no-keyword topic is 0 (thank you [Emma Ebowe](https://gov.harvard.edu/people/emma-ebowe) for pointing out this issue!)

# keyATM 0.1.0
### Major changes
* The first CRAN version
* Organize functions into a package
* Add keyATM Label
* Replace `hashmap` with `fastmap`
* Thank you [Santiago Olivella](https://github.com/solivella) for finding several bugs!

# keyATM 0.0.7
### Major changes
* We have a new syntax (this version does not support objects made in older keyATM)
* Faster read functions
* Memory efficiency

# keyATM 0.0.6
### Major changes
* Add keyATM Dynamic

# keyATM 0.0.5
### Major changes
* Add keyATM Covariate

# keyATM 0.0.4
### Major changes
* Faster estimation

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
