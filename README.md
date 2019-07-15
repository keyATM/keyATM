# keyATM

The topic dictionary code, broken out into a separate repository.

The document processing parts are in the `topicdict_model` function.
This also enables tokenization / stemming / stopwords, etc. 
from `quanteda`.  The C++ part of the model only stores indexes
from which the words, topic assignments, etc. can be reconstructed.

There are two document sets bundled in `extdata`: `macavity` 
(small, for testing) and `bara_paras` (abortion debate, slightly 
preprocessed).

Post processing functions for a fitted model are all in `posterior.R`
With the exception of some otherwise innocent unused variables in the C++, 
the package passes CRAN checks.

Will Lowe with C++ code from Tomo and Shusei. Jan 2018.

## `TOT` model
* Wang, Xuerui and Andrew McCallum. 2006. ``Topics over Time: A Non-Markov Continuous-Time Model of Topical Trends.'' Proceedings of the 12th ACM SIGKDD international conference on Knowledge discovery and data mining. [here](https://people.cs.umass.edu/~mccallum/papers/tot-kdd06s.pdf)

## Weighting
* To use non-weighting model, modify `int use_weight = 1` to `int use_weight = 0` [here](https://github.com/Shusei-E/keyATM/blob/weighting_lda/src/keyATM.h#L26). The model is in `weighting_lda` branch.
