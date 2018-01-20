# topicdict

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

## Remaining:

* ~~log likielihood is not yet right (nans and infs)~~
* slice sampler not working yet

Will Lowe. Jan 2018
