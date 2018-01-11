# topicdict

The topic dictionary code, broken out into a separate repository.
It's not clear to me whether I've faithfully rendered T&S's code.
They should check.  It's all in about 300 lines in `src/train.cpp`.
The user visible gibbs sampler function is `topicdict_train`.

The document processing parts are in the `topicdict_model` function.
This also enables tokenization / stemming / stopwords, etc. 
from quanteda.  The C++ part of the model only stores indexes
from which the words, topic assignments, etc. can be reconstructed.

There are two document sets bundled in `extdata`: `macavity` 
(small, for testing) and `bara_paras` (abortion debate, slightly 
preprocessed).  The vignette exercises Bara.

Post processing functions for a fitted model are all in `posterior.R`
With the exception of some unused variables in the C++, the package
passes CRAN checks.

The second vignette `data-structures` has some notes about the architecture
of the package and what might have gotten lost in translation. Maybe read this
before the help pages and code.

## Notes

* ~~On most machines it's been necessary to drop a recent boost into `/usr/local/include`~~

## Remaining:

* log likielihood is not yet right (nans and infs)

Will Lowe. Jan 2017
