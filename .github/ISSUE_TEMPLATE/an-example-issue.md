---
name: An Example Issue
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

## Problem summary (required):
This is an example to show how to report a bug.

## Environment information (required):
```
> sessionInfo()
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin18.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS/LAPACK: /usr/local/Cellar/openblas/0.3.6_1/lib/libopenblasp-r0.3.6.dylib

locale:
[1] C/UTF-8/C/C/C/C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] keyATM_0.1.0

loaded via a namespace (and not attached):
 [1] RcppEigen_0.3.3.5.0 xfun_0.8            tidyselect_0.2.5
 [4] remotes_2.1.0       rematch2_2.1.0      purrr_0.3.2
 [7] lattice_0.20-38     colorspace_1.4-1    vctrs_0.2.0
[10] testthat_2.2.1      htmltools_0.3.6     usethis_1.5.1
[13] yaml_2.2.0          rlang_0.4.0         pkgbuild_1.0.4
[16] pkgdown_1.3.0       pillar_1.4.2        glue_1.3.1
[19] withr_2.1.2         sessioninfo_1.1.1   lifecycle_0.1.0
[22] stringr_1.4.0       munsell_0.5.0       commonmark_1.7
[25] gtable_0.3.0        hashmap_0.2.2       devtools_2.1.0
[28] evaluate_0.14       codetools_0.2-16    memoise_1.1.0
[31] knitr_1.24          callr_3.3.1         ps_1.3.0
[34] curl_4.0            Rcpp_1.0.2          backports_1.1.5
[37] scales_1.0.0        desc_1.2.0          pkgload_1.0.2
[40] fs_1.3.1            ggplot2_3.2.1       digest_0.6.21
[43] stringi_1.4.3       processx_3.4.1      dplyr_0.8.3
[46] ggrepel_0.8.1       grid_3.6.1          rprojroot_1.3-2
[49] cli_1.1.0           tools_3.6.1         magrittr_1.5
[52] lazyeval_0.2.2      tibble_2.1.3        crayon_1.3.4
[55] tidyr_1.0.0         pkgconfig_2.0.3     zeallot_0.1.0
[58] MASS_7.3-51.4       Matrix_1.2-17       xml2_1.2.2
[61] prettyunits_1.0.2   rmarkdown_1.14      httr_1.4.1
[64] assertthat_0.2.1    roxygen2_6.1.1      rstudioapi_0.10
[67] R6_2.4.0            compiler_3.6.1
```

## Actual output (required):
```
Error: package or namespace load failed for 'keyATM' in dyn.load(file, DLLpath = DLLpath, ...):
 unable to load shared object '/usr/local/lib/R/3.6/site-library/keyATM/libs/keyATM.so':
  dlopen(/usr/local/lib/R/3.6/site-library/keyATM/libs/keyATM.so, 6): Symbol not found: __keyATM_keyATM_train_cov
  Referenced from: /usr/local/lib/R/3.6/site-library/keyATM/libs/keyATM.so
  Expected in: flat namespace
 in /usr/local/lib/R/3.6/site-library/keyATM/libs/keyATM.so
Error: loading failed
Execution halted
ERROR: loading failed
* removing '/usr/local/lib/R/3.6/site-library/keyATM'
Error: Command failed (1)
```

## Expected output (it helps us a lot):
Successful installation.

## Minimal and reproducible example with less than 50 lines (it helps us a lot):
```
devtools::install_github("keyATM/keyATM")
```
