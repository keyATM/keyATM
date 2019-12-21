---
name: An issue template
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

## Problem summary (required):
This is an example to show how to report a bug.

## Environment information (required):
```
> sessionInfo()  # please run this in R and copy&paste the output
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
