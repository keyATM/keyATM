# topicdict

The topic dictionary code, broken out into a separate repository.

(no longer necessary: On most machines it's been necessary to drop a recent boost into `/usr/local/include`)

## BUGS

* Not clear from the original code what num_topics means (must it really be the same number as the 
  seed topics?) crash is due to out of scope indexing in `train` at line 152.  

## TODO

* Figure out how T&S's `phi_s` code was meant to work.
* Complete the X updates
* Put the log likelihood calculation back in
* Add alpha back in
* Check for user stops in the code
* Confirm what's happening in the Z update for X==1

Will Lowe. Jan 2017
