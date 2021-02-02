---
layout: default
title: Asymptotic Limit - run_Limit.py
nav_order: 2
has_children: false
parent: How to run
---

# Asymptotic Limit - `run_Limit.py`

This tool is for calculating limits on the signal strength (`r`) of
a given simulated signal. The script takes the background only fit result
`fit_b` output from `run_MLfit.py`, morphs a copy of the pre-fit workspace
to this result (thus bypassing having to perform a likelihood fit again),
generates toys from the background estimate of the morphed pre-fit workspace
to create pseudo-data, and fits these toys to derive expected limits for the
simulated signal (blinding can be turned off so that the real data is also
fit to get an observed limit on `r`). The actual pseudo-data generation
and asymptotic limit calculation is handled by the AsymptoticLimit method
of Combine.

Because the background only fit result is used, the signal simulation used
in the configuration file for the `run_MLfit.py` does not matter. However,
the signal for limit setting obviously matters and in all likelihood
(no pun intended), the user will want to scan over several simulated signals to
generate limits as a function of their samples (perhaps each one was generated
at a different mass). Note that the model is rebuilt each time a signal 
(and the associated uncertainty templates) needs to be changed.

You do NOT have to write a configuration file for each simulated sample. At
any point in the python call to `run_Limit.py`, one can provide a string
swap specified by the syntax `oldString:newString` which will replace all
instances of `oldString` in the configuration files provided with `newString`.
In fact this can be used in `run_MLfit.py` as well.

Running `python run_Limit.py --help` returns the following:

```
Options:
  -h, --help            show this help message and exit
  -q <tag>, --tag=<tag>
                        Assigns a tag for this run
  -d <dir>, --projDir=<dir>
                        Points to the directory where the b-only fit result is
                        located
  --unblindData         Unblind the observation and calculate the observed
                        limit
  --recycleAll          Recycle everything from the previous run with this
                        tag. Note that this does not allow for string
                        substitution.
```