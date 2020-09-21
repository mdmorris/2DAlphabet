---
layout: default
title: How to run
nav_order: 3
has_children: true
---

# How To Run
The 2D Alphabet framework is designed to run by calling
a one of a set of python scripts all with the prefix
`run_`. The simplest and most important of these is
`run_MLfit.py`. 

The `run_MLfit.py`
script also creates the starting point for all other `run_`
scripts because it:
- Takes the model from pre-fit (where we know essentially nothing
  about the total background estimate) to post-fit where background
  estimate is meaningful,
- Creates a directory where all of the information of the run can
  be stored and accessed for later use.

Thus, `run_MLfit.py` should be run before the other scripts whenever
attempting to fit a new Alphabet pair or a new set of pairs.
