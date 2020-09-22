---
layout: default
title: Nuisance parameter impacts - run_Impacts.py
nav_order: 4
has_children: false
parent: How to run
---

# Nuisance parameter impacts - `run_Impacts.py`

This is the simplest of the four tools (it will probably be absorbed into
`run_Stats.py` in the future). It just wraps the commands to `combineTool.py`
that run the calculation of the nuisance parameter impacts on the value
of `r`. To produce the plot of the impacts, run the original command
with `--post` appended.

Running `python run_Impacts.py --help` returns the following:

```
Options:
  -h, --help            show this help message and exit
  -d <dir>, --projDir=<dir>
                        Home of the project - has the cards, fit results, etc
  --condor              Turn condor grid submission on
  -p, --post            Run the post processing to get impact plot
```