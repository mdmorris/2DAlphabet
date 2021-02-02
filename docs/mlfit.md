---
layout: default
title: Likelihood fit - run_MLfit.py
nav_order: 1
has_children: false
parent: How to run
---

# Likelihood fit - `run_MLfit.py`
This tool takes as input a list of JSON configuration
files (detailed in [configs](configs) via a command such as

```bash
python run_MLfit.py input_config1.json input_config2.json ... input_configN.json
```

where `input_config*.json` are the configurations for each
Alphabet pair. Based on the input from the JSON files, 2D
Alphabet grabs all distributions, generates the RooFit objects
from them, creates the pass and fail distributions of the
QCD/non-resonant background estimated from data, passes these
all to Combine, calls Combine (`-M FitDiagnostics`), interprets the fit result, and
plots the post-fit distributions. It also outputs additional
plots and read-out for debugging or reference. The `run_MLfit.py`
script also creates the starting point for all other `run_`
scripts because it:
- Takes the model from pre-fit (where we know essentially nothing
  about the total background estimate) to post-fit where background
  estimate is meaningful,
- Creates a directory where all of the information of the run can
  be stored and accessed for later use.

Thus, `run_MLfit.py` should be run before the other scripts whenever
attempting to fit a new Alphabet pair or a new set of pairs.

The main command line options are provided below with longer descriptions
along with the output of `python run_MLfit.py --help`.

- **`-q, --tag`** Assigns a tag for the run. The tag determines the 
naming of the folder (and some of the nested objects) the will be created
during processing. This folder is often referred to as the "project directory".
This option takes precedent
over tags defined in configuration files and is a good way to ignore the need
for configuration files to have the same tag and to make
sure you don't overwrite an old result on-the-fly.
- **`--rMin/--rMax`** These are the minimum and maximum bounds of the
signal strength (r) in the fit. They default to 0 and 5, respectively.
It may be useful to loosen or tighten the bounds if the fit is failing
near one of the boundaries.
- **`--recycleAll`** Recycles everything generated in the previous run
using the given tag. This means the construction of the workspace, rebinned histograms, and
other steps previous to the fit are skipped and the previous run versions
are loaded instead. Use this if you haven't changed your configuration
file and would like to speed things up.
- **`--skipFit`** Skips running the fit and goes directly to plotting.
- **`--skipPlots`** Skips running the plots. Sometimes this is the part
that takes the longest because a sampling method is used to estimate errors.

```
Options:
  -h, --help            show this help message and exit
  -q TAG, --tag=TAG     Assigns a tag for this run
  -s V1=1.0,V2=1.0..., --setParameter=V1=1.0,V2=1.0...
                        String of parameters to set pre-fit. Uses same comma
                        separated format as Combine (V1=1.0,V2=1.0...)
  --rMin=rMin           Minimum bound on r (signal strength)
  --rMax=rMax           Minimum bound on r (signal strength)
  --recycleAll          Recycle everything from the previous run with this tag
  --skipFit             Skip fit and go directly to plotting (WARNING: Will
                        use previous fit result if it exists and crash
                        otherwise)
  --skipPlots           Skip plotting
  --fullRun2            Plot sum of years 16, 17, 18
  --CL=CL               Command-line options to set for all configs
```