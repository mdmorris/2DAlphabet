---
layout: default
title: PROCESS
has_children: false
parent: Configuration files
nav_order: 1
---

# **`PROCESS`**

In this section, the user can define as many processes as they need. This 
includes data, background simulation, and signal simulation. Please note two 
important things: 
1. You should NOT define the non-resonant background that is meant to be 
   estimated from data ("multijet"). That background is naturally defined by
   the difference between data and the other backgrounds defined here.
2. Combine always requires that there be an observation (data), a background
   estimate, and a signal sample. This means that your configuration file
   must contain data (code 1) and signal (code 0). The codes classify the
   processes and are defined below. 

Each key in `PROCESS` is the name of each process of interest. Please name
your data as `data_obs`. This is a Combine convention that 2D Alphabet maintains.
Each process name is a key for a sub-dictionary that specifies
- `"TITLE":`
- `"FILE":` the path to the file containing the nominal pass and fail histograms for the process (string);
- `"HISTPASS":` the name of the histogram in the above file with the passing distribution (string);
- `"HISTFAIL":` the name of the histogram in the above file with the failing distribution (string);
- `"SYSTEMATICS":` a list of strings that correspond to the names of systematics
  in `SYSTEMATIC` (note `SYSTEMATIC` not `SYSTEMATICS`) that are applicable 
  to this process (list of strings);
- `"CODE":` a way to classify the treatment of the process: 0 (signal), 1 (data), 2 (background simulation)

Here is an example PROCESS dictionary section:

```json
"PROCESS": {
    "data_obs": {
        "FILE": "data_tau32medium_default.root",
        "HISTPASS": "MtwvMtPass",
        "HISTFAIL": "MtwvMtFail",
        "SYSTEMATICS":[],
        "CODE": 1
    },
    "ttbar": {
        "TITLE":"t#bar{t}",
        "FILE":"ttbar_tau32medium_default.root",
        "HISTPASS":"MtwvMtPass",
        "HISTFAIL":"MtwvMtFail",
        "SYSTEMATICS":["syst1","syst2","syst3","syst4"],
        "COLOR": 2,
        "CODE": 2
    },
    "signalLH2400": {
        "TITLE":"b*_{LH} (2.4 TeV)",
        "FILE": "signalLH2400_tau32medium_default.root",
        "HISTPASS": "MtwvMtPass",
        "HISTFAIL": "MtwvMtFail",
        "SYSTEMATICS":["syst1","syst3","syst5","jmr17","jms17"],
        "CODE": 0,
        "COLOR": 0
    }
```