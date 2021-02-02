---
layout: default
title: Installation
nav_order: 2
---

# Installation

2D Alphabet can only be used in a CMSSW environment. The Higgs Analysis
Combine Tool must be installed. Please follow the instructions below to
checkout a 2D Alphabet-friendly version of Combine. These instructions are
based off of the [Combine documentation](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) 
for 102X setup. Please cross-check the instructions here with the official
instructions. The only difference should be in the cloned repository and lack
of branch change. Note also that the CMSSW release is different from the
release recommended by the Combine tool documentation in order to maintain
compatibility with Fermilab's LPC changing to SL7 by September 1st, 2020.

For `csh`:
```sh
    set SCRAM_ARCH=slc7_amd64_gcc700
    cmsrel CMSSW_10_2_13
    cd CMSSW_10_2_13/src
    cmsenv
    git clone https://github.com/lcorcodilos/2DAlphabet.git
    git clone https://github.com/lcorcodilos/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    bash <(curl -s https://raw.githubusercontent.com/lcorcodilos/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh)
    scram b clean; scram b -j 10
    cmsenv
```

For `bash`:
```bash
    export SCRAM_ARCH=slc7_amd64_gcc700
    cmsrel CMSSW_10_6_14
    cd CMSSW_10_6_14/src
    cmsenv
    git clone https://github.com/lcorcodilos/2DAlphabet.git
    git clone https://github.com/lcorcodilos/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit/
    curl -s https://raw.githubusercontent.com/lcorcodilos/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh | bash
    scram b clean; scram b -j 10
    cmsenv
```

Run `combine --help` and check it returns the help menu to confirm you've successfully setup combine.

Try to call RooParametricHist2D in interactive python if you're feeling
uneasy and to ensure everything is working. 

```python
import ROOT
r = RooParametricHist()
```
