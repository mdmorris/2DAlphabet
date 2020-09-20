# 2DAlphabet
~~See Documentation folder for full instructions and explanations in main.pdf~~

While the documentation may be useful, it is woefully outdated. For install, please do the following:

```
cmsrel CMSSW_10_6_14
cd CMSSW_10_6_14/src
cmsenv
git clone https://github.com/lcorcodilos/2DAlphabet.git
git clone https://github.com/lcorcodilos/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit/
bash <(curl -s https://raw.githubusercontent.com/lcorcodilos/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh)
scram b -j 10
```

Note that the CMSSW release is different from the release recommended by the Combine tool documentation in order to maintain compatibility with Fermilab's LPC changing to SL7 by September 1st, 2020.