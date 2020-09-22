---
layout: default
title: Configuration files
has_children: true
nav_order: 3
---

# Configuration files

The goal of the JSON is to have an easily understandable and configurable input
that allows the
2D Alphabet software to read and organize analysis files and histograms to its
liking while also giving the user the ability to easily configure values like
bin sizes and ranges without lots of command line options. This means that while
the user is encouraged to add what they need, there are some keys and values
that must stay the same. These static strings are always in capital letters
to make them easy to distinguish. The six static sections are described below.

## Groups of configuration files

Multiple configuration files can be provided to 2D Alphabet to simultaneously
fit pass-fail pairs. There are several important notes to make about
this feature.

* The data between the configuration files **must** be statistically independent.
  It is the responsibility of the user to ensure this in their selection.
* Any commonly named systematic uncertainties across configuration files will
  automatically be correlated in the fit. For example, if you have `syst1`
  defined and used in `config1.json` and `config2.json`, a single nuisance parameter
  will be correlated across the `config1` and `config2` pass-fail pairs.
* Configuration files with a commonly named systematic still use the uncertainty
  templates provided in each corresponding configuration file. The templates
  are just mapped to the same nuisance  parameter generated (and tied to the Gaussian constraint).
* Provide unique names to systematic uncertainties across configuration files
  if you wish for them to be uncorrelated.