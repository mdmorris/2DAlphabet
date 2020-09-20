---
layout: default
title: Home
nav_order: 1
---

# Introduction

2D Alphabet is a wrapper to construct the workspace for a specific type of
background estimate, provide input to the Combine Statistical Tool, plot the
2D distributions from the fit result, and provide the infrastructure to test
this result. 

The name of the wrapper is derived from the data driven background estimate
performed to estimate combinatorial backgrounds that are otherwise poorly
modeled by Monte Carlo simulation. In many cases, the background being
modeled is QCD multijet production and as a default, "multijet" is the name of the output
background produced in plots and "qcd" in the workspace. However, depending on the selection,
there may be other backgrounds accounted for as well such as V+jets. 

As one might expect, the data driven background estimation method is a two
dimensional version of the Alphabet method. The 1D Alphabet method works on
the concept of a transfer function for non-resonant backgrounds as a function
of a variable, $$x$$, which is not the search variable. By defining a discriminator, $$D$$,
one can separate data events
based on whether they pass or fail the selection on $$D$$ (the sets $$D_{pass}$$
and $$D_{fail}$$). Then resonant backgrounds from MC can be subtracted from these
selections to get estimates of the non-resonant background. Since the
non-resonant background is smooth in $$x$$ in both the pass and fail, the
pass-to-fail ratio should also be smooth (and positive). 

Using this criteria, one can calculate the ratio per-bin as a function
in background enriched upper and lower sidebands of
$$x$$ and then fit for a function, 
$$R_{P/F}(x)$$, that interpolates the $$x$$ signal region. Then those events that
fail $$D$$ and exist in the $$x$$ signal selection can be multiplied by
$$R_{P/F}(x)$$ to estimate $$D_{pass}$$ in the signal region.

This 1D method has been used successfully in the past. However, if the signal
region is wide, the
background estimate is susceptible to shifts that do not accurately describe
the behavior of the background being estimated. Additionally, if the 
$$R_{P/F}$$ has a strong dependence on the measurement variable (for example,
resonance mass), 1D Alphabet will not be able to accurately predict the
background in this variable. The 2D Alphabet method attempts to remedy both
of these issues by constraining the shape of the $$R_{P/F}$$ with the
measurement variable, $$y$$, as the second dimension. 

While fitting for this surface, it's convenient to also fit for other
background contributions. For example, say an analysis has as its main
backgrounds QCD and $$Bkg_A$$. $$Bkg_A$$ is modeled well enough by simulation that
the MC can be used directly. However, there may be a recommendation to fit
for some MC-to-data corrections given a systematic uncertainty (for example,
this is the recommendation for top $$p_{T}$$ reweighting for $$t\bar{t}$$ MC).
This means that the contribution of $$Bkg_A$$ can change and therefore, the
amount subtracted from data to form the initial QCD estimate in the $$x$$
sidebands can change! Thus, QCD and $$Bkg_A$$ are fit simultaneously to derive
a total background estimate. The details of accounting for systematic
uncertainties are explained in the [systematics section](systematics).

A useful feature of the 2D Alphabet setup is that multiple pairs of pass and
fail distributions (called "Alphabet pairs") can be fit simultaneously. For
example, one can setup this estimate for 2016, 2017, and 2018 data and
selections (three pairs of pass and fail) and then fit all three at the same
time with nuisance parameters easily correlated across pairs. 

The statistical tool of choice to perform the fitting is the Higgs Analysis
Combine Tool. The documentation for the tool can be found
[here](https://cms-hcomb.gitbooks.io/combine/content/).
One significant addition is made to the central Combine release to make
Combine 2D Alphabet-friendly - the RooParametricHist2D class. This is
identical to the RooParametricHist class already provided by Combine but take
as input a TH2 instead of a TH1. The accompanying changes to accommodate this
class are made in the Combine Tool code. Additionally, code has been added
and modified in Combine Harvester (used for plotting) to plot the 2D
distributions. User's not interested in development do not need to worry
about the specifics of this since calls to RooParametricHist2D are all
internal. 