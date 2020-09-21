---
layout: default
title: BINNING
has_children: false
parent: Configuration files
nav_order: 3
---

# **`BINNING`**
This dictionary is the opportunity to define the axis binning of the user's space.
The binning values are split into x and y axis definitions where the x-axis describes
the variable whose signal region is blinded. Note that 2D Alphabet can rebin
and reduce the ranges of the input axes for you. The binning in each axis cannot be
beyond the range of the input histogram (hopefully, this is obvious) but it can be a
subset of the input range. The number of bins are restricted to be equal to or less
than the input number of bins (hopefully, this is also obvious). Additionally, any newly
defined bin edges must line up with input bin edges (ie. there's no bin splitting).

The signal bounds only apply to the `X` axis and must exist within the `X` axis range
defined in the configuration file so that there are always three regions defined
with non-zero width. Only the signal region will be blinded during the fitting
and plotting. The options to perform this blinding are described in the [`OPTIONS` section](OPTIONS).

For either axis, there are two options for binning: constant width and variable width.
To use variable bins, provided a list of bin edges to the `BINS` key. To have 
2D Alphabet calculate the bin edges with constant bin width, use the `MIN`, `MAX`, and
`NBINS` keys (as you would in a `TH1` definition).



An example formatting of the `BINNING` section is provided here.

```json
"X":
    "NAME": <name of internal variable on the x-axis>, # string
    "TITLE": <title to appear in plots with this axis>, # string
    "MIN": <lower bound of x-axis>, # float
    "MAX": <upper bound of x-axis>, # float
    "NBINS": <number of x bins from MIN to MAX>, # int
    "SIGSTART": <lower bound of signal region of x-axis>, # int
    "SIGEND": <upper bound of signal region of x-axis> # int
"Y"`
    "NAME": name of your variable on the y-axis, # string
    "TITLE": title that you'd like to appear in plots with this axis, # string
    "BINS": [<edge1>, <edge2>,...<edgeN>] # [float]
```

Because the `X` axis is split into three portions (`LOW`, `SIG`, `HIGH`), one can also
define binning per region of the `X` axis. `BINS` can be used as well as `MIN`, `MAX`, `NBINS` like so:

```json
"X":
    "NAME": <name of internal variable on the x-axis>, # string
    "TITLE": <title to appear in plots with this axis>, # string
    "LOW": {
        "MIN": <lower bound of x-axis>, # float
        "MAX": <upper bound of x-axis>, # float
        "NBINS": <number of x bins from MIN to MAX> # int
    }
    "SIG": {
        "MIN": <lower bound of x-axis>, # float
        "MAX": <upper bound of x-axis>, # float
        "NBINS": <number of x bins from MIN to MAX> # int
    }
    "HIGH": {
        "MIN": <lower bound of x-axis>, # float
        "MAX": <upper bound of x-axis>, # float
        "NBINS": <number of x bins from MIN to MAX> # int
    }
    "SIGSTART": <lower bound of signal region of x-axis>, # int
    "SIGEND": <upper bound of signal region of x-axis> # int
```



