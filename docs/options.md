---
layout: default
title: OPTIONS
has_children: false
parent: Configuration files
nav_order: 6
---

# `OPTIONS`

This section is dedicated to providing per-config options to 2D Alphabet. 
Unless noted, these options only affect `run_MLfit.py` since
that tool is responsible for building the model performing
the initial fit that is read by the other `run_*.py` tools.

Required among these options the `name` for the configuration file
so that the pass-fail pair can be uniquely identified among those provided
to 2D Alphabet.
Truly optional options include those for blinding, plotting, and changing 
fit functionality. The most common possible options are described in detail below.

- **`name`** (`string` - default=`None`): Required. Provides a unique name to the configuration file.
- **`tag`** (`string` - default=`None`): Provides a tag to group configuration files and 
  create an output directory of the same name ("project directory"). Overriden by
  the `-q/--tag` command-line option to `run_MLfit.py` which will set the provided tag
  to all configuration files.
- **`year`** (`int` - default=1): Default `1` indicates this is full Run II. Other options
  are `16`, `17`, and `18`. This option determines what luminosity will be included in the plots
  and helps group configuration files when plotting the full Run II together (when years have been fit
  separately).
- **`blindedFit`** (`bool` - default=`true`): Blinds the fit to the signal region defined by `SIGSTART` and `SIGEND`.
- **`blindedPlots`** (`bool` - default=`true`): Blind the plots to the signal region defined by `SIGSTART` and `SIGEND`.
- **`plotUncerts`** (`bool` - default=`true`): Turns on plotting of the `X` and `Y` projections of the
  shape uncertainty templates provided (after rebinning). Useful for comparing nominal templates
  with uncertainty variations.
- **`prerun`** (`bool` - default=`true`): Runs the ML fit for this configuration file
  by itself and feeds the resulting transfer function parameters into the simultaneous
  fit with the other configuration files. This option can be effective if the fit
  cannot converge to a minimum because of too few steps.
- **`plotPrefitSigInFitB`** (`bool` - default=`false`): By default, the background-only
  fit result plots will include no signal since `r = 0`. This option set to `true` will
  plot the pre-fit signal (normalized to `r = 1`) instead. The s+b fit results always
  plot with the post-fit value of `r`.
- **`plotTitles`** (`bool` - default=`true`): By default, plot titles are always included in the plotting
  (this is useful for determining the boundaries of slices that are projected onto the opposite axis).
  Set to `false` this option turns off plot titles so that they are more appropriate for publication.
- **`addSignals`** (`bool` - default=`true`): By default, this option adds the signals for the sake of plotting so that
  they are considered as "one signal" in the legend and plots. This is useful if you have defined a separate signal for each
  year in a Run II configuration file. When set to `false`, the plotting will consider the signals
  as separate and they will be drawn separately and with different legend entries. This is useful
  if you have multiple unique signals.

The following are "developer/advanced options" that would mainly be useful for anyone
working on the 2D Alphabet code or debugging a fit.

- **`draw`** (`bool -default=false`): This will live draw plots as they are created
  (sets `gROOT.SetBatch(kTRUE)`).
  Useful if one has pauses in the code to view the plots. Since 2D Alphabet most likely
  is running remotely, this will be very slow!
- **`verbosity`** (`int - default=0`): Sets Combine verbosity. Setting in the first of
  the provided configuration files will determine behavior of the group run.
- **`fitGuesses`** (dictionary - default=None): Not recommended for use since it has gone
  untested for quite some time. This is because it is largely unnecessary and will
  most likely be removed. For the sake of documentation, the option creates slices
  along `Y` and projects them onto `X` to perform 1D fits in the slices. The parameters
  of the 1D fits are then fit as a function of `Y` to determine the `Y` dependence.
  The output of this pseudo-2D fit is then fed into the true 2D fit with Combine as
  a better pre-fit starting point.
- **`overwrite`** (`bool` - default=`false`): Explicitly deletes the configuration file's
  corresponding sub-folder in the project directory before running. This is usually 
  not necessary unless you've done something such as remove a MC-based background
  from the model without changing the tag. The plots of that process from the previous
  run would remain since no new plots would be made to replace them.
- **`recycle`** - (`[string]` - default=[]): Not recommended for use since the introduction
  of the `--recycleAll` option of `run_MLfit.py`. Recycles named pieces of saveOut.p. 
  Available pieces are "name", "tag", "xVarName", "yVarName", "xVarTitle", "yVarTitle",
  "sigStart", "sigEnd", "freezeFail", "blindedFit", "blindedPlots", "newXbins", "full_x_bins",
  "newYbins", "rpf", "rpfVarNames", "organizedDict", "floatingBins". I cannot think
  of a good use for this option that would not also be dangerous or overcome by just re-running!
        