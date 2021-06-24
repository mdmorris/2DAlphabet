Configuration files
===================
.. contents::
    :local:

The goal of the JSON is to have an easily understandable and configurable input
that allows the
2D Alphabet software to read and organize analysis files and histograms to its
liking while also giving the user the ability to easily configure values like
bin sizes and ranges without lots of command line options. This means that while
the user is encouraged to add what they need, there are some keys and values
that must stay the same. These static strings are always in capital letters
to make them easy to distinguish. The six static sections are described below.

Note on groups of configuration files
-------------------------------------
Multiple configuration files can be provided to 2D Alphabet to separately configure
regions of space but simultaneously them. There are several important notes to make about
this feature.

* The data between the configuration files **must** be statistically independent.
  It is the responsibility of the user to ensure this in their selection.
* Any commonly named systematic uncertainties across configuration files will
  automatically be correlated in the fit. For example, if you have ``syst1``
  defined and used in ``config1.json`` and ``config2.json``, a single nuisance parameter
  will be correlated across the ``config1`` and ``config2``.
* Configuration files with a commonly named systematic still use the uncertainty
  templates provided in each corresponding configuration file. The templates
  are just mapped to the same nuisance parameter (and tied to the same Gaussian constraint).
* Provide unique names to systematic uncertainties across configuration files
  if you wish for them to be uncorrelated.

Section: ``GLOBAL``
-------------------
This section is designed to help users with large configuration file
names by allowing them to create JSON-wide variables. For example,
if all of your files are located in ``/long/path/to/my/files/``, you 
can store this string in the GLOBAL dictionary with a custom key 
(let's say ``path``). Now instead of having to write the full directory
path for every process and systematic, the user can just write ``path``.
This simplifies the JSON and also has the standard advantages of using
variables over several instances of the same object.

This functionality works by searching all strings in the JSON for instances
of each key in ``GLOBAL`` and replacing the key with its corresponding dictionary value.
Note that it does not work for booleans or integers/floats.

The user must be careful they don't accidentally use strings in the JSON
that are identical to keys in ``GLOBAL`` or strange substitutions could happen.
This means keys in ``GLOBAL`` should be at least partially descriptive 
(single character keys would be a very bad idea). 

Section: ``OPTIONS``
--------------------
This section is dedicated to providing per-config options to 2DAlphabet. 
Unless noted, these options only affect ``run_MLfit.py`` since
that tool is responsible for building the model performing
the initial fit that is read by the other ``run_*.py`` tools.

Required among these options is the ``name`` for the configuration file
so that the config can be uniquely identified among those provided
to 2DAlphabet.
Truly optional options include those for blinding, plotting, and changing 
fit functionality. The most common possible options are described in detail below.

* **name** (``string``, default=``None``): Required. Provides a unique name to the configuration file.
* **tag** (``string``, default=``None``): Provides a tag to group configuration files and 
  create an output directory of the same name ("project directory"). Overriden by
  the `-q/--tag` command-line option to `run_MLfit.py` which will set the provided tag
  to all configuration files.
* **year** (``int``, default=1): Default `1` indicates this is full Run II. Other options
  are `16`, `17`, and `18`. This option determines what luminosity will be included in the plots
  and helps group configuration files when plotting the full Run II together (when years have been fit
  separately).
* **blindedFit** (``bool``, default=``true``): Blinds the fit to the signal region defined by ``SIGSTART`` and ``SIGEND``.
* **blindedPlots** (``bool``, default=``true``): Blind the plots to the signal region defined by ``SIGSTART`` and ``SIGEND``.
* **plotUncerts** (``bool``, default=``true``): Turns on plotting of the ``X`` and ``Y`` projections of the
  shape uncertainty templates provided (after rebinning). Useful for comparing nominal templates
  with uncertainty variations.
* **prerun** (``bool``, default=``true``): Runs the ML fit for this configuration file
  by itself and feeds the resulting transfer function parameters into the simultaneous
  fit with the other configuration files. This option can be effective if the fit
  cannot converge to a minimum because of too few steps.
* **plotPrefitSigInFitB** (``bool``, default=``false``): By default, the background-only
  fit result plots will include no signal since `r = 0`. This option set to ``true`` will
  plot the pre-fit signal (normalized to ``r = 1``) instead. The s+b fit results always
  plot with the post-fit value of ``r``.
* **plotTitles** (``bool``, default=``true``): By default, plot titles are always included in the plotting
  (this is useful for determining the boundaries of slices that are projected onto the opposite axis).
  Set to ``false`` this option turns off plot titles so that they are more appropriate for publication.
* **addSignals** (``bool``, default=``true``): By default, this option adds the signals for the sake of plotting so that
  they are considered as "one signal" in the legend and plots. This is useful if you have defined a separate signal for each
  year in a Run II configuration file. When set to ``false``, the plotting will consider the signals
  as separate and they will be drawn separately and with different legend entries. This is useful
  if you have multiple unique signals.

The following are "developer/advanced options" that would mainly be useful for anyone
working on the 2D Alphabet code or debugging a fit.

* **draw** (``bool``, default=``false``): This will live draw plots as they are created
  (sets `gROOT.SetBatch(kTRUE)`).
  Useful if one has pauses in the code to view the plots. Since 2D Alphabet most likely
  is running remotely, this will be very slow!
* **verbosity** (``int``, default=``0``): Sets Combine verbosity. Setting in the first of
  the provided configuration files will determine behavior of the group run.
* **fitGuesses** (``dict``, default=``None``): Not recommended for use since it has gone
  untested for quite some time. This is because it is largely unnecessary and will
  most likely be removed. For the sake of documentation, the option creates slices
  along ``Y`` and projects them onto ``X`` to perform 1D fits in the slices. The parameters
  of the 1D fits are then fit as a function of ``Y`` to determine the ``Y`` dependence.
  The output of this pseudo-2D fit is then fed into the true 2D fit with Combine as
  a better pre-fit starting point.
* **overwrite** (``bool``, default=``false``): Explicitly deletes the configuration file's
  corresponding sub-folder in the project directory before running. This is usually 
  not necessary unless you've done something such as remove a MC-based background
  from the model without changing the tag. The plots of that process from the previous
  run would remain since no new plots would be made to replace them.
* **recycle** - (``[string]``, default=``[]``): Not recommended for use since the introduction
  of the ``--recycleAll`` option of ``run_MLfit.py``. Recycles named pieces of saveOut.p. 
  Available pieces are "name", "tag", "xVarName", "yVarName", "xVarTitle", "yVarTitle",
  "sigStart", "sigEnd", "freezeFail", "blindedFit", "blindedPlots", "newXbins", "full_x_bins",
  "newYbins", "rpf", "rpfVarNames", "organizedDict", "floatingBins". I cannot think
  of a good use for this option that would not also be dangerous or overcome by just re-running!
        
Section: ``PROCESS``
--------------------
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

Each key in ``PROCESS`` is the name of each process of interest. Please name
your data as ``data_obs``. This is a Combine convention that 2D Alphabet maintains.
Each process name is a key for a sub-dictionary that specifies
- ``"TITLE"``
- ``"FILE"`` the path to the file containing the nominal pass and fail histograms for the process (string);
- ``"HISTPASS"`` the name of the histogram in the above file with the passing distribution (string);
- ``"HISTFAIL"`` the name of the histogram in the above file with the failing distribution (string);
- ``"SYSTEMATICS"`` a list of strings that correspond to the names of systematics
in ``SYSTEMATIC`` (note ``SYSTEMATIC`` not ``SYSTEMATICS``) that are applicable 
to this process (list of strings);
- ``"CODE"`` a way to classify the treatment of the process: 0 (signal), 1 (data), 2 (background simulation)

Here is an example `PROCESS` dictionary section:

.. code-block:: json

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

Section: ``SYSTEMATIC``
-------------------------
Because it bears repeating, please note the difference between this section, ``SYSTEMATIC``,
and the list of ``SYSTEMATICS`` defined inside the ``PROCESS`` dictionary. The ``SYSTEMATIC``
dictionary is a place to define as many systematics as a user may need. Similar to the
processes, each key in ``SYSTEMATIC`` is the name of the systematic in the analysis and each
is classified by a code that determines how the systematic will be treated. However, the
dictionary the user defines for a given systematic is different depending on what type it is.
The self-explanatory codes are 1 and 2 which are log-normal uncertainties on the normalization.
Less obvious are codes 2 and 3 which are for shape based uncertainties (and thus have
corresponding histograms) and are either in the same file as the process's nominal histogram
(code 2) or in a separate file (code 3). Additionally, they have a scale value which allows
the user to change the Gaussian constraint on the shape. For no change in the constraint, use 1.0.
If you have templates representing a 2 :math:`\sigma` shift, use 0.5 to properly constrain
the associated nuisance parameter during the shape interpolation with Combine.

Symmetric, log-normal
^^^^^^^^^^^^^^^^^^^^^
.. code-block:: json

    {
      "CODE": 0,
      "VAL": <uncertainty> # float
    }

Asymmetric, log-normal
^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: json

    {
        "CODE": 1,
        "VALUP": <+1 sigma uncertainty>, # float
        "VALDOWN": <-1 sigma uncertainty> # float
    }

Shape based uncertainty - in same file as nominal histogram
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: json

    {
        "CODE: 2",
        "HISTPASS_UP": <name of hist (in same file as nominal hist) for +1 sig uncertainty in pass distribution>, # string
        "HISTPASS_DOWN": <name of hist (in same file as nominal hist) for -1 sig uncertainty in pass distribution>, # string
        "HISTFAIL_UP": <name of hist (in same file as nominal hist) for +1 sig uncertainty in fail distribution>, # string
        "HISTFAIL_DOWN": <name of the hist (in same file as nominal hist) for -1 sig uncertainty in fail distribution>, # string
        "SCALE": <scale value to change scale of nuisance constraint> # float
    }

Shape based uncertainty, in different file as nominal histogram.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is the more flexible but also more complicated option. The user can specify files three different ways. 

1. By using ``FILEUP`` and ``FILEDOWN`` to pick a file that *every* process can pull the shape
   uncertainty histograms from. 
2. Use keys of the form ``FILEUP_<proc>`` 
   where ``<proc>`` matches the name of a process that is defined in the ``PROCESS`` dictionary and
   has this shape uncertainty associated with it. This allows each systematic and process to come
   from a separate file. 
3. Use keys of the form ``FILEUP_*`` where the * acts
   as a wild card for the process and must also exist in the file name where the process would
   normally be written. For example, if a :math:`t\bar{t}` distribution with +1 :math:`\sigma` pileup
   uncertainty is stored in ``ttbar_pileup_up.root`` and the corresponding signal distributions
   are in ``signal_pileup_up.root``, one can use the key-value pair ``"FILE_UP_*":"*_pileup_up.root"``.

The user can also specify histogram names in four different ways.

1. ``HISTPASS`` and ``HISTFAIL`` which allows the user to specify only two histogram
   names if they don't change between "up" and "down" shapes. 
2. The second is if the "up" and "down" shapes *do* have different histogram names
   and uses the form ``HISTPASS_UP`` ``HISTFAIL_UP``. 
3. The totally generic way allows the user to use the form ``HISTPASS_UP_<proc>`` where ``<proc>``
   matches the name of a process that is defined in the ``PROCESS`` dictionary and has this
   shape uncertainty associated with it. 
4. The "*" wildcard can be used in place of ``<proc>`` just as with the file keys.
   
Below is an example of the totally generic way (3). The various options lead to flexibility but
the more organized you are, the easier it is to write the configuration file!

.. code-block:: json

    {
        "CODE": 3,
        "FILEUP_<proc>": </path/to/fileup_<proc>_up.root>, # string
        "FILEDOWN_<proc>": </path/to/filedown_<proc>_down.root>, # string
        "HISTPASS_UP": <name of hist (in same file as nominal hist) for +1 sig uncertainty in pass distribution>, # string
        "HISTPASS_DOWN": <name of hist (in same file as nominal hist) for -1 sig uncertainty in pass distribution>, # string
        "HISTFAIL_UP": <name of hist (in same file as nominal hist) for +1 sig uncertainty in fail distribution>, # string
        "HISTFAIL_DOWN": <name of the hist (in same file as nominal hist) for -1 sig uncertainty in fail distribution>, # string
        "SCALE": <scale value to change scale of nuisance constraint> # float
    }

Section: ``BINNING``
--------------------
This dictionary is the opportunity to define the axis binning of the user's space.
The binning values are split into x and y axis definitions where the x-axis describes
the variable whose signal region is blinded. Note that 2D Alphabet can rebin
and reduce the ranges of the input axes for you. The binning in each axis cannot be
beyond the range of the input histogram (hopefully, this is obvious) but it can be a
subset of the input range. The number of bins are restricted to be equal to or less
than the input number of bins (hopefully, this is also obvious). Additionally, any newly
defined bin edges must line up with input bin edges (ie. there's no bin splitting).

The signal bounds only apply to the ``X`` axis and must exist within the ``X`` axis range
defined in the configuration file so that there are always three regions defined
with non-zero width. Only the signal region will be blinded during the fitting
and plotting. The options to perform this blinding are described in the `OPTIONS` section.

For either axis, there are two options for binning: constant width and variable width.
To use variable bins, provided a list of bin edges to the `BINS` key. To have 
2D Alphabet calculate the bin edges with constant bin width, use the ``MIN``, ``MAX``, and
``NBINS`` keys (as you would in a ``TH1`` definition).

An example formatting of the ``BINNING`` section is provided here.

.. code-block:: json

    "X":
        "NAME": <name of internal variable on the x-axis>, # string
        "TITLE": <title to appear in plots with this axis>, # string
        "MIN": <lower bound of x-axis>, # float
        "MAX": <upper bound of x-axis>, # float
        "NBINS": <number of x bins from MIN to MAX>, # int
        "SIGSTART": <lower bound of signal region of x-axis>, # int
        "SIGEND": <upper bound of signal region of x-axis> # int
    "Y":
        "NAME": name of your variable on the y-axis, # string
        "TITLE": title that you'd like to appear in plots with this axis, # string
        "BINS": [<edge1>, <edge2>,...<edgeN>] # [float]

Because the ``X`` axis is split into three portions (``LOW``, ``SIG``, ``HIGH``), one can also
define binning per region of the ``X`` axis. ``BINS`` can be used as well as ``MIN``, ``MAX``, ``NBINS`` like so:

.. code-block:: json

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

Section: ``FIT``
----------------
This section defines the values of the fit parameters for the analytic 
portion of the transfer function. The
2D fit can accommodate any functional form. Each parameter in equation
should be marked with a @ symbol (RooFit convention) starting at 0. Each
parameter should also be specified with a range and (if desired) an error 
(which does not establish a constraint but just an initial step size 
in the minimization). 

An example section which uses a polynomial in both directions (of 
different orders) is provided below. In this example, the range of each parameter
is -10 to 10 with a starting (``NOMINAL``) value of 1.0 and step size
(``ERROR``) of 0.1.

To accommodate for the potential of a negative transfer function value
in all or part of the 2D space (which would be unphysical and make 
Combine very mad), 2DAlphabet autoomatically wraps the provided function
in ``max(1e-6, <your func>)``.
Thus, one should always check in the plots that the final parameters
do not produce a ~0 transfer function in large portions of the space.
The fit result will not be stable. This can be addressed by reducing
the ranges (``MIN``/``MAX``) of the parameters or changing the functional 
form.

.. code-block:: json

    "FIT": {
        "FORM":"(@0+@1*x**1+@2*x**2)*(1+@3*y)",
        "0": {
            "NOMINAL": 1.0,
            "MIN":-10.0,
            "MAX":10.0,
            "ERROR":0.1
        },
        "1": {
            "NOMINAL": 1.0,
            "MIN":-10.0,
            "MAX":10.0,
            "ERROR":0.1
        },
        "2": {
            "NOMINAL": 1.0,
            "MIN":-10.0,
            "MAX":10.0,
            "ERROR":0.1
        },
        "3": {
            "NOMINAL": 1.0,
            "MIN":-10.0,
            "MAX":10.0,
            "ERROR":0.1
        }
    }