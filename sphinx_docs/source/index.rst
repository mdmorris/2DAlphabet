.. 2DAlphabet documentation master file, created by
   sphinx-quickstart on Wed May 26 17:13:06 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to 2DAlphabet's documentation!
======================================
2D Alphabet is a framework to construct the workspace for a specific type of
background estimate, provide input to the Combine statistical tool, plot the
2D distributions from the fit result, and provide the infrastructure to test
this result. 

The documentation for the Higgs Analysis Combine Tool can be found
`here <https://cms-hcomb.gitbooks.io/combine/content/>`_.
One significant addition is made to the central Combine release to make
Combine 2D Alphabet-friendly - the RooParametricHist2D class. This is
identical to the RooParametricHist class already provided by Combine but take
as input a TH2 instead of a TH1. The accompanying changes to accommodate this
class are made in the Combine Tool code. Additionally, code has been added
and modified in Combine Harvester (used for plotting) to plot the 2D
distributions. User's not interested in development do not need to worry
about the specifics of this since calls to RooParametricHist2D are all
internal.

Overview of model building approach
------------------------------------
From Scratch: Multijet background estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The name of the framework is derived from its data-driven background estimate
of combinatorial backgrounds that are otherwise poorly
modeled by Monte Carlo simulation. In many cases, the background being
modeled is QCD multijet production and as a default, "multijet" is the name of the output
background produced in plots and "qcd" is the name in the workspace. However, depending on the selection,
there may be other backgrounds accounted for as well such as V+jets. 

As one might expect, the data driven background estimation method is a two
dimensional version of the Alphabet method which itself is named after the
ABCD method. All three use an analytic "transfer function" to "transfer" the
background contribution in a control region to the contribution in the signal
region. If the shapes of the background distributions in the control region and
signal region are identical, then the transfer function is just a constant factor
which just changes the normalization from one region to the other.

However, having a "flat" transfer function is unusual as there is typically
a shape dependence along the measurement variable. The ABCD method measures
data distributions in selection regions A, B, and D (see figure below) which
are enriched in background and depleted of signal. In the figure,
``var1`` and ``var2`` are the selection variables and the C region is the
signal region of the analysis. Binned distributions are created in some third variable
(``var3``, the measurement variable) for each of the four regions. The ratio
of A/B and C/D are assumed to be equal and therefore, :math:`A/B * D = C`.
Here, :math:`A/B` is the transfer function.

.. code-block:: none
   :caption: Figure: ABCD space
   :name: ABCD

   var1 |          |
        |    A     |     C
        |          |
        |----------------------
        |          |
        |    B     |     D
        |__________|___________ 
                              var2

Technically speaking, to create your QCD estimate, you measure the binned (or analytic depending what you are doing)
transfer function from A and B and then weight events in region D to get the estimate along your ``var3`` axis
of the QCD.

While this method has proven successful in the past, it has the disadvantage that
it extrapolates the shape of the background to region C by assuming that A/B and C/D
are equal. If a different var2 can be chosen such that the signal lies in the middle
of the axis (still region C of the next figure), one can instead interpolate the
background which is often more robust than extrapolation. It would be inconvenient
to call this the ABCDEF method so it is instead referred to
as "Alphabet" since it uses so many letters.

.. code-block:: none
   :caption: Figure: Alphabet space
   :name: alphabet

   var1 |          |           |
        |    A     |     C     |    E
        |          |           |
        |---------------------------------
        |          |           |
        |    B     |     D     |    F
        |__________|___________|__________
                                       var2

Of course, now one can't use the simple :math:`A/B * D = C` equality. Instead,
the ``var2`` axis is itself binned and two distributions along ``var2`` are created
- the upper and bottom portions of ``var1``.
One can visualize it as the following (where the `.` markers
are representing some histogram counts - D and C are empty since those are blinded).

.. code-block:: none
   :caption: Figure: Split Alphabet regions
   :name: alphabetsplit

   this is a test

    var1  |         |       |
   (upper)| .   A   |   C   |    E
          |_:_:_:_:_|_______|_._._._.__
                                     var2
     
    var1  |         |       |
   (lower)|    B    |   D   | .   F
          |_._:_:_._|_______|_:_:_:_:___
                                     var2

The ratio of the upper and lower histograms is calculated and you get some ratio-per-bin
along ``var2``.
Now these values can be fit with an analytic function that will interpolate through the 
middle signal region. The value of the analytic transfer function in the middle of the
``var2`` axis (call it :math:`f` (``var2``)) can then be used to do :math:`f(` ``var2`` :math:`)*D = C`.

Note that ``var2`` is not our measurement variable though. It's simply the variable in which we
measure the transfer function. The technical aspect of creating the QCD estimate in the C region is the same
as was done for the ABCD method but :math:`A/B` (which was a function of ``var3`` is now :math:`f`(``var2``).

We make the jump to 2DAlphabet by asking the interpolation method of Alphabet to measure :math:`f` in ``var3`` as well
so that it becomes :math:`f` (``var2``, ``var3``) and we get the best of both worlds. There are several advantages to this.

* There is typically a shape dependence in ``var3`` that vanilla Alphabet cannot handle.
* We still get the interpolation power of Alphabet that the ABCD method misses.
* The entire 2D space of ``var2`` vs ``var3`` can be used to constrain the analytic function.
* We get a complete model of the background in the dimensions we care about without needing to play games of weighting binning events along ``var3`` by the event's value in ``var2``.

The last point is the most powerful of the three because it allows us to build a model of the background without measuring
:math:`f` (``var2``, ``var3``) before fitting to data! Instead, it describes the binned PDF completely. We can then add it
into the total background model which consistents of backgrounds based on MC simulation
that can morph shape based on systematic uncertainties. Then :math:`f` can be measured simultaneous to the parameters determining the other backgrounds
in the likelihood fit to data (that is also performing the signal extraction)! This has nice benefits
like the fact that all correlations between :math:`f`'s parameters and the other nuisances will be accounted automatically.

One final note: :math:`f` (``var2``, ``var3``) does not have to be just the ratio between
the upper and lower regions of ``var1`` as described above (this method is commonly referred to as a "pass-fail ratio").
More terms can be added to the equation that determines the background model in the signal region, as long as they 
don't exist in a dimension outside of (``var2``, ``var3``). You could, for example, have 1D terms
in ``var2`` or factors coming from QCD simulation which, while not perfect, can still point you in the right direction.
To set this up requires writting a custom script using the 2DAlphabet API.

Other backgrounds and signal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As stated in the previous section, :math:`f` (``var2``, ``var3``) can be measured 
simultaneous to the other background contributions as well as the signal. 2DAlphabet handles this by
assuming all other backgrounds (and signal) are simulation based and that the simulation can (and is)
binned in the same (``var2``, ``var3``) space as the multijet/QCD model.
Uncertainties in the simulation can be represented by log-normal normalization terms
or shape templates which can morph the shape and normalization of the nominal MC template.
The details of accounting for systematic
uncertainties are explained in the :doc:`pages/config` section.

For more information, see the table of contents below.

Table of Contents
-----------------
.. toctree::
   :maxdepth: 2
   
   pages/how-to.rst
   pages/config.rst
   doctrees/modules.rst