---
layout: default
title: Configuration files
nav_order: 3
---

# Configuration files

While I've tried to keep the format of the JSON configuration file obvious,
it's not fair to assume that everyone else will feel the same way. I may use
phrases like, "The key 'TITLE' is defined inside 'PROCESS'." I'm aware this
is not correct since "PROCESS" is a key for a dictionary which has a key "TITLE"
and so "TITLE" can't be "inside" "PROCESS" but it simplifies the explanation
so it is the convention I'll use.  

The goal of the JSON is to have an easily configurable input that allows the
2D Alphabet software to read and organize analysis files and histograms to its
liking while also giving the user the ability to easily configure values like
bin sizes and ranges without lots of command line options. This means that while
the user is encouraged to add what they need, there are some keys and values
that must stay the same. These static strings are always in capital letters
to make them easy to distinguish. The six static sections are described below.

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
- `FILE:` the path to the file containing the nominal pass and fail histograms for the process (string);
- `HISTPASS:` the name of the histogram in the above file with the passing distribution (string);
- `HISTFAIL:` the name of the histogram in the above file with the failing distribution (string);
- `SYSTEMATICS:` a list of strings that correspond to the names of systematics
  in `SYSTEMATIC` (note `SYSTEMATIC` not `SYSTEMATICS`) that are applicable 
  to this process (list of strings);
- `CODE:` a way to classify the treatment of the process: 0 (signal), 1 (data), 2 (background simulation)

# **`SYSTEMATIC`**
Because it bears repeating, please note the difference between this section, `SYSTEMATIC`,
and the list of `SYSTEMATICS` defined inside the `PROCESS` dictionary. The `SYSTEMATIC`
dictionary is a place to define as many systematics as a user may need. Similar to the
processes, each key in `SYSTEMATIC` is the name of the systematic in the analysis and each
is classified by a code that determines how the systematic will be treated. However, the
dictionary the user defines for a given systematic is different depending on what type it is.
The self-explanatory types are:

- Symmetric, log-normal
    - `CODE: 0`
    - `VAL: <uncertainty>` (float)
- Asymmetric, log-normal
    - `CODE: 1`
    - `VALUP:` +1 $$\sigma$$ uncertainty (float)
    - `VALDOWN:` -1 $$\sigma$$ uncertainty (float)

Less obvious are codes 2 and 3 which are for shape based uncertainties (and thus have
corresponding histograms) and are either in the same file as the process's nominal histogram
(code 2) or in a separate file (code 3). Additionally, they have a scale value which allows
the user to change the gaussian constraint on the shape. For no change in the constraint, use 1.0.
If you have templates representing a 2 $$\sigma$$ shift, use 0.5 to properly constrain
the associated nuisance parameter during the shape interpolation with Combine.

## Shape based uncertainty, in same file as nominal histogram
- `CODE: 2`
- `HISTPASS_UP:` the name of the histogram (in the same file as the nominal histogram) 
    for +1$$\sigma$$ uncertainty in the pass distribution (string)
- `HISTPASS_DOWN:` the name of the histogram (in the same file as the nominal histogram)
    for -1$$\sigma$$ uncertainty in the pass distribution (string)
- `HISTFAIL_UP:` the name of the histogram (in the same file as the nominal histogram)
    for +1$$\sigma$$ uncertainty in the fail distribution (string)
- `HISTFAIL_DOWN:` the name of the histogram (in the same file as the nominal histogram)
    for -1$$\sigma$$ uncertainty in the fail distribution (string)
- `SCALE:` a scale value which allows the user to change the normalization of the shape (float)

## Shape based uncertainty, in different file as nominal histogram. 

This is the more flexible
but also more complicated option. The user can specify files three different ways. 

1. By using `FILEUP:` and `FILEDOWN:` to pick a file that *every* process can pull the shape
   uncertainty histograms from. 
2. Use keys of the form `FILEUP_<proc>:` 
   where `<proc>` matches the name of a process that is defined in the `PROCESS` dictionary and
   has this shape uncertainty associated with it. This allows each systematic and process to come
   from a separate file. 
3. Use keys of the form `FILEUP_*:` where the * acts
   as a wild card for the process and must also exist in the file name where the process would
   normally be written. For example, if a $$t\bar{t}$$ distribution with +1 $$\sigma$$ pileup
   uncertainty is stored in `ttbar_pileup_up.root` and the corresponding signal distribtutions
   are in `signal_pileup_up.root`, one can use the key-value pair `"FILE_UP_*":"*_pileup_up.root"`.

The user can also specify histogram names in four different ways.
1. `HISTPASS` and `HISTFAIL` which allows the user to specify only two histogram
   names if they don't change between "up" and "down" shapes. 
2. The second is if the "up" and "down" shapes \textit{do} have different histogram names
   and uses the form `HISTPASS_UP` `HISTFAIL_UP`. 
3. The totally generic way allows the user to use the form `HISTPASS_UP_<proc>` where `<proc>`
   matches the name of a process that is defined in the `PROCESS` dictionary and has this
   shape uncertainty associated with it. 
4. The "*" wildcard can be used in place of `<proc>` just as with the file keys.
   
Below is the layout of the totally generic way (3).

- "CODE: 3" (int)
- `FILEUP_<proc>:` `/path/to/fileup_<proc>.root` which contains the +1 $$\sigma$$
  uncertainty histogram for <proc> (string)
- `FILEDOWN_<proc>:` `/path/to/filedown_<proc>.root` which contains the -1 $$\sigma$$
  uncertainty histogram for <proc> (string)
- `HISTPASS_UP_<proc>:` the name of histogram for `<proc>` in `/path/to/fileup_<proc>.root`
  for +1 $$\sigma$$ uncertainty in the pass distribution (string)
- `HISTPASS_DOWN_<proc>:` the name of the histogram for `<proc>` in `/path/to/filedown_<proc>.root`
  for -1 $$\sigma$$ uncertainty in the pass distribution (string)
- `HISTFAIL_UP_<proc>:` the name of histogram for `<proc>` in `/path/to/fileup_<proc>.root`
  for +1$\sigma$ uncertainty in the fail distribution (string)
- `HISTFAIL_DOWN_<proc>:` the name of the histogram for `<proc>` in `/path/to/filedown_<proc>.root`
  for -1 $$\sigma$$ uncertainty in the fail distribution (string)
- `SCALE:` a scale value which allows the user to change the normalization of the shape (float)

This scheme is quite flexible. However, the more organized you are, the easier it is to write a configuration file.

# **`BINNING`**
One set of important user defined values is the 2D binning of the space being analyzed. This dictionary is the opportunity to define the axes binning of the user's space. The binning values are split into x and y axis definitions where the x-axis describes the variable whose signal region is blinded. Note that it \textit{is} possible to rebin and reduce the ranges of the input axes. However, this is mainly for quick tests to remove a bin or reduce the number of bins. For permanant situations, my recommendation is to remake the input histograms to the desired binning and have the configuration file match (it's one less thing that can go wrong!)
%the \verb"LOW" and \verb"HIGH" bin edges for the \verb"Y" axis that are defined here \textit{must} be consistent with the user's input histograms. Additionally, 
        
The binning in each axis cannot be beyond the range of the input histogram (hopefully, this is obvious) but it can be a subset of the input range. The number of bins are restricted to be equal to or less than the input number of bins (hopefully, this is also obvious). Additionally, the signal bounds only apply to the \verb"X" axis and must exist within the \verb"X" axis range defined in the configuration file. The user must specify if they want to fit their 2D space by blinding the signal region via the \verb"BLINDED" key which can take only boolean \verb"false" and \verb"true". 

        <!-- \begin{itemize}
            \item \verb"X"
            \begin{itemize}
                \item \verb"NAME:" name of your variable on the x-axis (string)
                \item \verb"TITLE:" title that you'd like to appear in plots with this axis (string)
                \item \verb"LOW:" lower bound of x-axis (int)
                \item \verb"HIGH:" upper bound of x-axis (int)
                \item \verb"NBINS:" number of x bins from LOW to HIGH (int)
                \item \verb"SIGSTART:" lower bound of signal region of x-axis (int)
                \item \verb"SIGEND:" upper bound of signal region of x-axis (int)
                \item \verb"BLINDED:" blinds or unblinds the signal region during the $R_{P/F}$ fit (bool)
            \end{itemize}
            \item \verb"Y"
            \begin{itemize}
                \item \verb"NAME:" name of your variable on the y-axis (string)
                \item \verb"TITLE:" title that you'd like to appear in plots with this axis (string)
                \item \verb"LOW:" lower bound of y-axis (int)
                \item \verb"HIGH:" upper bound of y-axis (int)
                \item \verb"NBINS:" number of x bins from \verb"LOW" to \verb"HIGH" (int)
            \end{itemize}
        \end{itemize}

    \subsubsection{FIT}
        The other set of important user defined values is the fit parameters for the transfer function from the fail sideband to the passing (or pass-fail ratio). The 2D fit of the transfer function assumes a polynomial in both directions. The order of the polynomials are up to the user (technically they are capped at order 10) as long as each parameter is defined. For example, if the highest order parameter defined is 2 in x and 1 in y, the user needs to define six parameters. Expanding the number of functions available to the user is currently in development \footnote{please let me know if there's something you'd like me to include and I can do my best to implement it}. 

        The user can convey the desired polynomial orders by separating the polynomial into \verb"XFORM" and \verb"YFORM" which are the shapes in the \verb"X" and \verb"Y" directions, respectively. For example, one could write \verb"(@1+@*2x)" for the \verb"XFORM" and \verb"(@1+@2*y)" for the \verb"YFORM" where \verb"@1, @2" in the \verb"XFORM" would correspond to \verb"X1" and \verb"X2" in the configuration file and \verb"@1, @2" in the \verb"YFORM" would correspond to \verb"Y1" and \verb"Y2" in the configuration file. Notice that the convention is that the numbering of the parameters in \verb"XFORM" and \verb"YFORM" must start at 1 - not 0.

        Finally, for each parameter, the user must specify \verb"NOMINAL", \verb"LOW", and \verb"HIGH" values that correspond to a guess of the initial value and range for the parameter while it floats in the fit. 

        When deciding on the nominal value and range for each parameter you should consider the following
        \begin{itemize}
        
        \end{itemize}

        An example is given below.
        \begin{itemize}
            \item \verb"FORM:" (@1+@2*x+@3*y+@4*x*y)
            \item \verb"X0Y0"
            \begin{itemize}
                \item \verb"NOMINAL:" nominal value (float)
                \item \verb"LOW:" lower bound (float)
                \item \verb"HIGH:" upper bound (float)
            \end{itemize}
            \item \verb"X1Y0"
            \begin{itemize}
                \item \verb"NOMINAL:" nominal value (float)
                \item \verb"LOW:" lower bound (float)
                \item \verb"HIGH:" upper bound (float)
            \end{itemize}
            \item \verb"X0Y1"
            \begin{itemize}
                \item \verb"NOMINAL:" nominal value (float)
                \item \verb"LOW:" lower bound (float)
                \item \verb"HIGH:" upper bound (float)
            \end{itemize}
            \item \verb"X1Y1"
            \begin{itemize}
                \item \verb"NOMINAL:" nominal value (float)
                \item \verb"LOW:" lower bound (float)
                \item \verb"HIGH:" upper bound (float)
            \end{itemize}
            \item ...
        \end{itemize}

    \subsubsection{GLOBAL}
        This dictionary is designed to help users with large configuration files by allowing them to create JSON-wide variables. For example, if all of your files are located in \verb"/long/path/to/my/variables/", you can store this string in the GLOBAL dictionary with a custom key (let's say \verb"dir/"). Now instead of having to write the full directory path for every process and systematic, the user can just write \verb"dir/". This simplifies the JSON and also has the standard advantages of using variables over several instances of the same object.

        This works by searching all strings in the JSON for instances of each key in \verb"GLOBAL" and replacing the key with its corresponding dictionary value.

        Thus, the user must be careful they don't accidentally use strings in the JSON that are identical to keys in \verb"GLOBAL" but that should be unchanged. This means keys in \verb"GLOBAL" should be descriptive (single character keys would be a bad idea). 

    \subsubsection{OPTIONS}
        This section is used to supply additional options to the fit that do not fit into any of the other five categories. A description of the possible options and the input they take are described below.

        \begin{itemize}
            \item 
        \end{itemize}}
 -->
