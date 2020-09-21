---
layout: default
title: Statistical tests - `run_Stats.py`
nav_order: 3
has_children: false
parent: How to run
---

# Statistical tests - `run_Stats.py`

The `run_Stats.py` tool hosts the infrastructure to run the majority
of necessary statistical tests on the model and fit. Each call
to `run_Stats.py` can only be run with a single "mode". The test
must be run on a project directory.
The available modes are represented by the following options.

- **`--gof`** Runs the goodness of fit test for the model in the 
  provided `--projDir`. If the fit in `--projDir` is blinded, the 
  goodness of fit test will be blinded as well. The saturated test
  statistic is used. The KS and AD test statistics are not available
  because they rely on CDFs which are not well defined for two dimensional
  distributions.
- **`--signalInjection <r>`** Injects the designated amount of signal
  `<r>` on top of the background-only model from `--projDir`.
- **`--ftest <option>`** Runs a Fischer test comparing the model in `--projDir`
  against a model with more parameters, `--altDir`. The process of running
  the test is split into several options so that pieces can be recycled 
  between comparisons:
  - **generate** Generates the toys `--projDir`. The `--altDir` does not
    need to be specified.
  - **fitMain** Fits the base/main model to the toys generated from itself
    and calculates the saturated test statistic for each.
  - **fitAlt** Fits the alternate model to the toys generated from the 
    base model and calculates the saturated test statistic for each.
  - **post** Collects the outputs and fitAlt and fitMain, plots, and calculates the 
    pvalue when comparing the fits to data against the fits to the toys.
  - **pvalue** Skips toys entirely and assumes toys follow an F-distribution.
    The pvalue is calculated when comparing data against the F-distribution. 
    This option is appropriate if you've first confirmed that the toys follow the
    F-distribution) which is plotted in the `post` option. 

* * *
### A note on toys

Statistical tests require that toys (pseudo-data) be generated and run. 
The number of toys to use can be specified with the  `-t/--toys` option.
Typically, 500 toys are a good number to generate and analyze.
The intricacies of the toy generation are handled by 2D Alphabet but 
in general, the toys are generated from a post-fit model (`fit_b` or a 
modified version depending on the "mode") and fit with the pre-fit model.
The seed for the toys can be changed with the `--seed` option.

* * *
### A note on using condor

For `run_Stats.py` specifically, there is the `--condor` option which
will run the Combine commands on condor batch nodes. When `--condor` 
is used, the jobs should be monitored by the user. Once the jobs are
finished, the same command should be run but with `--condor` replaced
by `--post` which will collect the returned tarballs, untar them 
temporarily, and generate the plots. 

The main advantage of using `--condor` is to split the toys among separate
jobs with the `--toyJobs` option which specifies into how many jobs to split
the toys. For example `-t 500 --toyJobs 50` would put 10 toys into each job.
The commands run in each job will have different seeds so that the toys
are not identical.

Note that even though the job outputs remain untarred (even after running
`--post`), the folders are not small. Depending on the model, they could be
between 1-2 GB. This can eat up space quickly on disk so be cognizant of how
many old tests you are saving!

Finally, the 2D Alphabet environment is already tarred and stored on
EOS for access to the jobs. If you make changes, you'll need to create
a new tarball of the environment and put it on EOS for the condor nodes
to access.

The `CondorHelper.py` can be used to submit jobs for other tools such as 
`run_Limit.py`. See [Condor Helper](condorhelper) for more information.
* * *
Running `python run_Stats.py --help` returns the following:

```
Options:
  -h, --help            show this help message and exit
  -d <dir>, --projDir=<dir>
                        Home of the project - has the cards, fit results, etc
  -a <dir>, --altDir=<dir>
                        Home of the alternative model that you'd like to
                        compare against the one in projDir
  -t <N>, --toys=<N>    Number of toys to generate - fed to combine
  --toyJobs=<N>         Number of jobs to split toy fitting into - fed to
                        condor job building
  --seed=<N>            Seed - fed to combine
  --rMin=<rMin>         Minimum bound on r (signal strength)
  --rMax=<rMax>         Minimum bound on r (signal strength)
  -w <name>             Model workspace (from text2workspace)
  --plotOnly            Only plot
  --dryrun              Dry run the combine commands to console
  --gof                 Perform goodness of fit test
  --signalInjection=<r>
                        Perform signal injection test
  --ftest=<option>      Perform F test. Options are 'generate', 'fitAlt',
                        'fitMain', 'post', 'pvalue'
  --post                Run in conjunction with diagnosticsWithToys or
                        signalInjection and condor jobs to process output
                        files
  --condor              Submit a condor job (only for signal injection and
                        diagnosticsWithToys currently
  --skipSnapshot        Use if you've already made a snapshot that you trust
                        for this project directory
```