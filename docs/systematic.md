---
layout: default
title: SYSTEMATIC
has_children: false
parent: Configuration files
nav_order: 2
---

# **`SYSTEMATIC`**
Because it bears repeating, please note the difference between this section, `SYSTEMATIC`,
and the list of `SYSTEMATICS` defined inside the `PROCESS` dictionary. The `SYSTEMATIC`
dictionary is a place to define as many systematics as a user may need. Similar to the
processes, each key in `SYSTEMATIC` is the name of the systematic in the analysis and each
is classified by a code that determines how the systematic will be treated. However, the
dictionary the user defines for a given systematic is different depending on what type it is.
The self-explanatory codes are 1 and 2 which are log-normal uncertainties on the normalization.
Less obvious are codes 2 and 3 which are for shape based uncertainties (and thus have
corresponding histograms) and are either in the same file as the process's nominal histogram
(code 2) or in a separate file (code 3). Additionally, they have a scale value which allows
the user to change the gaussian constraint on the shape. For no change in the constraint, use 1.0.
If you have templates representing a 2 $$\sigma$$ shift, use 0.5 to properly constrain
the associated nuisance parameter during the shape interpolation with Combine.

## Symmetric, log-normal
```json
"CODE": 0,
"VAL": <uncertainty> # float
```

## Asymmetric, log-normal

```json
"CODE": 1,
"VALUP": <+1 $$\sigma$$ uncertainty>, # float
"VALDOWN": <-1 $$\sigma$$ uncertainty> # float
```

## Shape based uncertainty, in same file as nominal histogram
```json
"CODE: 2",
"HISTPASS_UP": <name of hist (in same file as nominal hist) for +1 sig uncertainty in pass distribution>, # string
"HISTPASS_DOWN": <name of hist (in same file as nominal hist) for -1 sig uncertainty in pass distribution>, # string
"HISTFAIL_UP": <name of hist (in same file as nominal hist) for +1 sig uncertainty in fail distribution>, # string
"HISTFAIL_DOWN": <name of the hist (in same file as nominal hist) for -1 sig uncertainty in fail distribution>, # string
"SCALE": <scale value to change scale of nuisance constraint> # float
```

## Shape based uncertainty, in different file as nominal histogram. 

This is the more flexible
but also more complicated option. The user can specify files three different ways. 

1. By using `FILEUP` and `FILEDOWN` to pick a file that *every* process can pull the shape
   uncertainty histograms from. 
2. Use keys of the form `FILEUP_<proc>` 
   where `<proc>` matches the name of a process that is defined in the `PROCESS` dictionary and
   has this shape uncertainty associated with it. This allows each systematic and process to come
   from a separate file. 
3. Use keys of the form `FILEUP_*` where the * acts
   as a wild card for the process and must also exist in the file name where the process would
   normally be written. For example, if a $$t\bar{t}$$ distribution with +1 $$\sigma$$ pileup
   uncertainty is stored in `ttbar_pileup_up.root` and the corresponding signal distributions
   are in `signal_pileup_up.root`, one can use the key-value pair `"FILE_UP_*":"*_pileup_up.root"`.

The user can also specify histogram names in four different ways.
1. `HISTPASS` and `HISTFAIL` which allows the user to specify only two histogram
   names if they don't change between "up" and "down" shapes. 
2. The second is if the "up" and "down" shapes *do* have different histogram names
   and uses the form `HISTPASS_UP` `HISTFAIL_UP`. 
3. The totally generic way allows the user to use the form `HISTPASS_UP_<proc>` where `<proc>`
   matches the name of a process that is defined in the `PROCESS` dictionary and has this
   shape uncertainty associated with it. 
4. The "*" wildcard can be used in place of `<proc>` just as with the file keys.
   
Below is an example of the totally generic way (3).

```json
"CODE": 3,
"FILEUP_<proc>": </path/to/fileup_<proc>_up.root>, # string
"FILEDOWN_<proc>": </path/to/filedown_<proc>_down.root>, # string
"HISTPASS_UP": <name of hist (in same file as nominal hist) for +1 sig uncertainty in pass distribution>, # string
"HISTPASS_DOWN": <name of hist (in same file as nominal hist) for -1 sig uncertainty in pass distribution>, # string
"HISTFAIL_UP": <name of hist (in same file as nominal hist) for +1 sig uncertainty in fail distribution>, # string
"HISTFAIL_DOWN": <name of the hist (in same file as nominal hist) for -1 sig uncertainty in fail distribution>, # string
"SCALE": <scale value to change scale of nuisance constraint> # float
```

The various options lead to flexibility but the more organized you are, the easier it is to write the configuration file!