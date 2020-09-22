---
layout: default
title: FIT
has_children: false
parent: Configuration files
nav_order: 4
---

# `FIT`
This section defines the values of the fit parameters for the transfer
function from the fail region to the passing (or pass-fail ratio). The
2D fit can accommodate any functional form. Each parameter in equation
should be marked with a @ symbol (RooFit convention) starting at 0. Each
parameter should also be specified with range and if desired an error 
(which does not establish a constraint but just an initial step size 
in the minimization). 

An example section which uses a polynomial in both directions (of 
different orders) is provided below. The range of each parameter
is -10 to 10 with a starting (`NOMINAL`) value of 1.0 and step size
(`ERROR`) of 0.1.

To accommodate for the potential of a negative transfer function value
in all or part of the 2D space (which would be unphysical and make 
Combine very mad), the `FORM` provided is wrapped in `max(FORM, 1e-9)`.
Thus, one should always check in the plots that the final parameters
do not produce a ~0 transfer function in large portions of the space.
The fit result will not be stable. This can be addressed by reducing
the ranges (`MIN`/`MAX`) of the parameters or changing the functional 
form.

```json
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
```