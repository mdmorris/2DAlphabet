import ROOT
from ROOT import *

import array
from array import array

def fcn(npar, gin, f, par, iflag):
    lnL = 0
    failbin = 10000
    passbin = 1000
    rpf = par[0]

    passest = failbin*rpf

    lnL = passest - passbin*TMath.Log(passest) + TMath.Log(TMath.Factorial(passbin))#+ passbin*TMath.Log(passbin)-passbin

    print lnL
    print '\t 10000 * ' + str(rpf) + ' = ' + str(passest)

    f = 2*lnL

if __name__ == "__main__":

    minuit = TVirtualFitter.Fitter(0,1)
    minuit.SetFCN(fcn)
    minuit.SetParameter(0,'polX2Y2', 0.11, 0.01,0.,0.)
    arglist = array('d',2*[0.])
    arglist[0] = 10000
    minuit.ExecuteCommand("MIGRAD", arglist,0)
    up = 1.0
    minuit.SetErrorDef(1.0)
    arglist[0] = 0.
    minuit.ExecuteCommand("MINOS", arglist, 0);


