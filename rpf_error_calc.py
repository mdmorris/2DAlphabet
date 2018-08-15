###########################################################################################
# rpf_error_calc.py - Written 8/13/18 by Lucas Corcodilos 								  #
# ------------------------------------------------------- 								  #
# Takes as input the RooFitResult from a fit, the polynomial coefficient RRV names, 	  #
# the binning of the 2D space, and the number of samples to run. 						  #
#																						  #
# The parameters are randomized from the fit result and a 3D histogram is filled where 	  #
# the 3rd dimension is the value of the Rp/f from each sample. A TH1F is created from 	  #
# each 2D bin (thus the axis of the TH1F is the distribution of Rp/f values in that 2D 	  #
# bin) and the RMS is taken as the error on the Rp/f in that bin.						  # 
#																						  #
# The function returns a TH2F where each bin is the nominal value of the Rp/f with errors #
# as calculated as above.																  #
###########################################################################################


import ROOT
from ROOT import *

def main(fitresult,polycoeffs,binning,samples=100):
	
