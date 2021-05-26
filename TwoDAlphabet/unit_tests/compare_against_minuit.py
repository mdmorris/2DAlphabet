import ROOT
from ROOT import *

import sys

test_type = sys.argv[1]

combine_file = TFile.Open('../'+test_type+'100kv10k_nosig_ff/postfitshapes_b.root')
minuit_file = TFile.Open('result_minuit_'+test_type+'.root')

combine_estimate = combine_file.Get('pass_postfit/TotalBkg')
minuit_estimate = minuit_file.Get('final_est')

diff = combine_estimate.Clone('diff')
diff.Add(minuit_estimate,-1)

diff.Draw('surf')

raw_input('waiting')

# Now compare rpf paramters
# Minuit setup
minuit_fit = minuit_file.Get('newrpf')

# Combine setup
combine_fitresult_file = TFile.Open('../'+test_type+'100kv10k_nosig_ff/fitDiagnostics.root')
combine_fitresult = combine_fitresult_file.Get('fit_b')
param_final = combine_fitresult.floatParsFinal()
coeffs_final = param_final.selectByName('polyCoeff_*')
coeffIter_final = coeffs_final.createIterator()
coeff_final = coeffIter_final.Next()

if test_type == 'flat':
    O = 0
elif test_type == 'lin':
    O = 1
elif test_type == 'quad':
    O = 2

# Loop
print 'Name | Combine  | Minuit  '
print '----------------------------'
while coeff_final:
    coeffName = coeff_final.GetName().split('_')[1]
    oX = int(coeffName[1])
    oY = int(coeffName[3])

    minuit_param_index = (O+1)*oY+oX

    print 'X'+str(oX)+'Y'+str(oY)+' |',
    
    print "%.6f" % coeff_final.getValV() + ' | ',
    coeff_final = coeffIter_final.Next()

    print "%.6f" % minuit_fit.GetParameter(minuit_param_index)

combine_fitresult_file.Close()
minuit_file.Close()
combine_file.Close()
