'''Create a single canvas of the 2D data.
System argument is the directory from which you want
to get the data to plot (from organized_hists.root)
'''

import ROOT, sys, CMS_lumi, tdrstyle
import header
ROOT.gROOT.SetBatch(1)
loc = sys.argv[1]
if len(sys.argv) == 3: CMS_subtext = sys.argv[2]
else: CMS_subtext = ''
print ('CMS_subtext ',CMS_subtext)
f = ROOT.TFile.Open(loc+'/organized_hists.root')
if loc[-1] == '/':
    tag = loc.split('/')[-2]
else:
    tag = loc.split('/')[-1]
for k in f.GetListOfKeys():
    if 'data' in k.GetName():
        print (k.GetName())
h = f.Get('data_obs_pass_FULL_%s'%tag)
header.makeCan('data_pass_2D',loc+'plots/',[h], bkglist=[],totalBkg=None,signals=[],colors=[],
            titles=[''],dataName='Data',bkgNames=[],signalNames=[],logy=True,
            rootfile=False,xtitle='m_{t} [GeV]',ytitle='m_{tW} [GeV]',ztitle='Events / bin',dataOff=False,
            datastyle='pe',year=1, addSignals=True, extraText=CMS_subtext)
