from optparse import OptionParser
import subprocess
import array
from  array import array

import ROOT
from ROOT import *

import header
from header import WaitForJobs, make_smooth_graph, Inter

gStyle.SetOptStat(0)
gROOT.SetBatch(kTRUE)

parser = OptionParser()

# parser.add_option('-t', '--tag', metavar='FILE', type='string', action='store',
#                 default   =   'dataBsOff',
#                 dest      =   'tag',
#                 help      =   'Tag ran over')
parser.add_option('-s', '--signals', metavar='FILE', type='string', action='store',
                default   =   'bstar_signalsLH.txt',
                dest      =   'signals',
                help      =   'Text file containing the signal names and their corresponding cross sections')
parser.add_option('-P', '--plotOnly', action="store_true",
                default   =   False,
                dest      =   'plotOnly',
                help      =   'Only plots if True')
parser.add_option('-b', '--blind', action="store_true",
                default   =   True,
                dest      =   'blind',
                help      =   'Only plot observed limit if false')
parser.add_option('-l', '--lumi', metavar='F', type='string', action='store',
                default       =       '',
                dest          =       'lumi',
                help          =       'Luminosity option')
parser.add_option('-m', '--mod', metavar='F', type='string', action='store',
                default       =       '',
                dest          =       'mod',
                help          =       'Modification to limit. For example, different handedness of the signal')

(options, args) = parser.parse_args()

# tag = options.tag

# Open signal file
signal_file = open(options.signals,'r')
# Read in names of project spaces as a list of strings and strip whitespace
signal_names = signal_file.readline().split(',')
signal_names = [n.strip() for n in signal_names]
# Read in mass as a list of strings, strip whitespace, and convert to ints
signal_mass = signal_file.readline().split(',')
signal_mass = [int(m.strip()) for m in signal_mass]
# Read in xsecs as a list of strings, strip whitespace, and convert to floats
signal_xsecs = signal_file.readline().split(',')
signal_xsecs = [float(x.strip()) for x in signal_xsecs]

# Initialize arrays to eventually store the points on the TGraph
x_mass = array('d')
y_limit = array('d')
y_mclimit  = array('d')
y_mclimitlow68 = array('d')
y_mclimitup68 = array('d')
y_mclimitlow95 = array('d')
y_mclimitup95 = array('d')

# For each signal
for this_index, this_name in enumerate(signal_names):
    # Setup call for one of the signal
    this_xsec = signal_xsecs[this_index]
    this_mass = signal_mass[this_index]
    this_output = TFile.Open(this_name+'/higgsCombineTest.Asymptotic.mH120.root')
    this_tree = this_output.Get('limit')
    
    # Set the mass (x axis)
    x_mass.append(float(this_mass))

    # Grab the cross section limits (y axis)
    for ievent in range(int(this_tree.GetEntries())):
        this_tree.GetEntry(ievent)
        
        # Nominal expected
        if this_tree.quantileExpected == 0.5:
            y_mclimit.append(this_tree.limit*this_xsec)
        # -1 sigma expected
        if round(this_tree.quantileExpected,2) == 0.16:
            y_mclimitlow68.append(this_tree.limit*this_xsec)
        # +1 sigma expected
        if round(this_tree.quantileExpected,2) == 0.84:
            y_mclimitup68.append(this_tree.limit*this_xsec)
        # -2 sigma expected
        if round(this_tree.quantileExpected,3) == 0.025:
            y_mclimitlow95.append(this_tree.limit*this_xsec)
        # +2 sigma expected
        if round(this_tree.quantileExpected,3) == 0.975:
            y_mclimitup95.append(this_tree.limit*this_xsec)

        # Observed (plot only if unblinded)
        if this_tree.quantileExpected == -1: 
            if not options.blind:
                y_limit.append(this_tree.limit*this_xsec)
            else:
                y_limit.append(0.0)
    
# Make Canvas and TGraphs (mostly stolen from other code that formats well)
climits = TCanvas("climits", "climits",700, 600)
climits.SetLogy(True)
climits.SetLeftMargin(.18)
climits.SetBottomMargin(.18)  

# NOT GENERIC
# if options.hand == 'LH':
#     cstr = 'L'
# elif options.hand == 'RH':
#     cstr = 'R'
# elif options.hand == 'VL':
#     cstr = 'LR'
# else:
#     cstr = ''
cstr = options.mod


TPT = ROOT.TPaveText(.20, .22, .5, .27,"NDC")
TPT.AddText("All-Hadronic Channel") # NOT GENERIC
TPT.SetFillColor(0)
TPT.SetBorderSize(0)
TPT.SetTextAlign(12)

# Expected
g_mclimit = TGraph(len(x_mass), x_mass, y_mclimit)
g_mclimit.SetTitle("")
g_mclimit.SetMarkerStyle(21)
g_mclimit.SetMarkerColor(1)
g_mclimit.SetLineColor(1)
g_mclimit.SetLineStyle(2)
g_mclimit.SetLineWidth(3)
g_mclimit.SetMarkerSize(0.)

# Observed
if not options.blind:
    g_limit = TGraph(len(x_mass), x_mass, y_limit)
    g_limit.SetTitle("")
    g_limit.SetMarkerStyle(0)
    g_limit.SetMarkerColor(1)
    g_limit.SetLineColor(1)
    g_limit.SetLineWidth(3)
    g_limit.SetMarkerSize(0.5) #0.5
    g_limit.GetXaxis().SetTitle("m_{b*_{"+cstr+"}} (GeV)")  # NOT GENERIC
    g_limit.GetYaxis().SetTitle("Upper Limit #sigma_{b*_{"+cstr+"}} #times B(b*_{"+cstr+"}#rightarrowtW) [pb]") # NOT GENERIC
    g_limit.GetYaxis().SetRangeUser(0., 80.)
    g_limit.GetXaxis().SetRangeUser(1, 3.2)
    g_limit.SetMinimum(3.0e-3) #0.005
    g_limit.SetMaximum(7000.)
else:
    g_mclimit.GetXaxis().SetTitle("m_{b*_{"+cstr+"}} (GeV)")  # NOT GENERIC
    g_mclimit.GetYaxis().SetTitle("Upper Limit #sigma_{b*_{"+cstr+"}} #times B(b*_{"+cstr+"}#rightarrowtW) [pb]") # NOT GENERIC
    g_mclimit.GetYaxis().SetRangeUser(0., 80.)
    g_mclimit.GetXaxis().SetRangeUser(1, 3.2)
    g_mclimit.SetMinimum(3.0e-3) #0.005
    g_mclimit.SetMaximum(7000.)
# Expected
# g_mclimit = TGraph(len(x_mass), x_mass, y_mclimit)
# g_mclimit.SetTitle("")
# g_mclimit.SetMarkerStyle(21)
# g_mclimit.SetMarkerColor(1)
# g_mclimit.SetLineColor(1)
# g_mclimit.SetLineStyle(2)
# g_mclimit.SetLineWidth(3)
# g_mclimit.SetMarkerSize(0.)
# g_mclimit.GetXaxis().SetTitle("M_{b*} (GeV/c^{2})")
# g_mclimit.GetYaxis().SetTitle("Upper Limit #sigma_{b*_{"+cstr+"}} #times b (pb)")
# g_mclimit.GetYaxis().SetTitleSize(0.03)
# g_mclimit.Draw("l")
# g_mclimit.GetYaxis().SetRangeUser(0., 80.)

# Will later be 1 and 2 sigma expected
g_mcplus = TGraph(len(x_mass), x_mass, y_mclimitup68)
g_mcminus = TGraph(len(x_mass), x_mass, y_mclimitlow68)

g_mc2plus = TGraph(len(x_mass), x_mass, y_mclimitup95)
g_mc2minus = TGraph(len(x_mass), x_mass, y_mclimitlow95)

# Theory line
graphWP = ROOT.TGraph()
graphWP.SetTitle("")
graphWP.SetMarkerStyle(23)
graphWP.SetMarkerColor(4)
graphWP.SetMarkerSize(0.5)
graphWP.GetYaxis().SetRangeUser(0., 80.)
graphWP.GetXaxis().SetRangeUser(1, 3.2)
graphWP.SetMinimum(3.0e-3) #0.005
graphWP.SetMaximum(7000.)
q = 0
for index,mass in enumerate(signal_mass):
    xsec = signal_xsecs[index]
    graphWP.SetPoint(q,    mass ,   xsec    )
    q+=1

graphWP.SetLineWidth(3)
graphWP.SetLineColor(4)


# Theory up and down unnecessary if not splitting PDF uncertainty into shape and norm
#
# # Theory up
# graphWPup = ROOT.TGraph()
# graphWPup.SetTitle("")
# graphWPup.SetMarkerStyle(23)
# graphWPup.SetMarkerColor(4)
# graphWPup.SetLineColor(4)
# graphWPup.SetLineWidth(2)
# graphWPup.SetMarkerSize(0.5)

# q = 0
# for bsmass in masses:
#     rt_xsec = mult*xsec_sig_dict[str(int(bsmass))][0]
#     graphWPup.SetPoint(q,    bsmass/1000. ,   rt_xsec    )
#     q+=1

# # Theory down
# graphWPdown = ROOT.TGraph()
# graphWPdown.SetTitle("")
# graphWPdown.SetMarkerStyle(23)
# graphWPdown.SetMarkerColor(4)
# graphWPdown.SetLineColor(4)
# graphWPdown.SetLineWidth(2)
# graphWPdown.SetMarkerSize(0.5)

# q = 0
# for bsmass in masses:
#     rt_xsec = mult*xsec_sig_dict[str(int(bsmass))][1]
#     graphWPdown.SetPoint(q,    bsmass/1000. ,   rt_xsec    )
#     q+=1

# graphWPup.SetLineStyle(2 )
# graphWPdown.SetLineStyle(2 )

# WPunc = make_smooth_graph(graphWPdown, graphWPup)
# WPunc.SetFillColor(4)
# WPunc.SetFillStyle(3004)
# WPunc.SetLineColor(0)

# 2 sigma expected
g_error95 = make_smooth_graph(g_mc2minus, g_mc2plus)
g_error95.SetFillColor(kYellow)
g_error95.SetLineColor(0)

# 1 sigma expected
g_error = make_smooth_graph(g_mcminus, g_mcplus)
g_error.SetFillColor( kGreen)
g_error.SetLineColor(0)

if not options.blind:
    g_limit.Draw('al')
    g_error95.Draw("lf")
    g_error.Draw("lf")
    g_mclimit.Draw("l")
    g_limit.Draw("l")
    graphWP.Draw("l")
    g_limit.GetYaxis().SetTitleOffset(1.4)

else:
    g_mclimit.Draw("al")
    g_error95.Draw("lf")
    g_error.Draw("lf")
    g_mclimit.Draw("l")
    graphWP.Draw("l")
    g_mclimit.GetYaxis().SetTitleOffset(1.4)
# graphWP.Draw("l")

# WPunc.Draw("lf")
# graphWPup.Draw("l")
# graphWPdown.Draw("l")


legend = TLegend(0.5, 0.45, 0.86, 0.84, '')
if not options.blind:
    legend.AddEntry(g_limit, "Observed", "l")
legend.AddEntry(g_mclimit, "Expected (95% CL)","l")
legend.AddEntry(g_error, "#pm 1 #sigma Expected", "f")
legend.AddEntry(g_error95, "#pm 2 #sigma Expected", "f")
legend.AddEntry(graphWP, "Theory b*_{"+cstr+"}", "l")   # NOT GENERIC
# legend.AddEntry(graphWPup, "Theory b*_{"+cstr+"} 1 #sigma uncertainty", "l")

legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetLineColor(0)

legend.Draw("same")

text1 = ROOT.TLatex()
text1.SetNDC()
text1.SetTextFont(42)
text1.DrawLatex(0.2,0.84, "#scale[1.0]{CMS, L = "+options.lumi+" pb^{-1} at  #sqrt{s} = 13 TeV}") # NOT GENERIC

TPT.Draw()      
climits.RedrawAxis()
climits.SaveAs("limits_combine_"+options.lumi+"pb_"+options.signals[:options.signals.find('.')]+".pdf")

# Finally calculate the intercept
expectedLimit = Inter(g_mclimit,graphWP)[0]
upLimit = Inter(g_mcminus,graphWP)[0]
lowLimit = Inter(g_mcplus,graphWP)[0]

print 'Expected limit: '+str(expectedLimit/1000.) + ' +'+str(upLimit/1000.-expectedLimit/1000.) +' -'+str(expectedLimit/1000.-lowLimit/1000.) + ' TeV' # NOT GENERIC
