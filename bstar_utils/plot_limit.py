from optparse import OptionParser
import subprocess, sys, array, math, os
from array import array

sys.path.append(os.getcwd())

import ROOT
from ROOT import *

import header
from header import WaitForJobs, make_smooth_graph, Inter
import tdrstyle, CMS_lumi

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
parser.add_option('--unblind', action="store_false",
                default   =   True,
                dest      =   'blind',
                help      =   'Only plot observed limit if false')
parser.add_option('--drawIntersection', action="store_true",
                default   =   False,
                dest      =   'drawIntersection',
                help      =   'Draw intersection values')
parser.add_option('-l', '--lumi', metavar='F', type='string', action='store',
                default       =       '137.44',
                dest          =       'lumi',
                help          =       'Luminosity option')
parser.add_option('-m', '--mod', metavar='F', type='string', action='store',
                default       =       '',
                dest          =       'mod',
                help          =       'Modification to limit title on y-axis. For example, different handedness of the signal')

(options, args) = parser.parse_args()

# tag = options.tag

# Open signal file
signal_file = open(options.signals,'r')
# Read in names of project spaces as a list of strings and strip whitespace
signal_names = signal_file.readline().split(',')
signal_names = [n.strip() for n in signal_names]
# Read in mass as a list of strings, strip whitespace, and convert to ints
signal_mass = signal_file.readline().split(',')
signal_mass = [float(m.strip())/1000 for m in signal_mass]
# Read in xsecs as a list of strings, strip whitespace, and convert to floats
theory_xsecs = signal_file.readline().split(',')
theory_xsecs = [float(x.strip()) for x in theory_xsecs]
# Read in xsecs that samples were normalized to
signal_xsecs = signal_file.readline().split(',')

if 'tprime' in options.signals:
    signal_xsecs = [float(x.strip()) for x in signal_xsecs]
elif 'bprime' in options.signals:
    signal_xsecs = [float(x.strip()) for x in signal_xsecs]
else:
    signal_xsecs = [float(x.strip()) for x in signal_xsecs]

# Get theory PDF variations - BSTAR SPECIFIC
bs_path = '/uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_10_2_0/src/BStar13TeV/SFs/'

if 'LH' in options.signals:
    hand = 'LH'
elif 'RH' in options.signals:
    hand = 'RH'
elif 'VL' in options.signals:
    hand = 'VL'

if 'bprime' in options.signals:
    label = "B'"
    theory_var_file = TFile.Open(bs_path+'pdf_norm_uncertainties_TBprime.root')
    sigma_max = 1
    sigma_min = 1e-5
    mass_min = 1.4
    mass_max = 1.8
elif 'tprime' in options.signals:
    label = "T'"
    theory_var_file = TFile.Open(bs_path+'pdf_norm_uncertainties_TBprime.root')
    sigma_max = 10
    sigma_min = 5e-3
    mass_min = 1.4#0.59
    mass_max = 1.8#1.81
elif 'bstar' in options.signals:
    label = 'b*'
    theory_var_file = TFile.Open(bs_path+'pdf_norm_uncertainties_bstar.root')
    sigma_max = 20
    sigma_min = 0.0000008
    mass_min = 1.2
    mass_max = 4.2
else: label = 'X'

# Initialize arrays to eventually store the points on the TGraph
x_mass = array('d')
y_limit = array('d')
y_mclimit  = array('d')
y_mclimitlow68 = array('d')
y_mclimitup68 = array('d')
y_mclimitlow95 = array('d')
y_mclimitup95 = array('d')

tdrstyle.setTDRStyle()

# For each signal
for this_index, this_name in enumerate(signal_names):
    # Setup call for one of the signal
    this_xsec = signal_xsecs[this_index]
    this_mass = signal_mass[this_index]
    this_output = TFile.Open(this_name+'/higgsCombineTest.AsymptoticLimits.mH120.root')
    if not this_output: continue
    this_tree = this_output.Get('limit')
    print (this_mass)
    # Set the mass (x axis)
    x_mass.append(this_mass)
    # Grab the cross section limits (y axis)
    for ievent in range(int(this_tree.GetEntries())):
        this_tree.GetEntry(ievent)
        
        print '\t%s = %s'%(this_tree.quantileExpected,this_tree.limit)
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
climits.SetLeftMargin(.15)
climits.SetBottomMargin(.15)  
climits.SetTopMargin(0.1)
climits.SetRightMargin(0.05)

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

gStyle.SetTextFont(42)
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
    print 'Not blinded'
    g_limit = TGraph(len(x_mass), x_mass, y_limit)
    g_limit.SetTitle("")
    g_limit.SetMarkerStyle(7)
    g_limit.SetMarkerColor(1)
    g_limit.SetLineColor(1)
    g_limit.SetLineWidth(2)
    g_limit.SetMarkerSize(1) #0.5
    g_limit.GetYaxis().SetRangeUser(0., 80.)
    g_limit.GetXaxis().SetRangeUser(mass_min, mass_max)
    g_limit.SetMinimum(sigma_min) #0.005
    g_limit.SetMaximum(sigma_max)
else:
    print 'Blinded'
    g_mclimit.GetXaxis().SetTitle("m_{"+label+"_{"+cstr+"}} [TeV]")  # NOT GENERIC
    g_mclimit.GetYaxis().SetTitle("#sigma_{"+label+"_{"+cstr+"}} #times B("+label+"_{"+cstr+"}#rightarrow tW) (pb)") # NOT GENERIC
    g_mclimit.GetYaxis().SetRangeUser(0., 80.)
    g_mclimit.GetXaxis().SetRangeUser(mass_min, mass_max)
    g_mclimit.SetMinimum(sigma_min) #0.005
    g_mclimit.SetMaximum(sigma_max)
# Expected
# g_mclimit = TGraph(len(x_mass), x_mass, y_mclimit)
# g_mclimit.SetTitle("")
# g_mclimit.SetMarkerStyle(21)
# g_mclimit.SetMarkerColor(1)
# g_mclimit.SetLineColor(1)
# g_mclimit.SetLineStyle(2)
# g_mclimit.SetLineWidth(3)
# g_mclimit.SetMarkerSize(0.)
# g_mclimit.GetXaxis().SetTitle("M_{b*} (TeV/c^{2})")
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
graphWP.GetXaxis().SetRangeUser(mass_min, mass_max)
graphWP.SetMinimum(sigma_min) #0.005
graphWP.SetMaximum(sigma_max)
for index,mass in enumerate(signal_mass):
    xsec = theory_xsecs[index]
    graphWP.SetPoint(index,    mass,   xsec    )

graphWP.SetLineWidth(1)
graphWP.SetLineColor(4)


# Theory up and down unnecessary if not splitting PDF uncertainty into shape and norm
#
# Theory up
graphWPup = ROOT.TGraph()
graphWPup.SetTitle("")
graphWPup.SetMarkerStyle(23)
graphWPup.SetMarkerColor(4)
graphWPup.SetLineColor(4)
graphWPup.SetLineWidth(2)
graphWPup.SetMarkerSize(0.5)

# Theory down
graphWPdown = ROOT.TGraph()
graphWPdown.SetTitle("")
graphWPdown.SetMarkerStyle(23)
graphWPdown.SetMarkerColor(4)
graphWPdown.SetLineColor(4)
graphWPdown.SetLineWidth(2)
graphWPdown.SetMarkerSize(0.5)

if label == 'b*': years = ['16','17','18']
else: years = ['16']

for i,bsmass in enumerate(signal_mass):
    pdf_up = 0
    pdf_down = 0
    for y in years:
        if label == 'b*': this_signal_name = 'signal%s%s_%s'%(hand,int(bsmass*1000),y)
        elif label == "B'": this_signal_name = 'Bprime%s%s_%s'%(hand,int(bsmass*1000),y)
        elif label == "T'": this_signal_name = 'Tprime%s%s_%s'%(hand,int(bsmass*1000),y)

        if hand != 'VL':
            pdf_hist = theory_var_file.Get(this_signal_name)
            pdf_up += (pdf_hist.GetBinContent(1)-1)**2
            pdf_down += (1-pdf_hist.GetBinContent(2))**2
        else:
            pdf_hist_LH = theory_var_file.Get(this_signal_name.replace('VL','LH'))
            pdf_hist_RH = theory_var_file.Get(this_signal_name.replace('VL','RH'))
            pdf_up += (pdf_hist_LH.GetBinContent(1)-1)**2 + (pdf_hist_RH.GetBinContent(1)-1)**2
            pdf_down += (1-pdf_hist_LH.GetBinContent(2))**2 + (1-pdf_hist_RH.GetBinContent(2))**2
    
    theory_up = (1+math.sqrt(pdf_up))*theory_xsecs[i]
    theory_down = (1-math.sqrt(pdf_down))*theory_xsecs[i]
    print ('%s %s +%s -%s'%(bsmass,theory_xsecs[i],theory_up,theory_down))
    graphWPup.SetPoint(i, bsmass, theory_up)
    graphWPdown.SetPoint(i, bsmass, theory_down)

# graphWPup.SetLineStyle(2 )
# graphWPdown.SetLineStyle(2 )

WPunc = make_smooth_graph(graphWPdown, graphWPup)
WPunc.SetFillColor(4)
WPunc.SetFillStyle(3004)
WPunc.SetLineColor(4)
WPunc.SetLineWidth(1)

# 2 sigma expected
g_error95 = make_smooth_graph(g_mc2minus, g_mc2plus)
g_error95.SetFillColor(kOrange)
g_error95.SetLineColor(0)

# 1 sigma expected
g_error = make_smooth_graph(g_mcminus, g_mcplus)
g_error.SetFillColor( kGreen+1)
g_error.SetLineColor(0)

if not options.blind:

    g_limit.GetXaxis().SetTitle("m_{"+label+"_{"+cstr+"}} [TeV]")  # NOT GENERIC
    g_limit.GetYaxis().SetTitle("#sigma_{"+label+"_{"+cstr+"}} #times B("+label+"_{"+cstr+"}#rightarrow tW) (pb)") # NOT GENERIC
    g_limit.GetXaxis().SetTitleSize(0.055)
    g_limit.GetYaxis().SetTitleSize(0.05)
    g_limit.Draw('ap')
    g_error95.Draw("lf")
    g_error.Draw("lf")
    g_mclimit.Draw("l")
    g_limit.Draw("lp")
    # graphWP.Draw("l")
    # WPunc.Draw("f")
    # graphWPup.Draw("l")
    # graphWPdown.Draw("l")

    g_limit.GetYaxis().SetTitleOffset(1.5)
    g_limit.GetXaxis().SetTitleOffset(1.25)

else:
    g_mclimit.GetXaxis().SetTitle("m_{"+label+"_{"+cstr+"}} [TeV]")  # NOT GENERIC
    g_mclimit.GetYaxis().SetTitle("#sigma_{"+label+"_{"+cstr+"}} B("+label+"_{"+cstr+"}#rightarrow tW) (pb)") # NOT GENERIC
    g_limit.GetXaxis().SetTitleSize(0.055)
    g_limit.GetYaxis().SetTitleSize(0.05)
    g_mclimit.Draw("al")
    g_error95.Draw("lf")
    g_error.Draw("lf")
    g_mclimit.Draw("l")
    # graphWP.Draw("l")
    g_mclimit.GetYaxis().SetTitleOffset(1.5)
    g_mclimit.GetXaxis().SetTitleOffset(1.25)
    
graphWP.Draw("l")
WPunc.Draw("fl")
graphWPup.Draw("l")
graphWPdown.Draw("l")

# Finally calculate the intercept
expectedMassLimit,expectedCrossLimit = Inter(g_mclimit,graphWP) #if len(Inter(g_mclimit,graphWP)) > 0 else -1.0
upLimit,trash = Inter(g_mcminus,graphWP) if len(Inter(g_mcminus,graphWP)) > 0 else -1.0
lowLimit,trash = Inter(g_mcplus,graphWP) if len(Inter(g_mcplus,graphWP)) > 0 else -1.0

expLine = TLine(expectedMassLimit,g_mclimit.GetMinimum(),expectedMassLimit,expectedCrossLimit)
expLine.SetLineStyle(2)
expLine.Draw()

if options.drawIntersection:
    expLineLabel = TPaveText(expectedMassLimit-300, expectedCrossLimit*2, expectedMassLimit+300, expectedCrossLimit*15, "NB")
    expLineLabel.SetFillColorAlpha(kWhite,0)
    expLineLabel.AddText(str(int(expectedMassLimit))+' TeV')
    expLineLabel.Draw()

print 'Expected limit: '+str(expectedMassLimit) + ' +'+str(upLimit-expectedMassLimit) +' -'+str(expectedMassLimit-lowLimit) + ' TeV' # NOT GENERIC
if not options.blind:
    obsMassLimit,obsCrossLimit = Inter(g_limit,graphWP) if len(Inter(g_limit,graphWP)) > 0 else -1.0
    print 'Observed limit: '+str(obsMassLimit) + ' TeV'

    obsLine = TLine(obsMassLimit,g_mclimit.GetMinimum(),obsMassLimit,obsCrossLimit)
    obsLine.SetLineStyle(2)
    obsLine.Draw()

    if options.drawIntersection:
        obsLineLabel = TPaveText(obsMassLimit-300, obsCrossLimit*3, obsMassLimit+300, obsCrossLimit*12,"NB")
        obsLineLabel.SetFillColorAlpha(kWhite,0)
        obsLineLabel.AddText(str(int(obsMassLimit))+' TeV')
        obsLineLabel.Draw()

# Legend and draw
gStyle.SetLegendFont(42)
legend = TLegend(0.6, 0.5, 0.91, 0.87, '')
legend.SetHeader("95% CL upper limits")
if not options.blind:
    legend.AddEntry(g_limit, "Observed", "l")
legend.AddEntry(g_mclimit, "Median expected","l")
legend.AddEntry(g_error, "68% expected", "f")
legend.AddEntry(g_error95, "95% expected", "f")
legend.AddEntry(WPunc, "Theory "+label+"_{"+cstr+"}", "lf")   # NOT GENERIC
# legend.AddEntry(WPunc, "Theory "+label+"_{"+cstr+"} PDF uncertainty", "f")

legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetLineColor(0)

legend.Draw("same")

# text1 = ROOT.TLatex()
# text1.SetNDC()
# text1.SetTextFont(42)
# text1.DrawLatex(0.17,0.88, "#scale[1.0]{CMS, L = "+options.lumi+" fb^{-1} at  #sqrt{s} = 13 TeV}") # NOT GENERIC

# TPT.Draw()      
climits.RedrawAxis()

CMS_lumi.extraText = 'Preliminary'
CMS_lumi.lumiTextSize     = 0.5

CMS_lumi.cmsTextSize      = 0.8
CMS_lumi.CMS_lumi(climits, 1, 11)

climits.SaveAs("prelimits_combine_"+options.lumi.replace('.','p')+"fb_"+options.signals[options.signals.find('/')+1:options.signals.find('.')]+'_'+cstr+".pdf")


