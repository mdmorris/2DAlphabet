#####################################################################################################
# INCOMPLETE AND NOT GENERIC - bad lines marked NOT GENERIC                                         #
#####################################################################################################
# do_full_limits.py - written by Lucas Corcodilos, 5/31/18                                          #
# --------------------------------------------------------                                          #
# This is the wrapper of 2DAlphabet.py that can run limits for multiple signals in the current      #
# session (condor support is in development)                                                        #
# Afterwards, the limits are plotted on a traditional 'Brazilian flag' limit plot.                  #
# This is NOT generic and requires some work from the user. The user must use the command line      #
# options to define a config file template, a file containing the list of signal names on the top   #
# line and the list of cross sections for each signal on the next line (in the same order as the    #
# signal name listand both comma separated), and the placeholder in the config template that will   #
# be replaced with the signal names (default is TEMPLATE).                                          #
#                                                                                                   #
#                                                                                                   #
# WARNING AFTER CONDOR IS SUPPORTED: In order to use condor, this script requires some prep by the  #
# user. In order for condor to run, it needs to have the 2D Alphabet environment setup WITH the     #
# workspaces that are created in the middle of this script. There are obviously several ways to do  #
# this but the simplest is to tarball the entire CMSSW release that the user is working in.         #
# However, it's entirely possible that you have other things in this CMSSW that condor does NOT     #
# So below is an explanation of what I've setup to get the cleanest and smallest environment in     #
# condor without too much confusion. Users could obviously approach this problem in different ways. #
# This script will need to be edited to reflect those changes.                                      #
#                                                                                                   #
# 2DAlpha_env Setup (not working)                                                                   #
# -----------------                                                                                 #
# - Copy your entire CMSSW release to a nice corner where you won't accidentally play with it.      #
#   I put mine in a folder named 2DAlpha_env                                                        #
# - Remove anything in the release that 2D Alphabet does not need (this should just be 2DAlphabet/, #
#   HiggsAnalysis/, and any folders with your 2D preselected root files).                           #
# - Setup an alias `rsync2Denv` that uses rsync to move changes in your working 2DAlphabet          #
#   directory to 2DAlpha_env, makes a .tgz archive, and then uses xrdcp to send it to EOS. We don't #
#   need to worry about this for the root file folder and HiggsAnalysis/ because these won't be     #
#   changing frequently.                                                                            #    
# --- In your .tcshrc `alias rsync2Denv 'rsync -azP --delete /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet /uscms_data/d3/lcorcodi/BStar13TeV/2DAlpha_env/CMSSW_7_4_7_patch2/src/2DAlphabet; tar -climitszf /uscms_data/d3/lcorcodi/BStar13TeV/2DAlpha_env/2DAlpha_env.tgz /uscms_data/d3/lcorcodi/BStar13TeV/2DAlpha_env/CMSSW_7_4_7_patch2; ; xrdcp /uscms_data/d3/lcorcodi/BStar13TeV/2DAlpha_env/2DAlpha_env.tgz root://cmseos.fnal.gov//store/user/lcorcodi/'`
# --- The first time you run this, it will sync everything over and will take some time but will    #
#     only sync changes afterwards                                                                  #
# --- Make sure the EOS location matches where you copy from in --envString                         #
# - This script will call `rsync2Denv` before running on condor to ensure that the environment is   #
#   up to date with all of the latest workspaces.                                                   #
# - If the workspaces do not need to be rerun and the user would just like to rerun the condor jobs,#
#   the grid_sub_<tag>.csh script can just be called.                                               #
#####################################################################################################
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

parser.add_option('-t', '--template', metavar='FILE', type='string', action='store',
                default   =   '',
                dest      =   'template',
                help      =   'JSON template file to be changed for each signal. Name should be "input_<tag>_template.json" where tag will be used to organize outputs')
parser.add_option('-d', '--dummy', metavar='FILE', type='string', action='store',
                default   =   'TEMPLATE',
                dest      =   'dummy',
                help      =   'Dummy name to be replaced in the template with the signal names')
parser.add_option('-s', '--signals', metavar='FILE', type='string', action='store',
                default   =   '',
                dest      =   'signals',
                help      =   'Text file containing the signal names and their corresponding cross sections')
parser.add_option('-P', '--plotOnly', action="store_true",
                default   =   False,
                dest      =   'plotOnly',
                help      =   'Only plots if True')
parser.add_option('-e', '--envString', metavar='FILE', type='string', action='store',
                default   =   'xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/2DAlpha_env.tgz ./ ; tar -xvzf 2DAlpha_env.tgz ; cd CMSSW_7_4_7_patch2/src/ ; scram b ProjectRename ; scram b -j 10 ; cmsenv ; cd 2DAlphabet ;',
                dest      =   'envString',
                help      =   'Same argument given to `./development/runManySections.py --createCommandFile` to setup the condor environment')
parser.add_option('-u', '--unblind', action="store_true",
                default   =   False,
                dest      =   'unblind',
                help      =   'Only plot observed limit if true')
parser.add_option('-l', '--lumi', metavar='F', type='string', action='store',
                default       =       '35851',
                dest          =       'lumi',
                help          =       'Luminosity option')
parser.add_option('-H', '--hand', metavar='F', type='string', action='store',
                default       =       '',
                dest          =       'hand',
                help          =       'LH,RH,VL')

(options, args) = parser.parse_args()

# Open signal file
signal_file = open(options.signals,'r')
# Read in names as a list of strings and strip whitespace
signal_names = signal_file.readline().split(',')
signal_names = [n.strip() for n in signal_names]
# Read in mass as a list of strings, strip whitespace, and convert to ints
signal_mass = signal_file.readline().split(',')
signal_mass = [int(m.strip()) for m in signal_mass]
# Read in xsecs as a list of strings, strip whitespace, and convert to floats
signal_xsecs = signal_file.readline().split(',')
signal_xsecs = [float(x.strip()) for x in signal_xsecs]


# Open template
input_template_file = open(options.template,'r')

# Derive the tag
tag = options.template[options.template.find('input_')+len('input_'):options.template.find('_template.json')]

# If you want to do more than plot
if not options.plotOnly:

    # Try to create a folder for this <tag>
    try:
        # Make directory
        subprocess.call(['mkdir ' + tag], shell=True)
    except:
        print 'dir ' + tag + '/ already exists'
    
    # For each signal mass
    for name in signal_names:
        # Setup a signal specific config
        config_file_location = options.template.replace("template",name)
        config_file_name = config_file_location[config_file_location.find('input_'):config_file_location.find('.json')+len('.json')]
        print 'Running 2DAlphabet for ' + name

        # Execute some shell commands to setup 2DAlphabet for this config (notice we AREN'T running combine here)
        commands = ['sed s/' + options.dummy+'/'+name+'/g ' + options.template + ' > '+config_file_location,  # find and replace the dummy
                    'mkdir ' + tag + '/'+name,
                    'mv ' + config_file_location + ' ' + tag + '/'+name+'/' ,                             # move to folder <tag>/<tag>_signalname/
                    'python 2DAlphabet.py -i '+tag+'/'+name+'/'+config_file_name + ' -l -p -b']      # Run 2DAlphabet.py without the Combine step
        for c in commands:
            print '\t Executing ' + c
            subprocess.call([c],shell=True)

    # Now run combine with the listOfJobs generated
    subprocess.call(['source ' + tag+ '/listOfJobs.csh'],shell=True)
    
    for name in signal_names:
        subprocess.call(['mv '+'higgsCombine'+tag+'_'+name+'.Asymptotic.mH120.root ' + tag+'/'+name+'/'])

    # IN DEVELOPMENT
    # Now equipped with jobs and a list of them, we can submit to condor

    # # Edit the grid submission script
    # commands = ['''sed 's@ENVSTRING@"'''+options.envString+'''"@g' grid_sub_template.csh > grid_sub_'''+tag+'''.csh''', # using @ as sed delimiter
    #             "sed -i s@TEMPDIR@"+tag+"@g grid_sub_"+tag+".csh"]
    # for c in commands:
    #     print '\t Executing ' + c
    #     subprocess.call([c],shell=True)

    # # Sync the working directory to the tgz on EOS and submit
    # subprocess.call(["rsync2Denv"],shell=True)
    # subprocess.call(["source grid_sub_"+tag+".csh"],shell=True)

    # WaitForJobs(tag+'/listOfJobs.txt')


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
    this_output = TFile.Open(tag+'/'+this_name+'/'+'higgsCombine'+tag+'_'+this_name+'.Asymptotic.mH120.root')
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
            if options.unblind:
                y_limit.append(this_tree.limit*this_xsec)
            else:
                y_limit.append(0.0)
    
# Make Canvas and TGraphs (mostly stolen from other code that formats well)
climits = TCanvas("climits", "climits",700, 600)
climits.SetLogy(True)
climits.SetLeftMargin(.18)
climits.SetBottomMargin(.18)  

# NOT GENERIC
if options.hand == 'LH':
    cstr = 'L'
elif options.hand == 'RH':
    cstr = 'R'
elif options.hand == 'VL':
    cstr = 'LR'
else:
    cstr = ''


TPT = ROOT.TPaveText(.20, .22, .5, .27,"NDC")
TPT.AddText("All-Hadronic Channel") # NOT GENERIC
TPT.SetFillColor(0)
TPT.SetBorderSize(0)
TPT.SetTextAlign(12)

# Observed
g_limit = TGraph(len(x_mass), x_mass, y_limit)
g_limit.SetTitle("")
g_limit.SetMarkerStyle(0)
g_limit.SetMarkerColor(1)
g_limit.SetLineColor(1)
g_limit.SetLineWidth(3)
g_limit.SetMarkerSize(0.5) #0.5
g_limit.GetXaxis().SetTitle("M_{b*_{"+cstr+"}} (GeV)")  # NOT GENERIC
g_limit.GetYaxis().SetTitle("Upper Limit #sigma_{b*_{"+cstr+"}} #times B(b*_{"+cstr+"}#rightarrowtW) [pb]") # NOT GENERIC
g_limit.GetYaxis().SetRangeUser(0., 80.)
g_limit.GetXaxis().SetRangeUser(1, 3.2)
g_limit.SetMinimum(3.0e-3) #0.005
g_limit.SetMaximum(7000.)


# Expected
g_mclimit = TGraph(len(x_mass), x_mass, y_mclimit)
g_mclimit.SetTitle("")
g_mclimit.SetMarkerStyle(21)
g_mclimit.SetMarkerColor(1)
g_mclimit.SetLineColor(1)
g_mclimit.SetLineStyle(2)
g_mclimit.SetLineWidth(3)
g_mclimit.SetMarkerSize(0.)
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

g_limit.Draw('al')
g_error95.Draw("lf")
g_error.Draw("lf")
g_mclimit.Draw("l")
g_limit.Draw("l")
graphWP.Draw("l")
# graphWP.Draw("l")

# WPunc.Draw("lf")
# graphWPup.Draw("l")
# graphWPdown.Draw("l")


legend = TLegend(0.5, 0.45, 0.86, 0.84, '')
if options.unblind:
    legend.AddEntry(g_limit, "Observed", "l")
legend.AddEntry(g_mclimit, "Expected (95% CL)","l")
legend.AddEntry(g_error, "#pm 1 #sigma Expected", "f")
legend.AddEntry(g_error95, "#pm 2 #sigma Expected", "f")
legend.AddEntry(graphWP, "Theory b*_{"+cstr+"}", "l")   # NOT GENERIC
# legend.AddEntry(graphWPup, "Theory b*_{"+cstr+"} 1 #sigma uncertainty", "l")
g_limit.GetYaxis().SetTitleOffset(1.4)

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
climits.SaveAs(tag+"/limits_combine_"+options.lumi+"pb_"+options.signals[:options.signals.find('.')]+".pdf")

# Finally calculate the intercept
expectedLimit = Inter(g_mclimit,graphWP)[0]
upLimit = Inter(g_mcminus,graphWP)[0]
lowLimit = Inter(g_mcplus,graphWP)[0]

print str(expectedLimit/1000.) + ' +'+str(upLimit/1000.-expectedLimit/1000.) +' -'+str(expectedLimit/1000.-lowLimit/1000.) + ' TeV' # NOT GENERIC
