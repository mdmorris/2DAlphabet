#####################################################################################################
# 2DAlphabet.py - written by Lucas Corcodilos, 3/7/18                                               #
# ---------------------------------------------------                                               #
# This is the wrapper and callable script that runs the full 2DAlphabet workflow. The only inputs   #
# are a properly formatted JSON file (see input_example.txt for an example) and the options         #
# which specify which parts of the workflow to run (default is all).                                #
#####################################################################################################

#########################################################
#                       Imports                         #
#########################################################
from optparse import OptionParser
import json
import subprocess
import os
import pprint
pp = pprint.PrettyPrinter(indent = 2)

import make_card
import get_fit_guesses
import input_organizer
import build_fit_workspace
import plot_fit_results
import header
from header import ascii_encode_dict

import ROOT
from ROOT import *

gStyle.SetOptStat(0)

#########################################################
#                       Options                         #
#########################################################
parser = OptionParser()

parser.add_option('-i', '--input', metavar='FILE', type='string', action='store',
                  default   =   '',
                  dest      =   'input',
                  help      =   'JSON file to be imported. Name should be "input_<tag>.json" where tag will be used to organize outputs')
parser.add_option('-p', '--pseudo2D', action="store_true",
                  default   =   False,
                  dest      =   'pseudo2D',
                  help      =   'Recalculate the fit guesses using pseudo2D method')
parser.add_option('-P', '--plotOnly', action="store_true",
                  default   =   False,
                  dest      =   'plotOnly',
                  help      =   'Only plots if True')
parser.add_option('-d', '--draw', action="store_true",
                  default   =   False,
                  dest      =   'draw',
                  help      =   'Draws canvases live')
parser.add_option('-s', '--signalOff', action="store_true",
                  default   =   False,
                  dest      =   'signalOff',
                  help      =   'Turns off signal by using --expectSignal=0 option in Combine')
parser.add_option('-l', '--runLimits', action="store_true",
                  default   =   False,
                  dest      =   'runLimits',
                  help      =   'Runs Combine limits and plots outputs')
parser.add_option('-f', '--runFit', action="store_true",
                  default   =   False,
                  dest      =   'runFit',
                  help      =   'Runs Combine Rp/f fit and plots outputs')
parser.add_option('-b', '--batch', action="store_true",
                  default   =   False,
                  dest      =   'batch',
                  help      =   'Runs limits in batch mode for multiple signals')
# parser.add_option('-D', '--globalDir', type='string', action="store",
#                   default   =   '',
#                   dest      =   'globalDir',
#                   help      =   'Only used by the do_full_limits.py wrapper to do multiple signals')


(options, args) = parser.parse_args()

if options.draw == False:
    gROOT.SetBatch(kTRUE)

if options.runLimits and options.runFit:
    print 'ERROR: Cannot run limits and fit at same time. Quitting...'
    quit()

#########################################################
#                     JSON Processing                   #
#########################################################
with open(options.input) as fInput_config:
    # try:
        input_config = json.load(fInput_config, object_hook=ascii_encode_dict)  # Converts most of the unicode to ascii

        for process in [proc for proc in input_config['PROCESS'].keys() if proc != 'HELP']:
            for index,item in enumerate(input_config['PROCESS'][process]['SYSTEMATICS']):           # There's one list that also
                input_config['PROCESS'][process]['SYSTEMATICS'][index] = item.encode('ascii')       # needs its items converted

    # except:
    #     print "Input configuration is not in JSON format. Quitting..."
    #     quit()

sig_option = ' --rMin -5 --rMax 5'
sig_tag = ''
if options.signalOff:
    sig_option = '  --rMin 0 --rMax 0'
    sig_tag = '_nosig'

# If input has form input_<tag>_<sig>.json...
if len(options.input.split('_')) == 3:
    maindir = options.input.split('_')[1] + '/'
    subdir = options.input.split('_')[2]
    subdir = subdir[:subdir.find('.')] + '/'
    tag = maindir[:-1] + '_' + subdir[:-1] + sig_tag
    try:
        subprocess.call(['mkdir ' + maindir  + subdir[:-1]+sig_tag], shell=True)
        subprocess.call(['mkdir ' + maindir + subdir[:-1]+sig_tag + '/plots'], shell=True)
    except:
        print 'dir ' + maindir  + subdir[:-1]+sig_tag + '/ already exists'

# Otherwise treat it as input_<tag>.json
else:
    maindir = options.input[options.input.find('input_')+len('input_'):options.input.find('.')] + '/'
    tag = maindir[:-1] + sig_tag
    subdir = ''
    try:
        subprocess.call(['mkdir ' + tag], shell=True)
        subprocess.call(['mkdir ' + tag + '/plots'], shell=True)
    except:
        print 'dir ' + tag + '/ already exists'

print 'maindir: ' + maindir
print 'subdir: ' + subdir

#####################################
# Do GLOBAL variable substitution   #
# --------------------------------- #
# Relies on certain JSON structure. #
# Anything marked as 'changeable'   #
# needs to be checked for GLOBAL    #
# variable.                         # 
#                                   #
# - HELP is unchangeable            #
# - mainkeys are unchangeable       #
# - subkeys are changeable in       #
#   PROCESS, SYSTEMATIC, and FIT    #
# - subsubkeys are changeable       #
# - subsubkey values are changeable #
#                                   #
# CURRENTLY ONLY SUPPORTS STRING    #
# REPLACEMENT                       #
#####################################

print "Doing GLOBAL variable replacement in input json...",
for glob_var in input_config['GLOBAL'].keys():
    if glob_var != "HELP":                                          # For each key in GLOBAL that is not HELP

        for mainkey in input_config.keys():                         # For each main (top level) key in input_config that isn't GLOBAL
            if mainkey != 'GLOBAL':                                 # Mainkeys are all unchangeable (uppercase) so no check performed
                
                for subkey in input_config[mainkey].keys():         # For each subkey of main key dictionary
                    if subkey.find(glob_var) != -1:                  # check subkey
                        subkey = subkey.replace(glob_var,input_config['GLOBAL'][glob_var])   # replace it

                    if subkey == 'HELP':                                        # If the subkey isn't HELP, the key value is a dict
                        continue
                    elif subkey.find('FORM') != -1:
                        if subkey.find(glob_var) != -1:                  # check subkey
                            subkey = subkey.replace(glob_var,input_config['GLOBAL'][glob_var])   # replace it
                    else:
                        for subsubkey in input_config[mainkey][subkey].keys():  # so loop though subsubkeys
                            if subsubkey.find(glob_var) != -1:                                   # check subsubkey
                                subsubkey = subsubkey.replace(glob_var,input_config['GLOBAL'][glob_var])    # replace it

                            try:
                                if input_config[mainkey][subkey][subsubkey].find(glob_var) != -1:                               # check subsubkey val
                                    input_config[mainkey][subkey][subsubkey] = input_config[mainkey][subkey][subsubkey].replace(glob_var,input_config['GLOBAL'][glob_var]) # replace it

                            except:                 # Exception if not of type string
                                continue

# Save out the json with the GLOBAL var replacement
fInput_config_vars_replaced = open(maindir + subdir + '/input_'+tag+'_vars_replaced.json', 'w')
json.dump(input_config,fInput_config_vars_replaced,indent=2, sort_keys=True)
fInput_config_vars_replaced.close()

print 'Done'

# A quick flag to check for blinding
if input_config['BINNING']['X']['BLINDED'] == True:
    if options.runLimits == False:
        print 'Background estimate is blinded'
        blinded = True
    else:
        print 'Background estimate is NOT blinded'
        blinded = False
else:
    print 'Background estimate is NOT blinded'
    blinded = False



#########################################################
#             Get new fit parameter guesses             #
#########################################################
if options.pseudo2D == True:
    input_config = get_fit_guesses.main(input_config,blinded, maindir+subdir)
    fInput_config_new_fit_guesses = open(maindir+subdir + 'input_'+tag+'_new_fit_guesses.json', 'w')
    json.dump(input_config,fInput_config_new_fit_guesses,indent=2, sort_keys=True)
    fInput_config_new_fit_guesses.close()

#########################################################
#               Organize input histograms               #
#########################################################
# input_organizer.main() returns a dictionary with all of the TH2s organized
print 'Organizing histograms into single file...'
organized_dict = input_organizer.main(input_config,blinded)
print 'Done'

# Plot only
if options.plotOnly:
    plot_fit_results.main(input_config,organized_dict,blinded,tag,maindir+subdir)
    quit()


#########################################################
#             Make the workspace for Combine            #
#########################################################
# Make the RooWorkspace - creates workspace name 'w_2D' in base.root
workspace = build_fit_workspace.main(organized_dict,input_config,blinded,tag)
subprocess.call(['mv base_'+tag+'.root ' + maindir + subdir+'/'], shell=True)
print 'Workspace built'

#########################################################
#             Make the data card for Combine            #
#########################################################
# Make the data card - makes a text file named card_2D.txt, return 0
print 'Making Combine card...'
make_card.main(input_config, blinded, tag, maindir+subdir)
print 'Done'

syst_option = ''
# Check if signal has systematics
for proc in input_config['PROCESS'].keys():
    if type(input_config['PROCESS'][proc]) == dict:
        if input_config['PROCESS'][proc]['CODE'] == 0:
            if len(input_config['PROCESS'][proc]['SYSTEMATICS']) == 0:
                syst_option = ' -S 0'

if options.runFit:
    
    # Run Combine
    print 'Executing combine -M MaxLikelihoodFit '+maindir+subdir + '/card_'+tag+'.txt --saveWithUncertainties --saveWorkspace' + syst_option + sig_option 
    subprocess.call(['combine -M MaxLikelihoodFit '+maindir+subdir + '/card_'+tag+'.txt --saveWithUncertainties --saveWorkspace' + syst_option + sig_option], shell=True)

    # Test that Combine ran successfully 
    if not os.path.isfile('MaxLikelihoodFitResult.root'):
        print "Combine failed and never made MaxLikelihoodFitResult.root. Quitting..."
        quit()

    subprocess.call(['mv MaxLikelihoodFitResult.root ' + maindir+subdir + '/'], shell=True)
    subprocess.call(['mv higgsCombineTest.MaxLikelihoodFit.mH120.root ' + maindir+subdir + '/'], shell=True)
    subprocess.call(['mv mlfit.root ' + maindir+subdir + '/'], shell=True)

    plot_fit_results.main(input_config,organized_dict,blinded,tag,maindir+subdir)



if options.runLimits:
    # If running with the wrapper that processes multiple signals (for limits), we'll write out the command for the wrapper script
    if options.batch:
        limit_command = open(maindir+'listOfJobs.csh','a+')
        limit_command.write('combine -M Asymptotic '+maindir+subdir +'card_'+tag+'.txt --saveWorkspace --name '+tag+ ' ' + syst_option + sig_option+'\n')
        limit_command.close()
    
    else:
        print 'Executing combine -M Asymptotic '+maindir+subdir +'/card_'+tag+'.txt --saveWorkspace --name '+tag+ ' ' + syst_option + sig_option 
        subprocess.call(['combine -M Asymptotic '+maindir+subdir +'/card_'+tag+'.txt --saveWorkspace --name '+tag+ ' ' + syst_option + sig_option], shell=True)
        subprocess.call(['mv higgsCombine'+tag+'.Asymptotic.mH120.root ' + maindir+subdir +'/'], shell=True)


# subprocess.call(['mv base_'+tag+'.root ' + tag + '/'], shell=True)
# subprocess.call(['mv card_'+tag+'.txt ' + tag + '/'], shell=True)