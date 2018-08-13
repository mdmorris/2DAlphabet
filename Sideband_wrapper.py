#####################################################################################################
# Sideband_wrapper.py - Written by Lucas Corcodilos, 7/6/18                                         #
# ---------------------------------------------------------                                         #
# This wrapper is designed to run 2DAlphabet.py without the Combine fit for a standard config       #
# and the config for background enriched sidebands (thus creating the workspaces for the signal     #
# region and the sidebands). It then runs a simultaneous fit after the cards have been combined.    #
#                                                                                                   #
# This method allows the template morphing fit of a background in an enriched sideband to constrain #
# the nusiance parameters of the fit and thus the background's shape in the signal region.          #
#####################################################################################################

from optparse import OptionParser
import subprocess
import json
import pprint
pp = pprint.PrettyPrinter(indent = 2)

import ROOT
from ROOT import *

import plot_fit_results
import input_organizer
import header
from header import ascii_encode_dict

gStyle.SetOptStat(0)
gROOT.SetBatch(kTRUE)

parser = OptionParser()

parser.add_option('-i', '--input', metavar='FILE', type='string', action='store',
                  default   =   '',
                  dest      =   'input',
                  help      =   'JSON file to be imported. Name should be "input_<tag>.json" where tag will be used to organize outputs')
parser.add_option('-s', '--sidebandInputs', metavar='FILE', type='string', action='store',
                  default   =   '',
                  dest      =   'sidebandInputs',   # Don't have file names with more than two underscores!
                  help      =   'Comma separated list of JSON files. Name should be "input_<tag>_<bkg>.json" where tag will be used to organize outputs and bkg corresponds to the enriched background process in the config file')
parser.add_option('-b', '--blinded', action="store_true",
                  default   =   False,
                  dest      =   'blinded',
                  help      =   'Turns off data points. Different from blinding the signal region during the fit')
parser.add_option('-p', '--plotOnly', action="store_true",
                  default   =   False,
                  dest      =   'plotOnly',
                  help      =   'Only plots. Does not run anything else (assumes that has already been done')


(options, args) = parser.parse_args()

# Basic exiting if the inputs are not formatted correctly
if options.sidebandInputs == '':
    print 'No sideband input configuration files given. Consider using just 2DAlphabet.py - Quitting...'
    quit()
elif len(options.sidebandInputs.split('_')) != 3:
    print 'Input sideband configuration file names are not of the correct format. They should be of the form "input_<tag>_<bkg>.json" where tag matches that given in the main configuration file and bkg corresponds to the enriched background process in the config file'
    quit()

# Get tag and bkg names
tag = options.input.split('_')[1][:options.input.split('_')[1].find('.')]

sideband_cfg_strings = options.sidebandInputs.split(',')
sideband_dirs = {}          
for cfg in sideband_cfg_strings:
    bkgName = cfg[cfg.find(tag)+len(tag)+1:cfg.find('.')]
    sideband_dirs[bkgName] = tag+'/'+bkgName+'/'

if not options.plotOnly:
    # Run 2D Alphabet setup for combine
    combine_calls = ['python 2DAlphabet.py -i '+options.input+' -p']    # Only call input and to run pseudo2D first

    for cfg in sideband_cfg_strings:
        combine_calls.append('python 2DAlphabet.py -i '+cfg+' --batch -p')


    for c in combine_calls:
        print 'Executing ' + c
        subprocess.call([c],shell=True)


    # Combine the cards
    commands = ['mv '+tag+'/card_'+tag+'.txt ./']
    card_call = 'combineCards.py card_'+tag+'.txt '                 # Main card
    for bkgName in sideband_dirs.keys():
        commands.append('mv '+sideband_dirs[bkgName]+'card_'+tag+'_'+bkgName+'.txt ./')
        card_call+='card_'+tag+'_'+bkgName+'.txt '   # Each sideband card
    card_call+='> card_master.txt'                                          # Send it to master

    commands.append(card_call)

    # Need to make a change to card_master.txt. combineCards.py will rename the categories with prefix ch#_
    # where the # is an integer (starting with 1) and corresponds to the order of the inputs to combineCards.py
    # The problem with this is that it's totally non-descriptive and there's no easy way to know what you're
    # getting if you don't know the order. So I've renamed the sideband categories like pass_<bkgName> and will
    # remove all instances of ch#_
    for num in range(1,len(sideband_cfg_strings)+2):
        commands.append("sed -i 's/ch"+str(num)+"_//g' card_master.txt")

    commands.append('mv card_'+tag+'.txt '+tag+'/')
    for bkgName in sideband_dirs.keys():
        commands.append('mv card_'+tag+'_'+bkgName+'.txt '+tag+'/')

    commands.append('mv card_master.txt '+tag+'/')

    for c in commands:
        print 'Executing ' + c
        subprocess.call([c],shell=True)

    # Run Combine
    print 'Executing ' + 'combine -M MaxLikelihoodFit '+ tag+'/card_master.txt --saveWithUncertainties --saveWorkspace --rMin -5 --rMax 5'
    subprocess.call(['combine -M MaxLikelihoodFit '+tag + '/card_master.txt --saveWithUncertainties --saveWorkspace --rMin -5 --rMax 5'],shell=True)

    subprocess.call(['mv MaxLikelihoodFitResult.root '+tag+'/'],shell=True)
    subprocess.call(['mv mlfit.root '+tag+'/'],shell=True)
    subprocess.call(['mv higgsCombineTest.MaxLikelihoodFit.mH120.root '+tag+'/'],shell=True)

########
# Plot #
########

# Some requirements first
# JSON -> Dictionary
with open(tag+'/input_'+tag+'_new_fit_guesses.json') as fInput_config:
    input_config = json.load(fInput_config, object_hook=ascii_encode_dict)  # Converts most of the unicode to ascii

    for process in [proc for proc in input_config['PROCESS'].keys() if proc != 'HELP']:
        for index,item in enumerate(input_config['PROCESS'][process]['SYSTEMATICS']):           # There's one list that also
            input_config['PROCESS'][process]['SYSTEMATICS'][index] = item.encode('ascii')       # needs its items converted

# Make organized Dict
organized_dict = input_organizer.main(input_config,True)

# Plot
plot_fit_results.main(input_config,organized_dict,True,tag,tag+'/',options.blinded)



# Now do the same but for the sidebands
for string in sideband_cfg_strings:
    bkgName = cfg[cfg.find(tag)+len(tag)+1:cfg.find('.')]
    with open(tag+'/'+bkgName+'/input_'+tag+'_'+bkgName+'_new_fit_guesses.json') as fInput_config:
        this_input_config = json.load(fInput_config, object_hook=ascii_encode_dict)  # Converts most of the unicode to ascii

        for process in [proc for proc in this_input_config['PROCESS'].keys() if proc != 'HELP']:
            for index,item in enumerate(this_input_config['PROCESS'][process]['SYSTEMATICS']):           # There's one list that also
                this_input_config['PROCESS'][process]['SYSTEMATICS'][index] = item.encode('ascii')       # needs its items converted

    # Make organized Dict
    this_organized_dict = input_organizer.main(this_input_config,True)

    # Plot
    plot_fit_results.main(this_input_config,this_organized_dict,True,tag,tag+'/'+bkgName+'/',options.blinded,suffix='_'+bkgName,mlfitresultDir=tag+'/')
