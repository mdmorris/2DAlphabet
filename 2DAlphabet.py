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
import input_organizer
import build_workspace
import plot_fit_results
import header
from header import ascii_encode_dict


#########################################################
#                       Options                         #
#########################################################
parser = OptionParser()

parser.add_option('-i', '--input', metavar='FILE', type='string', action='store',
                  default   =   '',
                  dest      =   'input',
                  help      =   'JSON file to be imported. Name should be "input_<tag>.json" where tag will be used to organize outputs')
parser.add_option('-m', '--move', action="store_true",
                  default   =   False,
                  dest      =   'move',
                  help      =   'Moves files to folder <tag> if True')


(options, args) = parser.parse_args()

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

tag = options.input[options.input.find('_')+1:options.input.find('.')]

if options.move:
    try:
        subprocess.call(['mkdir ' + tag], shell=True)
    except:
        print 'dir ' + tag + '/ already exists'

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
                    if subkey.find(glob_var) != 1:                  # check subkey
                        subkey = subkey.replace(glob_var,input_config['GLOBAL'][glob_var])   # replace it

                    if subkey != 'HELP':                                        # If the subkey isn't HELP, the key value is a dict
                        for subsubkey in input_config[mainkey][subkey].keys():  # so loop though subsubkeys
                            if subsubkey.find(glob_var) != 1:                                   # check subsubkey
                                subsubkey = subsubkey.replace(glob_var,input_config['GLOBAL'][glob_var])    # replace it

                            try:
                                if input_config[mainkey][subkey][subsubkey].find(glob_var) != -1:                               # check subsubkey val
                                    input_config[mainkey][subkey][subsubkey] = input_config[mainkey][subkey][subsubkey].replace(glob_var,input_config['GLOBAL'][glob_var]) # replace it

                            except:                 # Exception if not of type string
                                continue

# Save out the json with the GLOBAL var replacement
fInput_config_vars_replaced = open('input_'+tag+'_vars_replaced.json', 'w')
json.dump(input_config,fInput_config_vars_replaced,indent=2, sort_keys=True)
fInput_config_vars_replaced.close()
if options.move:
    subprocess.call(['mv ' +'input_'+tag+'_vars_replaced.json ' + tag + '/'], shell=True)
print 'Done'


# A quick flag to check for blinding
if 'SIGSTART' in input_config['BINNING']['X'] and 'SIGEND' in input_config['BINNING']['X']:
    blinded = True
elif 'SIGSTART' in input_config['BINNING']['X'] and 'SIGEND' not in input_config['BINNING']['X']:
    print 'Error! SIGSTART is defined but SIGEND is not. Quitting.'
    quit()
elif 'SIGSTART' not in input_config['BINNING']['X'] and 'SIGEND' in input_config['BINNING']['X']:
    print 'Error! SIGEND is defined but SIGSTART is not. Quitting.'
    quit()
else:
    blinded = False

#########################################################
#               Organize input histograms               #
#########################################################
# input_organizer.main() returns a dictionary with all of the TH2s organized
print 'Organizing histograms into single file...'
organized_dict = input_organizer.main(input_config,blinded)
print 'Done'

#########################################################
#             Make the workspace for Combine            #
#########################################################
# Make the RooWorkspace - creates workspace name 'w_2D' in base.root
workspace = build_workspace.main(organized_dict,input_config,blinded,tag)
print 'Workspace built'

#########################################################
#             Make the data card for Combine            #
#########################################################
# Make the data card - makes a text file named card_2D.txt, return 0
print 'Making Combine card...'
make_card.main(input_config, blinded, tag)
print 'Done'

# Run Combine
print 'Executing combine -M MaxLikelihoodFit card_'+tag+'.txt --saveWithUncertainties --saveWorkspace -v 2 --rMin -50 --rMax 50' 
subprocess.call(['combine -M MaxLikelihoodFit card_'+tag+'.txt --saveWithUncertainties --saveWorkspace -v 2 --rMin -50 --rMax 50'], shell=True)

# Test that Combine ran successfully 
if not os.path.isfile('MaxLikelihoodFitResult.root'):
    print "Combine failed and never made MaxLikelihoodFitResult.root. Quitting..."
    quit()

if options.move:
    subprocess.call(['mv base_'+tag+'.root ' + tag + '/'], shell=True)
    subprocess.call(['mv card_'+tag+'.txt ' + tag + '/'], shell=True)

# Plot
plot_fit_results.main(input_config,organized_dict,blinded,tag)

