from TwoDAlphabetClass import TwoDAlphabet, runMLFit
import sys, traceback
from optparse import OptionParser
import subprocess
import header
import ROOT
from ROOT import *

parser = OptionParser()

parser.add_option('-q', '--tag', metavar='F', type='string', action='store',
                default =   '',
                dest    =   'quicktag',
                help    =   'Assigns a tag for this run')
parser.add_option('--rMin', metavar='F', type='string', action='store',
                default =   '0',
                dest    =   'rMin',
                help    =   'Minimum bound on r (signal strength)')
parser.add_option('--rMax', metavar='F', type='string', action='store',
                default =   '5',
                dest    =   'rMax',
                help    =   'Minimum bound on r (signal strength)')
parser.add_option("--recycleAll", action="store_true", 
                default =   False,
                dest    =   "recycleAll",
                help    =   "Recycle everything from the previous run with this tag")
parser.add_option("--skipFit", action="store_true", 
                default =   False,
                dest    =   "skipFit",
                help    =   "Skip fit and go directly to plotting (WARNING: Will use previous fit result if it exists and crash otherwise)")

(options, args) = parser.parse_args()

inputConfigs = args

print 'Setting on-fly parameters:'
print '\ttag\t\t = '+options.quicktag
print '\trecycleAll\t = '+str(options.recycleAll)
print '\tskipFit\t\t = '+str(options.skipFit)
print 'Remaining arguments:'
for i in inputConfigs:
    print '\t'+i

twoDinstances = []

# If simultaneous fit
if len(inputConfigs) > 1:
    # Instantiate all class instances
    for i in inputConfigs:
        instance = TwoDAlphabet(i,options.quicktag,options.recycleAll)
        twoDinstances.append(instance)

    # For each instance, check tags match and if they don't, ask the user for one
    for t in twoDinstances:
        if t.tag != twoDinstances[0].tag:
            print 'ERROR: tags in configuration files do not match. '+t.tag+' does not match '+twoDinstances[0].tag+'. Please make sure they match and try again. Quitting...'
            quit()
    thistag = twoDinstances[0].tag


    # Combine the cards
    print 'cd ' + thistag
    with header.cd(thistag):
        card_combination_command = 'combineCards.py'
        for i in twoDinstances:
            card_combination_command += ' '+i.name+'/card_'+i.name+'.txt'
        card_combination_command += ' > card_'+thistag+'.txt'

        print 'Executing ' + card_combination_command
        subprocess.call([card_combination_command],shell=True)
        for num in range(1,len(twoDinstances)+1):
            subprocess.call(["sed -i 's/ch"+str(num)+"_//g' card_"+thistag+".txt"],shell=True)

    if not options.skipFit:
        runMLFit(twoDinstances,options.rMin,options.rMax)

    # Plot
    for t in twoDinstances:
        try:
            t.plotFitResults('b',simfit=True)
        except Exception as exc:
            print traceback.format_exc()
            print exc
            print 'Failed to run b plots for '+t.name
        try:
            t.plotFitResults('s',simfit=True)
        except Exception as exc:
            print traceback.format_exc()
            print exc
            print 'Failed to run s plots for '+t.name

# If single fit
else:
    instance = TwoDAlphabet(inputConfigs[0],options.quicktag,options.recycleAll)
    
    if not skipFit:
        runMLFit([instance],options.rMin,options.rMax)
    thistag = instance.projPath

    # Plot
    try:
        instance.plotFitResults('b')
    except Exception as exc:
        print traceback.format_exc()
        print exc
        print 'Failed to run b plots for '+instance.name
    try:
        instance.plotFitResults('s')
    except Exception as exc:
        print traceback.format_exc()
        print exc
        print 'Failed to run s plots for '+instance.name
    