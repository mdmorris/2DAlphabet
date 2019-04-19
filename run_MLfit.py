from TwoDAlphabetClass import TwoDAlphabet, runMLFit
import sys, traceback
import subprocess
import header
import ROOT
from ROOT import *

inputConfigs = sys.argv[1:]
quicktag = False
for i,c in enumerate(inputConfigs):
    if 'tag=' in c:
        quicktag = c.split('=')[1]
        inputConfigs.pop(i)

twoDinstances = []

# If simultaneous fit
if len(inputConfigs) > 1:
    # Instantiate all class instances
    for i in inputConfigs:
        instance = TwoDAlphabet(i,quicktag)
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

    runMLFit(twoDinstances)

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
    instance = TwoDAlphabet(inputConfigs[0],quicktag)
    runMLFit([instance])
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
    