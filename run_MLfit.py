from TwoDAlphabetClass import TwoDAlphabet, runMLFit, runLimit
import sys
import subprocess
import header

inputConfigs = sys.argv[1:]

twoDinstances = []

# If simultaneous fit
if len(inputConfigs) > 1:
    # Instantiate all class instances
    for i in inputConfigs:
        instance = TwoDAlphabet(i)
        twoDinstances.append(instance)

    # For each instance, check tags match and if they don't, ask the user for one
    quicktag = False
    for t in twoDinstances:
        if t.tag != twoDinstances[0].tag:
            print 'ERROR: tags in configuration files do not match. '+t.tag+' does not match '+twoDinstances[0].tag+'. Please make sure they match and try again. Quitting...'
            quit()

    # Combine the cards
    print 'cd ' + twoDinstances[0].tag
    with header.cd(twoDinstances[0].tag):
        card_combination_command = 'combineCards.py'
        for i in twoDinstances:
            card_combination_command += ' '+i.name+'/card_'+i.name+'.txt'
        card_combination_command += ' > card_'+twoDinstances[0].tag+'.txt'

        print 'Executing ' + card_combination_command
        subprocess.call([card_combination_command],shell=True)
        for num in range(1,len(twoDinstances)+1):
            subprocess.call(["sed -i 's/ch"+str(num)+"_//g' card_"+twoDinstances[0].tag+".txt"],shell=True)

    runMLFit(twoDinstances)

    for t in twoDinstances:
        try:
            t.plotFitResults('b',simfit=True)
        except:
            print 'Failed to run b plots for '+t.name
        try:
            t.plotFitResults('s',simfit=True)
        except:
            print 'Failed to run s plots for '+t.name

# If single fit
else:
    instance = TwoDAlphabet(inputConfigs[0])
    runMLFit([instance])
    instance.plotFitResults('b')
    instance.plotFitResults('s')
    