# Requires that you've run a b-only fit via run_MLfit.py at least once! Need to use the output.

from TwoDAlphabetClass import TwoDAlphabet, runLimit
import sys, traceback, os
from optparse import OptionParser
import subprocess
import header

parser = OptionParser()

parser.add_option('-q', '--tag', metavar='F', type='string', action='store',
                default =   '',
                dest    =   'quicktag',
                help    =   'Assigns a tag for this run')
parser.add_option('-d', '--projDir', metavar='F', type='string', action='store',
                default =   '',
                dest    =   'projDir',
                help    =   'Points to the directory where the b-only fit result is located')
parser.add_option("--unblindData", action="store_true", 
                default =   False,
                dest    =   "unblindData",
                help    =   "Unblind the observation and calculate the observed limit")


(options, args) = parser.parse_args()

inputConfigsAndArgs = args

# Assign and summarize input
stringSwaps = {}
inputConfigs = []
postfitWorkspaceDir = options.projDir
for i,c in enumerate(inputConfigsAndArgs):
    if ':' in c: # specify string swap
        stringSwaps[c.split(':')[0]] = c.split(':')[1]
        print c.split(':')[0] +' = '+c.split(':')[1]
    elif '.json' in c:
        inputConfigs.append(c)


print 'tag                              = ' + str(options.quicktag)
print 'Location of b-only fit workspace = ' + postfitWorkspaceDir 
print 'Unblind data points                = '+ str(options.unblindData)
print 'Config Replacements:'
for s in stringSwaps.keys():
    print '\t'+ s + ' -> ' + stringSwaps[s]
print 'Configs:'
for c in inputConfigs:
    print '\t'+c
    
# Check for workspace to load
if os.path.isfile(postfitWorkspaceDir +'/fitDiagnostics.root'):
    print 'Loading '+postfitWorkspaceDir +'/fitDiagnostics.root'
else:
    print "ERROR: post-fit '"+postfitWorkspaceDir+"/fitDiagnostics.root' could not be found. Please specify the location of the post-fit workspace (postfitb_workspace.root) that you would like to load. Quitting..."
    quit()


# Initialize
twoDinstances = []

# If simultaneous fit
if len(inputConfigs) > 1:
    # Instantiate all class instances
    for i in inputConfigs:
        instance = TwoDAlphabet(i,options.quicktag,stringSwaps=stringSwaps)
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

    if 'condor' in os.getcwd(): runLimit(twoDinstances,postfitWorkspaceDir,blindData=(not bool(options.unblindData)),location='condor')
    else: runLimit(twoDinstances,postfitWorkspaceDir,blindData=(not bool(options.unblindData)),location='local')

# If single fit
else:
    instance = TwoDAlphabet(inputConfigs[0],options.quicktag,stringSwaps=stringSwaps)
    if 'condor' in os.getcwd(): runLimit([instance],postfitWorkspaceDir,blindData=(not bool(options.unblindData)),location='condor')
    else: runLimit([instance],postfitWorkspaceDir,blindData=(not bool(options.unblindData)),location='local')
