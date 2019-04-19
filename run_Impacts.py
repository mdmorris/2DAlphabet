import sys, os
import subprocess
import header

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--projDir", dest="projDir",
                  help="Home of the project - has the cards, fit results, etc")
parser.add_option("-c", "--crab",
                  action="store_true", dest="crab", default=False,
                  help="Turn crab grid submission on")
parser.add_option("-p", "--post",
                  action="store_true", dest="post", default=False,
                  help="Run the post processing to get impact plot")
parser.add_option("-s", "--storage",
                  dest="storage", default='T3_US_FNALLPC',
                  help="Crab3 storage site (config.Site.storageSite)")

(options, args) = parser.parse_args()

projDir = options.projDir # home of the workspace - has the cards, fit results, etc
taskName = projDir.split('/')[0]
if projDir.split('/')[-1] != '': card_tag = projDir.split('/')[-1]
else: card_tag = projDir.split('/')[-2]

if taskName == '':
    print 'ERROR in project directory name (where your workspace and data card lives). Did you accidentally provide a leading slash? (ie /projDir/) Quitting...'
    quit()
if options.crab:
    print 'Crab task name = '+taskName

if not os.path.isdir(projDir): 
    print projDir +' is not a directory. Quitting...'
    quit()

# By default, this calculates the impacts for every RooRealVar in your workspace.
# That would mean EVERY FAIL BIN would need to be scanned (100s or even 1000s of parameters).
# So instead, we'll be a list of only the nuisance parameters and ask to just fit those

print 'cd '+projDir
with header.cd(projDir):
    twoDAlphaConfig = header.openJSON('rerunConfig.json')
    impactNuisanceString = '--named '
    for s in twoDAlphaConfig['SYSTEMATIC'].keys():
        if s !='HELP':
            impactNuisanceString+=s+','
    # Cut off the trailing comma
    impactNuisanceString = impactNuisanceString[:-1]

    commands = []

    if not options.post:
        # Remove old runs if they exist
        # print 'rm *_paramFit_*.root *_initialFit_*.root'
        # subprocess.call(['rm *_paramFit_*.root *_initialFit_*.root'],shell=True)
        header.executeCmd('text2workspace.py -b card_'+card_tag+'.txt -o impactworkspace.root')
        initialfit_cmd = 'combineTool.py -M Impacts -n '+taskName+' -d impactworkspace.root --doInitialFit --robustFit 1 -m 120 '+impactNuisanceString
        header.executeCmd(initialfit_cmd)

        if options.crab:
            # Need to write a custom crab config for the storage site
            crabConfig = open('custom_crab.py','w')
            crabConfig.write('def custom_crab(config):')
            print '>> Customising the crab config'
            crabConfig.write("    config.Site.storageSite = '"+options.storage+"'")
            crabConfig.close()

            print 'Executing dry-run of crab jobs'
            impact_cmd = 'combineTool.py -M Impacts -n '+taskName+' -d impactworkspace.root --robustFit 1 --doFits -m 120 --job-mode crab3 --task-name Impacts'+taskName+' --custom-crab custom_crab.py '+impactNuisanceString
            subprocess.call([impact_cmd+' --dry-run'],shell=True)
            proceed = raw_input('Please examine this command and confirm it is correct before submitting to crab. Do you wish to proceed? [Y/N]')
            if proceed.lower() == 'y':
                commands.append('source /cvmfs/cms.cern.ch/crab3/crab.sh; cmsenv')
                commands.append(impact_cmd)
                
            else:
                print 'Quitting...'
                quit()

        else:     
            commands.append('combineTool.py -M Impacts -n '+taskName+' -d impactworkspace.root --robustFit 1 --doFits -m 120 '+impactNuisanceString)

    elif options.post:
        # Grab the crab output, untar it, and put it in the main directory (now it looks like it was run locally)
        if options.crab:
            commands.append('crab getoutput -d crab_Impacts'+taskName)
            for filename in os.listdir('crab_Impacts'+taskName+'/results/'):
                commands.append('tar -xvf crab_ImpactsunitTest/results/'+filename)

        commands.append('combineTool.py -M Impacts -n '+taskName+' -d impactworkspace.root -m 120 '+impactNuisanceString+' -o impacts.json')
        commands.append('plotImpacts.py -i impacts.json -o impacts')

    # Run commands
    for c in commands:
        print 'Executing: '+c
        subprocess.call([c],shell=True)