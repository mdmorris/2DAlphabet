import sys, os
import subprocess
import header
import CombineHarvester.CombineTools.ch as ch

def SystematicParser(cardname):
    systs = []
    f = open(cardname,'r')
    for l in f.readlines():
        if 'lnN' in l or 'shape' in l or 'rpf' in l:
            syst_name = l.split(' ')[0]
            if syst_name != 'shapes': systs.append(syst_name.rstrip())

    return systs

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--projDir", dest="projDir",
                  help="Home of the project - has the cards, fit results, etc")
parser.add_option("--condor",
                  action="store_true", dest="condor", default=False,
                  help="Turn condor grid submission on")
parser.add_option("-p", "--post",
                  action="store_true", dest="post", default=False,
                  help="Run the post processing to get impact plot")
# parser.add_option("-s", "--storage",
#                   dest="storage", default='T3_US_FNALLPC',
#                   help="Crab3 storage site (config.Site.storageSite)")

(options, args) = parser.parse_args()

projDir = options.projDir # home of the workspace - has the cards, fit results, etc
taskName = projDir.split('/')[0]
if projDir.split('/')[-1] != '': card_tag = projDir.split('/')[-1]
else: card_tag = projDir.split('/')[-2]

if taskName == '':
    raise NameError('ERROR in project directory name (where your workspace and data card lives). Did you accidentally provide a leading slash? (ie /projDir/) Quitting...')
if options.condor:
    print 'Condor task name = '+taskName

if not os.path.isdir(projDir): 
    raise TypeError(projDir +' is not a directory. Quitting...')

# By default, this calculates the impacts for every RooRealVar in your workspace.
# That would mean EVERY FAIL BIN would need to be scanned (100s or even 1000s of parameters).
# So instead, we'll be a list of only the nuisance parameters and ask to just fit those

if options.condor: 
    header.executeCmd('tar --exclude="*.tgz" --exclude="*.std*" --exclude="run_combine*.sh" --exclude="*GoodnessOfFit*" --exclude="*.png" --exclude="*.pdf" --exclude="*.log" -czvf tarball.tgz '+projDir)
    header.executeCmd('mv tarball.tgz '+projDir)
print 'cd '+projDir
with header.cd(projDir):
    systs = SystematicParser('card_'+card_tag+'.txt')
    impactNuisanceString = '--named '
    for s in systs:
        impactNuisanceString+=s+','
    impactNuisanceString = impactNuisanceString[:-1]# Cut off the trailing comma

    if not options.post:
        # Remove old runs if they exist
        header.executeCmd('rm *_paramFit_*.root *_initialFit_*.root')
        # Build a post-fit workspace
        #header.executeCmd('text2workspace.py -b card_'+card_tag+'.txt -o impactworkspace.root')
        header.setSnapshot()
        initialfit_cmd = 'combineTool.py -M Impacts -n '+taskName+' --rMin -15 --rMax 15 -d initialFitWorkspace.root --snapshotName initialFit --doInitialFit --cminDefaultMinimizerStrategy 0 -m 2000 '+impactNuisanceString
        header.executeCmd(initialfit_cmd)
        impact_cmd = 'combineTool.py -M Impacts -n '+taskName+' --rMin -15 --rMax 15 -d initialFitWorkspace.root --snapshotName initialFit --doFits --cminDefaultMinimizerStrategy 0 -m 2000 '+impactNuisanceString
        if options.condor:
            JOB_PREFIX = """#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/10XwithNano.tgz ./
export SCRAM_ARCH=slc6_amd64_gcc700
scramv1 project CMSSW CMSSW_10_2_13
tar -xzf 10XwithNano.tgz
rm 10XwithNano.tgz

mkdir tardir; cp tarball.tgz tardir/; cd tardir
tar -xzf tarball.tgz
cp -r * ../CMSSW_10_2_13/src/2DAlphabet/
cd ../CMSSW_10_2_13/src/2DAlphabet/
eval `scramv1 runtime -sh`
scramv1 b clean; scramv1 b

cd %s
                """ % (projDir)
            job_prefix_out = open('impact_prefix.txt','w')
            job_prefix_out.write(JOB_PREFIX)
            job_prefix_out.close()

            impact_cmd = impact_cmd+' --job-mode condor --dry-run --prefix-file impact_prefix.txt --sub-opts "transfer_input_files = tarball.tgz" --task-name Impacts'+taskName
        else:
            header.executeCmd(impact_cmd)

    elif options.post:
        # Grab the output
        header.executeCmd('combineTool.py -M Impacts -n '+taskName+' --rMin -5 --rMax 5 -d initialFitWorkspace.root --snapshotName initialFit -m 2000 '+impactNuisanceString+' -o impacts.json')
        header.executeCmd('plotImpacts.py -i impacts.json -o impacts')

    # # Run commands
    # for c in commands:
    #     print 'Executing: '+c
    #     subprocess.call([c],shell=True)
