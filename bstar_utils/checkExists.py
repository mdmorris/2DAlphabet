'''Check for any files or histograms that don't exist in config.
System argument is the configuration file to check. Assumes
that if a "path" variable exists, it's called `path`.
'''

import ROOT, sys, os

from header import openJSON

config = openJSON(sys.argv[1])

path = config['GLOBAL']['path']

for p in config['PROCESS'].keys():
    if 'help' in p.lower(): continue

    print ('Checking '+p)
    full_path = config['PROCESS'][p]['FILE']
    if 'path' in full_path:
        full_path = full_path.replace('path',path)

    if not os.path.exists(full_path):
        print ('%s: %s does not exist'%(p,full_path))

    else:
        f = ROOT.TFile.Open(full_path)
        hpassName = config['PROCESS'][p]['HISTPASS']
        hfailName = config['PROCESS'][p]['HISTFAIL']
        if not f.Get(hpassName):
            print ('%s: %s does not exist'%(p,hpassName)) 
        if not f.Get(hfailName):
            print ('%s: %s does not exist'%(p,hfailName))
