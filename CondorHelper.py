import subprocess
from optparse import OptionParser
import time


parser = OptionParser()

parser.add_option('-r', '--runscript', metavar='FILE', type='string', action='store',
                default   =   '',
                dest      =   'runscript',
                help      =   'Run template to use')
parser.add_option('-i', '--inputs', metavar='FILE', type='string', action='store',
                default   =   '',
                dest      =   'inputs',
                help      =   'Inputs to send along')
parser.add_option('-a', '--args', metavar='FILE', type='string', action='store',
                default   =   '',
                dest      =   'args',
                help      =   'Text file with python arguments')
parser.add_option('-t', '--tarball', metavar='FILE', type='string', action='store',
                default   =   'tarball.tgz',
                dest      =   'tarball',
                help      =   'Name of tarball to make and send.')

(options, args) = parser.parse_args()

commands = []

print options.args

#if 'bstar' in options.runscript: commands.append('source condor/prep.csh')

# Tar stuff
if options.inputs != '':
    commands.append('rm '+options.tarball)
    commands.append("tar czvf "+options.tarball+" "+options.inputs+' --exclude="*.tgz" --exclude="*.png" --exclude="higgsCombinegof_*.root" --exclude="*.sh" --exclude="*.std*"')

# Make JDL from template
timestr = time.strftime("%Y%m%d-%H%M%S")
out_jdl = 'temp_'+timestr+'_jdl'
commands.append("sed 's$TEMPSCRIPT$"+options.runscript+".tmp$g' condor/templates/jdl_template > "+out_jdl)
commands.append("sed -i 's$TEMPARGS$"+options.args+"$g' "+out_jdl)
commands.append("sed -i 's$TEMPTAR$"+options.tarball+"$g' "+out_jdl)
commands.append("sed 's$TEMPTAR$"+options.tarball+"$g' "+options.runscript+" > "+options.runscript+'.tmp')
commands.append("condor_submit "+out_jdl)
commands.append("mv "+out_jdl+" condor/jdls/")
# commands.append("condor_q lcorcodi")

for s in commands:
    print s
    subprocess.call([s],shell=True)
