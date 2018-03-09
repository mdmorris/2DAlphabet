#####################################################################################################
# 2DAlphabet.py - written by Lucas Corcodilos, 3/7/18												#
# --------------------------------------------------- 												#
# This is the wrapper and callable script that runs the full 2DAlphabet workflow. The only inputs 	#
# are a properly formatted JSON file (see input_example.txt for an example) and the options		#
# which specify which parts of the workflow to run (default is all). 								#
#####################################################################################################

#########################################################
# 						Imports							#
#########################################################
from optparse import OptionParser
import 2DAlphabet/make_card.py
import 2DAlphabet/input_organizer.py
import 2DAlphabet/build_workspace.py
import 2DAlphabet/plot_fit_results.py
import 2DAlphabet/header.py


#########################################################
#						Options 						#
#########################################################
parser = OptionParser()

parser.add_option('-i', '--input', metavar='F', type='int', action='store',
                  default   =   '',
                  dest      =   'input',
                  help      =   'JSON file to be imported')

(options, args) = parser.parse_args()

