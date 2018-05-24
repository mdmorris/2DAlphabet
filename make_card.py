#####################################################################################################
# make_card.py - written by Lucas Corcodilos, 3/19/18                                               #
# ---------------------------------------------------------                                         #
# This script builds the data card used by Combine from the JSON input by the user. You'll see      #
# several spots where it looks like I have extra spaces in my strings. This is deliberate. Every    #
# line is collimated by separating each 'column' with a space. The function colliMate then splits   #
# the line based on the position of spaces. Thus if you want to skip a column, you put two spaces.  #
#####################################################################################################

import header
from header import colliMate

def main(inputConfig, blinded, tag):
    # Recreate file
    card_new = open('card_'+tag+'.txt','w')
    
    #######################################################
    # imax (bins), jmax (backgrounds), kmax (systematics) #
    #######################################################
    imax = '2'                      # pass, fail
    channels = ['pass', 'fail']

    # Get the length of the list of all process that have CODE 2 or 3 (and ignore "HELP" key) and add 1 for qcd (which won't be in the inputConfig)
    jmax = str(len([proc for proc in inputConfig['PROCESS'].keys() if proc != 'HELP' and inputConfig['PROCESS'][proc]['CODE'] >= 2]) + 1)
    # Get the length of the lsit of all systematics (and ignore "HELP" key)
    kmax = str(len([syst for syst in inputConfig['SYSTEMATIC'].keys() if syst != 'HELP']))

    card_new.write('imax '+imax+'\n')      
    card_new.write('jmax '+jmax+'\n')
    card_new.write('kmax '+kmax+'\n')
    card_new.write('-'*120+'\n')

    ##########
    # Shapes #
    ##########
    procs_with_systs = [proc for proc in inputConfig['PROCESS'].keys() if proc != 'HELP' and len(inputConfig['PROCESS'][proc]['SYSTEMATICS']) != 0]
    procs_without_systs = [proc for proc in inputConfig['PROCESS'].keys() if proc != 'HELP' and len(inputConfig['PROCESS'][proc]['SYSTEMATICS']) == 0]
    procs_without_systs.append('qcd')   # Again, qcd not in the input JSON but needs to be in the combine card!

    for proc in procs_without_systs:
        card_new.write(colliMate('shapes  '+proc+' * base_'+tag+'.root w_2D:'+proc+'_$CHANNEL\n'))
    for proc in procs_with_systs:
        card_new.write(colliMate('shapes  '+proc+' * base_'+tag+'.root w_2D:'+proc+'_$CHANNEL w_2D:'+proc+'_$CHANNEL_$SYSTEMATIC\n'))

    card_new.write('-'*120+'\n')

    ####################################
    # Set bin observation values to -1 #
    ####################################
    tempString = 'bin  '
    for chan in channels:
        tempString += (chan+' ')
    tempString += '\n'
    card_new.write(colliMate(tempString))

    tempString = 'observation  '
    for ichan in range(int(imax)):
        tempString += '-1 '
    tempString += '\n'
    card_new.write(colliMate(tempString))

    card_new.write('-'*120+'\n')

    ######################################################
    # Tie processes to bins and rates and simultaneously #
    # create the systematic uncertainty rows             #
    ######################################################
    bin_line = 'bin  '
    processName_line = 'process  '
    processCode_line = 'process  '
    rate_line = 'rate  '
    syst_lines = {}

    # Fill syst_lines with keys to initialized strings
    for syst in [systematic for systematic in inputConfig['SYSTEMATIC'].keys() if systematic != 'HELP']:
        if inputConfig['SYSTEMATIC'][syst]['CODE'] == 0 or inputConfig['SYSTEMATIC'][syst]['CODE'] == 1:        # lnN
            syst_type = 'lnN'
        elif inputConfig['SYSTEMATIC'][syst]['CODE'] == 2 or inputConfig['SYSTEMATIC'][syst]['CODE'] == 3:      # shape
            syst_type = 'shape'
        else:
            print 'Systematic ' + syst + ' does not have one of the four allowed codes (0,1,2,3). Quitting.'
            quit()
        
        syst_lines[syst] = syst + ' ' + syst_type + ' '

    signal_procs = [proc for proc in inputConfig['PROCESS'].keys() if proc != 'HELP' and proc != 'data_obs' and inputConfig['PROCESS'][proc]['CODE'] == 0]
    MC_bkg_procs = [proc for proc in inputConfig['PROCESS'].keys() if proc != 'HELP' and proc != 'data_obs' and (inputConfig['PROCESS'][proc]['CODE'] == 2 or inputConfig['PROCESS'][proc]['CODE'] == 3)]
    all_procs = [proc for proc in inputConfig['PROCESS'].keys() if proc != 'HELP' and proc != 'data_obs']
    all_procs.append('qcd')

    for chan in channels:
        for proc in all_procs:
            # Start lines
            bin_line += (chan+' ')
            processName_line += (proc+' ')

            # If signal
            if proc in signal_procs:
                processCode_line += (str(0-signal_procs.index(proc))+' ')
                rate_line += ('-1 ')

            # If bkg
            if proc in MC_bkg_procs:
                processCode_line += (str(MC_bkg_procs.index(proc)+1)+' ')
                if inputConfig['PROCESS'][proc]['CODE'] == 2:       # No floating normalization
                    rate_line += '-1 '                                            
                elif inputConfig['PROCESS'][proc]['CODE'] == 3:     # Floating normalization        
                    rate_line += '1 '

            # If qcd
            if proc == 'qcd':
                processCode_line += (str(len(MC_bkg_procs)+2)+' ')
                rate_line += '1 '

            # Fill systematic lines
            for syst_line_key in syst_lines.keys():
                # If the systematic is applicable to the process
                if proc != 'qcd':
                    if syst_line_key in inputConfig['PROCESS'][proc]['SYSTEMATICS']:
                        # If symmetric lnN...
                        if inputConfig['SYSTEMATIC'][syst_line_key]['CODE'] == 0:
                            thisVal = str(inputConfig['SYSTEMATIC'][syst_line_key]['VAL'])
                        # If asymmetric lnN...
                        elif inputConfig['SYSTEMATIC'][syst_line_key]['CODE'] == 1:
                            thisVal = str(inputConfig['SYSTEMATIC'][syst_line_key]['VALUP']) + '/' + str(inputConfig['SYSTEMATIC'][syst_line_key]['VALDOWN'])
                        # If shape...
                        else:
                            thisVal = str(inputConfig['SYSTEMATIC'][syst_line_key]['SCALE'])
                    # Otherwise place a '-'
                    else:
                        thisVal = '-'  
                else:
                    thisVal = '-' 

                syst_lines[syst_line_key] += (thisVal+' ')



    card_new.write(colliMate(bin_line+'\n'))
    card_new.write(colliMate(processName_line+'\n'))
    card_new.write(colliMate(processCode_line+'\n'))
    card_new.write(colliMate(rate_line+'\n'))
    card_new.write('-'*120+'\n')

    ############################
    # Systematic uncertainties #
    ############################
    for line_key in syst_lines.keys():
        card_new.write(colliMate(syst_lines[line_key]+'\n'))


    #################################################
    # Mark floating values as flatParams            # 
    # In general we float                           #
    # - fitParam*_*                                 #
    #                                               #
    # If there's a renormalized process we float    #
    # - process_norm                                #
    # - Fail_bin_x-y_process_nominal                #
    # - Fail_bin_x-y_init                           #
    #                                               #
    # Otherwise we float                            #
    # - Fail_bin_x-y                                #
    #################################################
    if 'XFORM' in inputConfig['FIT'].keys() and 'YFORM' in inputConfig['FIT'].keys():
        nxparams = max([int(param[1:]) for param in inputConfig['FIT'].keys() if param.find('X') != -1 and param != 'XFORM'])
        nyparams = max([int(param[1:]) for param in inputConfig['FIT'].keys() if param.find('Y') != -1 and param != 'YFORM'])

        for nparams in [nxparams, nyparams]:
            if nparams == nxparams:
                thisVar = 'X'
            else:
                thisVar = 'Y'
            for ip in range(1,nparams+1):
                card_new.write(colliMate('fitParam'+thisVar+'_'+str(ip)+' flatParam\n',22))
    elif 'FORM' in inputConfig['FIT'].keys():
        for coeff in [key for key in inputConfig['FIT'].keys() if key != 'HELP' and key.find('FORM') == -1]:
            if 'LOW' in inputConfig['FIT'][coeff].keys() and 'HIGH' in inputConfig['FIT'][coeff].keys():
                lower_coeff = coeff.lower()
                card_new.write(colliMate('polyCoeff_'+lower_coeff+' flatParam\n',22))


    # Check if we have any renormalized MCs, store which processes, and declare the _norm as a flatParam
    renormFlag = False
    renormProcesses = []
    for process in inputConfig['PROCESS']:
        if process != 'HELP' and inputConfig['PROCESS'][process]['CODE'] == 3:
            card_new.write(colliMate(process + '_norm flatParam\n',22))
            renormProcesses.append(process)
            renormFlag = True

    # Clearer code if I grab all of this
    xbins_low = inputConfig['BINNING']['X']['LOW']
    xbins_high = inputConfig['BINNING']['X']['HIGH']
    xbins_n = inputConfig['BINNING']['X']['NBINS']
    xbins_width = float((xbins_high - xbins_low)/xbins_n)

    ybins_n = inputConfig['BINNING']['Y']['NBINS']


    # Now float the failing bins
    for ybin in range(1,ybins_n+1):
        # Setup xbins with flatParams needed based on if blinded or not
        if blinded:
            # Calculate the number of bins in low and high x
            xbins_sigstart = inputConfig['BINNING']['X']['SIGSTART']
            xbins_sigend = inputConfig['BINNING']['X']['SIGEND']
            xbins_n_low = float(xbins_sigstart - xbins_low)/xbins_width
            xbins_n_high = float(xbins_high - xbins_sigend)/xbins_width
            #       Low                           High
            xbins = range(1,int(xbins_n_low+1)) + range(int(xbins_n-xbins_n_high+1),int(xbins_n+1))
        else:
            xbins = range(1,xbins_n+1)

        # Write the flatParams
        for xbin in xbins:
            # if renormFlag:
            #     card_new.write(colliMate('Fail_bin_'+str(xbin)+'-'+str(ybin)+'_init flatParam\n',22))
            #     # for process in renormProcesses:
            #     #     card_new.write(colliMate('Fail_bin_'+str(xbin)+'-'+str(ybin)+'_'+process+'_nominal flatParam\n',22))
            # else:
            card_new.write(colliMate('Fail_bin_'+str(xbin)+'-'+str(ybin)+' flatParam\n',22))
            
       
    card_new.close()

