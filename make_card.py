# Finally do some minor card prep
card_template = open('card_2Dv3p1.tmp','r')
card_new = open('card_2Dv3p1.txt','w')
card_new.write(card_template.read())
card_template.close()
card_new.write('\n')
flatParamLines = []
for sband in ['Low','High']:
    for ix in range(1,dDists['qcd_fail']['xbins_'+sband]+1):
        for iy in range(1,dDists['qcd_fail']['ybins_'+sband]+1):
            flatParamLines.append('Fail'+sband+'_bin_'+str(ix)+'-'+str(iy)+' flatParam \n')
card_new.writelines(flatParamLines)
card_new.close()

