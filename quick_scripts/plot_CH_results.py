import ROOT
from ROOT import *

import header
from header import makeCan

gStyle.SetOptStat(0)

gen = 'dataSidebandFewerbins_tt'

fd = TFile.Open('../'+gen+'/fitDiagnostics.root')
post_file = TFile.Open('../'+gen+'/postfitshapes.root')

data_hists_y = []
bkg_hists_y = []
data_hists_x = []
bkg_hists_x = []
total_hists_y = []
total_hists_x = []

for thisCat in ['fail','pass']:
    data = post_file.Get(thisCat + '_postfit/data_obs')
    qcd = post_file.Get(thisCat + '_postfit/qcd')
    ttbar = post_file.Get(thisCat + '_postfit/ttbar')
    total = post_file.Get(thisCat+ '_postfit/TotalBkg')

    projs_y = {
        # "full":{
        #     "xlow":0,
        #     "xhigh":-1
        # },
        "signal":{
            "xlow":6,
            "xhigh":7
        },
        "SB_low":{
            "xlow":0,
            "xhigh":5
        },
        "SB_high":{
            "xlow":8,
            "xhigh":-1
        }
    }


    # Full projection
    for key in ['SB_low','signal','SB_high']:
        data_proj_y = data.ProjectionY('data_'+key+'_'+thisCat+'_py',projs_y[key]['xlow'],projs_y[key]['xhigh'],'e')
        qcd_proj_y = qcd.ProjectionY('qcd_'+key+'_'+thisCat+'_py',projs_y[key]['xlow'],projs_y[key]['xhigh'],'e')
        ttbar_proj_y = ttbar.ProjectionY('ttbar_'+key+'_'+thisCat+'_py',projs_y[key]['xlow'],projs_y[key]['xhigh'],'e')
        total_proj_y = total.ProjectionY('total_'+key+'_'+thisCat+'_py',projs_y[key]['xlow'],projs_y[key]['xhigh'],'e')

        data_proj_y.SetTitle(key+' - '+thisCat)
        qcd_proj_y.SetTitle(key+' - '+thisCat)
        ttbar_proj_y.SetTitle(key+' - '+thisCat)

        qcd_proj_y.SetMinimum(0)
        data_proj_y.SetMinimum(0)
        ttbar_proj_y.SetMinimum(0)

        data_hists_y.append(data_proj_y)
        bkg_hists_y.append([ttbar_proj_y,qcd_proj_y])
        total_hists_y.append(total_proj_y)

    projs_x = {
        # "full":{
        #     "xlow":0,
        #     "xhigh":-1
        # },
        "turn-on":{
            "ylow":0,
            "yhigh":4
        },
        "falling":{
            "ylow":5,
            "yhigh":10
        },
        "tail":{
            "ylow":11,
            "yhigh":-1
        }
    }


    # Full projection
    for key in ['turn-on','falling','tail']:
        data_proj_x = data.ProjectionX('data_'+key+'_'+thisCat+'_px',projs_x[key]['ylow'],projs_x[key]['yhigh'],'e')
        qcd_proj_x = qcd.ProjectionX('qcd_'+key+'_'+thisCat+'_px',projs_x[key]['ylow'],projs_x[key]['yhigh'],'e')
        ttbar_proj_x = ttbar.ProjectionX('ttbar_'+key+'_'+thisCat+'_px',projs_x[key]['ylow'],projs_x[key]['yhigh'],'e')
        total_proj_x = total.ProjectionX('total_'+key+'_'+thisCat+'_px',projs_x[key]['ylow'],projs_x[key]['yhigh'],'e')

        data_proj_x.SetTitle(key+' - '+thisCat)
        qcd_proj_x.SetTitle(key+' - '+thisCat)
        ttbar_proj_x.SetTitle(key+' - '+thisCat)

        qcd_proj_x.SetMinimum(0)
        data_proj_x.SetMinimum(0)
        ttbar_proj_x.SetMinimum(0)

        data_hists_x.append(data_proj_x)
        bkg_hists_x.append([qcd_proj_x,ttbar_proj_x])
        total_hists_x.append(total_proj_x)



# for i in data_hists:
#     print i.GetName(),
# print '\n'
# for j in bkg_hists:
#     for k in j:
#         print k.GetName()
#     print '\n'

makeCan('Y_Projections','',data_hists_y,bkg_hists_y,[kRed,kYellow],xtitle='M_{tt}')
makeCan('X_Projections','',data_hists_x,bkg_hists_x,[kYellow,kRed],xtitle='M_{t}')

makeCan('totalBkg_x_projections','',total_hists_x,xtitle='M_{tt}')
makeCan('totalBkg_y_projections','',total_hists_y,xtitle='M_{t}')

    # Signal projection

    # Sideband projection

    # Sideband low projection

    # Sideband high projection
