import ROOT,sys
sys.path.append('../')
from ROOT import TFile,TGraph,gPad
from collections import OrderedDict
import header

version = sys.argv[1]
out = TFile.Open('significance_%s.root'%version+"_Run2",'RECREATE')
tableDict = OrderedDict()
masses = range(1400,4200,200) 
handNames = {'LH':'Left-handed','RH':'Right-handed','VL':'Vector-like'}
for h in ['LH','RH','VL']:
    graph = TGraph(len(masses))

    for i,m in enumerate(masses):
        if (m/100)%2 != 0:
            val = 0
        else:
            f = TFile.Open('Limits%s/Full%s%s%s/higgsCombineTest.Significance.mH120.root'%(version,h,m,version+"_Run2"))
            t = f.Get('limit')
            t.GetEntry(0)
            val = t.limit

        graph.SetPoint(i,m,val)
        row_name = 'm = %s GeV'%(m)
        if row_name not in tableDict.keys():
            tableDict[row_name] = {}
        tableDict[row_name][handNames[h]] = '%.2f'%val

    graph.GetXaxis().SetTitle('Signal mass (GeV)')
    graph.GetYaxis().SetTitle('Significance')
    graph.SetLineColor(ROOT.kBlue)
    graph.SetLineWidth(2)
    graph.SetTitle('')
    graph.SetMinimum(0)
    graph.SetMaximum(5.1)
    graph.Draw('C A *')

    out.cd()
    graph.Write()
    gPad.Print('significance_%s_%s.pdf'%(h,version+"_Run2"))

    header.dictToLatexTable(tableDict,'significance_%s.tex'%version+"_Run2")
