from collections import OrderedDict
from TwoDAlphabet.helpers import roofit_form_to_TF1
from numpy.lib.function_base import piecewise
from ROOT import RooRealVar, RooFormulaVar, RooArgList, RooParametricHist2D, RooConstVar, TFormula
from binning import copy_hist_with_new_bins
import numpy as np

_subspace = ['LOW','SIG','HIGH']
class Generic2D:
    '''Wraps various input distributions in a common type so that
    distributions can easily be manipulated and compared.
    '''
    def __init__(self,name,binning):
        self.name = name
        self.binning = binning
        self.nuisances = []
        self.binVars = OrderedDict()
        self.rph = {}
        self._varStorage = []

    def _manipulate(self,name,other,operator=''):
        out = {c:Generic2D(name+'_'+c,self.binning) for c in _subspace}
        for cat in _subspace:
            new_cat_name = name+'_'+cat
            bin_pars = RooArgList()
            for ybin in range(1,len(self.binning.ybinList)):
                for xbin in range(1,len(self.binning.xbinList)):
                    new_bin_name   = '%s_bin_%s-%s'%(new_cat_name,xbin,ybin)
                    self_bin_name  = new_bin_name.replace(new_cat_name, self.name+'_'+cat)
                    other_bin_name = new_bin_name.replace(new_cat_name, other.name+'_'+cat)
                    out.binVars[new_bin_name] = RooFormulaVar(
                                                    new_bin_name, new_bin_name, '@0%s@1'%operator,
                                                    RooArgList(
                                                        self.binVars[self_bin_name],
                                                        other.binVars[other_bin_name]))
                    bin_pars.add(out.binVars[new_bin_name])

            out[cat].rph = RooParametricHist2D(
                            new_cat_name, new_cat_name,
                            self.binning.xVars[cat],
                            self.binning.yVar,
                            bin_pars,
                            self.binning.CreateHist(new_cat_name+'_template'))

        for nuisance in set(self.nuisances+other.nuisances):
            out.nuisances.append(nuisance)

        return out

    def Add(self,name,other,factor='1'):
        if factor.startswith('-'):
            op = '%s*'%factor
        elif factor == '1':
            op = '+'
        else:
            op = '+%s*'%factor
        return self._manipulate(name,other,op)

    def Multiply(self,name,other):
        return self._manipulate(name,other,'*')
    def Divide(self,name,other):
        return self._manipulate(name,other,'/')

    def Clone(self):
        out = Generic2D(self.name,self.binning)
        out.nuisances = self.nuisances
        out.binVars = self.binVars
        out.rph = self.rph
        return out

    def RooParametricHist(self):
        out = {}
        for cat in _subspace:
            cat_name = self.name+'_'+cat
            cat_hist = self.binning.CreatHist(cat_name+'_temp')
            bin_pars = RooArgList()
            for binVar in self.binVars.values():
                bin_pars.add(binVar)
            out[cat] = RooParametricHist2D(
                        cat_name, cat_name,
                        self.binning.xVars[cat],
                        self.binning.yVar,
                        bin_pars, cat_hist
            )
        return out
            

class ParametricFunction(Generic2D):
    def __init__(self,name,binning,formula,constraints={}):
        '''[summary]

        Args:
            formula ([type]): Must reference by ordinal with @. Use "x" and "y" to represent
                the "x" and "y" axes of the space. All other terms are indexed starting at 0. Ex. "@0 + x*@1 +y*@2".
            constraints (dict, optional): Map of formula parameters to constraint information. Defaults to {} in which
                case the constraint will be flat, the starting value of the parameter will be 0 with a step size of 0.1,
                and the range of the parameter will be [-1000,1000].
        '''
        super().__init__(name,binning)
        self.formula = formula
        self.nuisances = self.createFuncVars(constraints)

        for cat in _subspace:
            cat_name = name+'_'+cat
            for ybin in range(1,len(self.binning.ybinList)):
                for xbin in range(1,len(self.binning.xbinList)):
                    bin_name = '%s_bin_%s-%s'%(cat_name,xbin,ybin)
                    xConst,yConst = self.mappedBinCenter(xbin,ybin)
                    final_formula = "max(1e-9,%s)"%(self.replaceXY(xConst,yConst))
                    self.binVars[bin_name] = RooFormulaVar(
                        bin_name, bin_name,
                        final_formula,
                        RooArgList([self.nuisances[n]['obj'] for n in self.nuisances])
                    )


    def replaceXY(self,x,y):
        f = self.formula.replace(' ','')
        f = f.replace('+x','+%s'%x).replace('+y','+%s'%y)
        f = f.replace('*x','+%s'%x).replace('*y','+%s'%y)
        f = f.replace('-x','+%s'%x).replace('-y','+%s'%y)
        f = f.replace('/x','+%s'%x).replace('/y','+%s'%y)
        f = f.replace('(x','+%s'%x).replace('(y','+%s'%y)
        return f

    def getNparams(self,f):
        return TFormula('t',roofit_form_to_TF1(self.formula)).GetNpars()

    def createFuncVars(self,constraints):
        out = OrderedDict()
        for i in range(self.getNparams()):
            name = '%s_par%s'%(self.name,i)
            constraint = 'flatParam'; MIN = -1000; MAX = 1000; NOM = 0, ERROR = 0.1
            if i in constraints:
                constraint = constraints['constraint']
                MIN = constraints['MIN']
                MAX = constraints['MAX']
                NOM = constraints['NOM']
                ERROR = constraints['ERROR']

            out[name] = {'name':name, 'obj': RooRealVar(name,name,NOM,MIN,MAX), 'constraint': constraint}
            out[name]['obj'].setError(ERROR)
        return out
    
    def mappedBinCenter(self,xbin,ybin):
        '''Convert global x and y bin to RooConstVars in center of bins

        Args:
            xbin (int): X axis bin number
            ybin (int): Y axis bin number

        Returns:
            tuple of RooConstVar: First item is X const var and second is Y
        '''

        # Get this bin center 
        x_center = self.binning.GetBinCenterX(xbin)
        y_center = self.binning.GetBinCenterX(ybin)

        # Get min and range from dummy_TH2
        x_min = self.binning.xbinList[0]
        y_min = self.binning.ybinList[0]
        x_range = self.binning.xbinList[-1] - x_min
        y_range = self.binning.ybinList[-1] - y_min

        # Remap to [-1,1]
        x_center_mapped = (x_center - x_min)/x_range
        y_center_mapped = (y_center - y_min)/y_range

        return x_center_mapped,y_center_mapped

    def setFuncParam(self,varname,value):
        self.nuisances[varname].setVal(value)

    def getBinVal(self,xbin,ybin):
        return self.getBinVar(xbin,ybin).getValV()

    def getBinVar(self,xbin,ybin,c=''):
        if c == '': # using a global xbin that needs to be translated
            c, xbin = self.binning.xcatFromGlobal(xbin)
        formula_name = '%s_bin_%s-%s'%(self.name+'_'+c,xbin,ybin)
        return self.binVars[formula_name]
        
class BinnedDistribution(Generic2D):
    def __init__(self,name,inhist,binning,constant=False):
        '''[summary]

        Args:
            inhist (TH2): Input 2D histogram to build RooParametricHist2D.
            binning (Binning): Binning object used to create LOW, SIG, HIGH regions along X axis.
            constant (bool, optional): If true, use RooConstVars for bins. Defaults to False and RooRealVars are used.
        '''
        super().__init__(name,binning)
        for cat in _subspace:
            cat_name = name+'_'+cat
            cat_hist = copy_hist_with_new_bins(cat_name,'X',inhist,self.binning.xbinByCat[cat])
            for ybin in range(1,cat_hist.GetNbinsY()+1):
                for xbin in range(1,cat_hist.GetNbinsX()+1):
                    bin_name = '%s_bin_%s-%s'%(cat_name,xbin,ybin)
                    if constant:
                        self.binVars[bin_name] = RooConstVar(bin_name, bin_name, cat_hist.GetBinContent(xbin,ybin))
                    else:
                        self.binVars[bin_name] = RooRealVar(bin_name, bin_name, cat_hist.GetBinContent(xbin,ybin))
                        self.nuisances.append({'name':bin_name, 'constraint':'flatParam', 'obj': self.binVars[bin_name]})
                    self._varStorage.append(self.binVars[bin_name]) # For safety if we add shape templates            
                     
    def AddShapeTemplates(self,nuis_name,up_shape,down_shape,constraint="param"):
        nuisance_par = RooRealVar(nuis_name,nuis_name,0,-5,5)
        self.nuisances.append({'name':nuis_name, 'constraint':constraint, 'obj': nuisance_par})

        for cat in _subspace:
            rph_name = self.rph[cat].GetName()
            cat_hist_up =   copy_hist_with_new_bins(up_shape.GetName()+'_'+cat,  'X', up_shape,   self.binning.xbinByCat[cat])
            cat_hist_down = copy_hist_with_new_bins(down_shape.GetName()+'_'+cat,'X', down_shape, self.binning.xbinByCat[cat])
            for ybin in range(1,cat_hist_up.GetNbinsY()+1):
                for xbin in range(1,cat_hist_up.GetNbinsX()+1):
                    bin_name = '%s_%s_x%sy%s'%(rph_name,nuis_name,xbin,ybin)
                    self.binVar[bin_name] = singleBinInterpCombine( # change to singleBinInterpQuad to change interpolation method
                                                bin_name, self.binVars[bin_name], nuisance_par,
                                                cat_hist_up.GetBinContent(xbin,ybin),
                                                cat_hist_down.GetBinContent(xbin,ybin)
                    )

class SmoothedDistribtuion(BinnedDistribution):
    def __init__(self,name,binning,inhist):
        super().__init__(name,binning)


def InitQCDHist(data,bkgList):
    qcd = data.Clone(data.GetName().replace('data_obs','qcd'))
    for bkg in bkgList:
        qcd.Add(bkg,-1)
    return qcd


# Probably junk ------------------
def singleBinInterpQuad(name, nuis, binVar, upVal, downVal):
    nomVal = binVar.getValV()
    a,b,c = solve_quad([(-1,downVal),(0,nomVal),(1,upVal)])
    m_0, b_0 = linear_extrap((-1,downVal),[a,b,c])
    m_1, b_1 = linear_extrap((1,upVal),[a,b,c])
    pieces = [
        '(@0 < -1)*(%s*@0+%s)'%(m_0,b_0),
        '(@0 > -1 && @0 < 1)*(%s*@0**2+%s*@0+%s)'%(a,b,c),
        '(@0 > 1)*(%s*@0+%s)'%(m_1,b_1)]
    
    return RooFormulaVar(name, name, '@1*(({0})-{1})/{1}'.format('+'.join(pieces),nomVal), RooArgList(nuis,binVar))

def singleBinInterpCombine(name, nuis, binVar, upVal, downVal):
    nomVal = binVar.getValV()
    pieces = '(abs(@0) < 1)*(0.125*@0*(pow(@0,2)*(3*pow(@0,2) - 10) + 15))-(@0 < -1)+(@0 > 1)'
    alpha = '@1*( (@0/2) * (({up}-{down})+({up}+{down}-2*{nom})*({piecewise})) )/{nom}'.format(up=upVal,down=downVal,nom=nomVal,piecewise=pieces)
    return RooFormulaVar(name,name,alpha,RooArgList(nuis,binVar))

class quadratic:
    def __init__(self, params):
        self.a = params[0]
        self.b = params[1]
        self.c = params[2]
    def y(self,x):
        return self.a * x**2 + self.b * x + self.c
    def slope(self,x):
        return 2*self.a*x + self.b

class line:
    def __init__(self, slope, intercept):
        self.m = slope
        self.b = intercept
    def y(self,x):
        return self.m*x + self.b

def solve_quad(threepts):
    # A = [[x_1^2, x_1, 1],
    #      [x_2^2, x_2, 1],
    #      [x_3^2, x_3, 1]]
    # X = [[a],[b],[c]]
    # B = [[y_1], [y_2], [y_3]]
    A, B = [], []
    for pt in threepts:
        x, y = pt[0], pt[1]
        A.append([x**2, x, 1])
        B.append([y])
    A = np.array(A)
    B = np.array(B)
    out = [x[0] for x in np.linalg.solve(A,B)]
    return out[0], out[1], out[2]

def linear_extrap(point,quad_params):
    q = quadratic(quad_params)
    slope = q.slope(point[0])
    const = -1*slope*point[0]+point[1]
    return slope,const

print (solve_quad([(0,0),(-1,1),(1,1)]))
