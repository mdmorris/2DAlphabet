import ROOT
from ROOT import *

import header, os, subprocess, copy

class BinnedParametricFunc():
    def __init__(self,fitDict,name,dummy_TH2=None,funcname=''):
        self.fitDict = fitDict
        self.funcname = funcname
        self.name = name
        self.dummy_TH2 = dummy_TH2
        self.funcVars = {}
        self.binVars = {} # RooFormulaVars for each bin for later evaluation
        self.allVars = []
        self.formula = self.getRooFunctionForm()

        # Initialize vars 
        for c in sorted(self.fitDict.keys()):
            input_param_vals = self.fitDict[c]
            if c.isdigit():
                this_nom = input_param_vals['NOMINAL']
                if 'MIN' in input_param_vals.keys() and 'MAX' in input_param_vals.keys():
                    this_min = input_param_vals['MIN']
                    this_max = input_param_vals['MAX']
                    varname = self.funcname+c+'_'+self.name
                    self.funcVars[varname] = RooRealVar(varname, varname, this_nom, this_min, this_max)

                    if 'ERROR_UP' in input_param_vals.keys() and 'ERROR_DOWN' in input_param_vals.keys():
                        self.funcVars[varname].setAsymError(input_param_vals['ERROR_DOWN'],input_param_vals['ERROR_UP'])
                    elif 'ERROR' in self.fitDict[c].keys():
                        self.funcVars[varname].setError(input_param_vals['ERROR'])

                else:
                    input_break = raw_input('WARNING: Upper and lower bounds on fit parameter ' + c + ' not defined. Please type "N" if you would like to quit or enter if you would like to treat the parameter as constant')
                    if input_break == 'N':
                        quit()
                    else:
                        self.funcVars[varname] = RooConstVar(varname,varname,this_nom)

    def getRooFunctionForm(self):
        params = [int(param) for param in self.fitDict.keys() if param != 'FORM' and param != 'HELP' and param != "START"]
        if len(params) > 0:
            nCoeffs = max(params)+1
        else: nCoeffs = 0
        
        gFormula = self.fitDict['FORM'].replace('+x','+@'+str(nCoeffs)).replace('+y','+@'+str(nCoeffs+1))
        gFormula = gFormula.replace('*x','*@'+str(nCoeffs)).replace('*y','*@'+str(nCoeffs+1))
        gFormula = gFormula.replace('-x','-@'+str(nCoeffs)).replace('-y','-@'+str(nCoeffs+1))
        gFormula = gFormula.replace('/x','/@'+str(nCoeffs)).replace('/y','/@'+str(nCoeffs+1))
        gFormula = gFormula.replace('(x','(@'+str(nCoeffs)).replace('(y','(@'+str(nCoeffs+1))
        gFormula = gFormula.replace(' x',' @'+str(nCoeffs)).replace(' y',' @'+str(nCoeffs+1))
        return gFormula
    
    def Eval(self,xbin,ybin,name=''):
        '''Evaluate and store the value of the parametric function in provided global bins

        Args:
            xbin (int): X axis bin number
            ybin (int): Y axis bin number

        Returns:
            RooFormulaVar: The RooFormulaVar for this bin of the parametric function
        '''
        xConst,yConst = self.binToConst(xbin,ybin)

        # Get category name for naming
        c = xConst.GetName().split('_')[2]

        if name != '':
            full_formula_name = name
        else:
            full_formula_name = 'formula_bin_'+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name

        formula_list = RooArgList()
        for c in self.getFuncVarNames():
            formula_list.add(self.funcVars[c])

        formula_list.add(xConst)
        formula_list.add(yConst)
        func_val = RooFormulaVar(full_formula_name,full_formula_name,"max(1e-9,"+self.formula+")",formula_list)

        self.storeFuncBin(func_val,full_formula_name)         

        return func_val
    
    def binToConst(self,xbin,ybin):
        '''Convert global x and y bin to RooConstVars in center of bins

        Args:
            xbin (int): X axis bin number
            ybin (int): Y axis bin number

        Returns:
            tuple of RooConstVar: First item is X const var and second is Y
        '''

        # Get this bin center 
        x_center = self.dummy_TH2.GetXaxis().GetBinCenter(xbin)
        y_center = self.dummy_TH2.GetYaxis().GetBinCenter(ybin)

        # Get min and range from dummy_TH2
        x_min = self.dummy_TH2.GetXaxis().GetBinLowEdge(1)
        y_min = self.dummy_TH2.GetYaxis().GetBinLowEdge(1)
        x_range = self.dummy_TH2.GetXaxis().GetBinUpEdge(self.dummy_TH2.GetNbinsX()) - x_min
        y_range = self.dummy_TH2.GetYaxis().GetBinUpEdge(self.dummy_TH2.GetNbinsY()) - y_min

        # Remap to [-1,1]
        x_center_mapped = (x_center - x_min)/x_range
        y_center_mapped = (y_center - y_min)/y_range

        # And assign it to a RooConstVar 
        x_const = RooConstVar("ConstVar_x_"+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+str(xbin)+'-'+str(ybin)+'_'+self.name,x_center_mapped)
        y_const = RooConstVar("ConstVar_y_"+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+str(xbin)+'-'+str(ybin)+'_'+self.name,y_center_mapped)
        
        self.allVars.append(x_const)
        self.allVars.append(y_const)

        return x_const,y_const

    def storeFuncBin(self,val,name):
        self.binVars[name] = val

    def setFuncParam(self,varname,value):
        self.funcVars[varname].setVal(value)

    def getFuncBinVal(self,c,xbin,ybin):
        formula_name = 'formula_bin_'+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name
        return self.binVars[formula_name].getValV()

    def getFuncBinRRV(self,c,xbin,ybin):
        formula_name = 'formula_bin_'+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name
        return self.binVars[formula_name]

    def getFuncVarNames(self):
        return sorted(self.funcVars.keys())

class RpfHandler(BinnedParametricFunc):
    def __init__ (self,fitDict,name,dummy_TH2=None,tag=''):
        BinnedParametricFunc.__init__(self,fitDict,name,dummy_TH2,'rpf')

    # Removes everything but basic values for storage
    def getReducedCopy(self):
        newcopy = RpfHandler(self.fitDict,self.name,self.dummy_TH2)
        newcopy.funcVars = copy.deepcopy(self.funcVars)
        return newcopy
