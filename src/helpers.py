import argparse, subprocess, json, ROOT, os
from collections import defaultdict
from contextlib import contextmanager

# Function stolen from https://stackoverflow.com/questions/9590382/forcing-python-json-module-to-work-with-ascii
def open_json(f,twoDconfig=True):
    '''Open a JSON file. Specify twoDconfig to true if this is a 2DAlphabet 
    configuration file.

    Function adapted from https://stackoverflow.com/questions/9590382/forcing-python-json-module-to-work-with-ascii

    @param f (str): File name and path.
    @param twoDconfig (bool, optional): Set to True if the JSON file is a 2DAlphabet
        configuration file. Defaults to True.

    Returns:
        dict: JSON opened as a python dictionary.
    '''
    with open(f) as fInput_config:
        input_config = json.load(fInput_config, object_hook=ascii_encode_dict)  # Converts most of the unicode to ascii

        if twoDconfig:
            for process in [proc for proc in input_config['PROCESS'].keys() if proc != 'HELP']:
                for index,item in enumerate(input_config['PROCESS'][process]['SYSTEMATICS']):           # There's one list that also
                    input_config['PROCESS'][process]['SYSTEMATICS'][index] = item.encode('ascii')  

    return input_config

def ascii_encode_dict(data): 
    '''Convert a unicode encoded dictionary into ascii.

    @param data (dict): Input dictionary.

    Returns:
        dict: Dict encoded with ascii instead of unicode.
    '''
    ascii_encode = lambda x: x.encode('ascii') if isinstance(x, unicode) else x 
    return dict(map(ascii_encode, pair) for pair in data.items())

def arg_dict_to_list(indict):
    '''Convert dictionary of arguments into a list of
    `<key>=<value>` strings.

    @param indict (dict): Dictionary of arguments where keys are
        the argument name and values are the value for the argument to take.

    Returns:
        list(str): List of `<key>=<value>` strings
    '''
    return ['%s=%s'%(k,v) for k,v in indict.items()]

def parse_arg_dict(parser,indict):
    '''Use ArgumentParser to parse arguments in an input dictionary
    where the key is the arg name and the value is the value for the arg to take.

    @param parser (ArgumentParser): Object to update with indict key-values.
    @param indict (dict): Dictionary holding argument names (keys) and values (values).

    Returns:
        Namespace: Object storing each argument as an attribute of the Namespace.
    '''
    default_args_dict = vars(parser.parse_args([]))
    default_args_dict.update(indict)
    new_namespace = argparse.Namespace()
    new_namespace.__dict__.update(default_args_dict)
    return parser.parse_args(namespace=new_namespace)

def execute_cmd(cmd,dryrun=False):
    '''Print and execute a command-line command via subprocess.call().
    If dryrun==True, only print.

    @param cmd (str): Command to execute as a subprocess.
    @param dryrun (bool, optional): If True, only print. Defaults to False.
    '''
    print ('Executing: '+cmd)
    if not dryrun:
        subprocess.call([cmd],shell=True)

@contextmanager
def cd(newdir):
    '''Change active directory so that the local directory becomes newdir.
    This affects everything from subprocess calls to ROOT Print() and SaveAs()
    commands. Usage:

    with cd('/path/to/stuff/'):
        ...
        <code acting in /path/to/stuff/ directory>
        ...

    @param newdir (str): Directory to cd to within the `with` statement.
    '''
    print ('cd '+newdir)
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def make_RDH(myTH2,RAL_vars,altname=''):
    '''Create a RooDataHist from the input TH2 and RooArgList of
    axis variables.

    @param myTH2 (TH2): Histogram to turn into RooDataHist.
    @param RAL_vars (RooArgList): List of RooRealVars representing the axes.
    @param altname (str, optional): Alternate name. Defaults to '' in which case the TH2 name is used.

    Returns:
        RooDataHist
    '''
    name = myTH2.GetName() if altname == '' else altname
    thisRDH = ROOT.RooDataHist(name,name,RAL_vars,myTH2)
    return thisRDH

# def make_RHP(myRDH,RAL_vars):
#     name = myRDH.GetName()
#     thisRAS = ROOT.RooArgSet(RAL_vars)
#     thisRHP = ROOT.RooHistPdf(name,name,thisRAS,myRDH)
#     return thisRHP

def dict_copy(inDict,structureOnly=False):
    '''Recursively copy a dictionary with the option to only copy the
    dictionary structure.

    @param inDict (dict): Input dictionary.
    @param structureOnly (bool, optional): If True, only copy the structure of the dictionary
        with stand-in value of 0. Defaults to False.

    Returns:
        dict: Copy.
    '''
    newDict = {}
    for k1,v1 in inDict.items():
        if isinstance(v1,dict):
            newDict[k1] = dict_copy(v1)
        else:
            if structureOnly: newDict[k1] = 0
            else: newDict[k1] = v1
    return newDict

def nested_dict(level,t):
    '''Create an empty dictionary nested to <level> levels
    and with type <t> as the default stand-in in the inner-most dictionary.

    @param level (int): Number of nested levels (ie. max number of keys that can be used).
    @param t (type): Type to provide to defaultdict to determine the stand-in values.

    Returns:
        [type]: [description]
    '''
    if level > 1:
        out = defaultdict(lambda: nested_dict(level-1,t))
    else:
        out = defaultdict(t)
    return out

def roofit_form_to_TF1(RFVform,shift=0): # shift tells function how much to shift the indices of the coefficients by
    '''Convert a RooFit formula (using @) to a ROOT.TF1 type formula (using []).

    @param RFVform (str): RooFit formula string.
    @param shift (int, optional): Number of indices to shift the coefficients in the case of abormal
        indexing. Defaults to 0.

    Returns:
        str: ROOT.TF1 formula string.
    '''
    TF1form = ''
    lookingForDigit = False
    for index,char in enumerate(RFVform):
        if char == '@':
            TF1form+='['
            lookingForDigit = True
        elif lookingForDigit:
            if char.isdigit():
                if RFVform[index+1].isdigit() == False:     # if this is the last digit
                    TF1form+=str(int(char)+shift)
            else:
                TF1form+=']'+char
                lookingForDigit = False
        else:
            TF1form+=char

    return TF1form

def set_hist_maximums(histList,factor=1.1):
    '''Take in a list of histograms and set the maximum of each
    to the maximum of the group multiplied by `factor`.

    @param histList (list(TH1)): List of histograms to compare and set.
    @param factor (float, optional): Defaults to 1.1.

    Returns:
        list(TH1): Histograms with maximum set.
    '''
    # Set the max of the range so we can see all three histograms on the same plot
    out = []
    yMax = histList[0].GetMaximum()
    for h in range(1,len(histList)):
        if histList[h].GetMaximum() > yMax:
            yMax = histList[h].GetMaximum()
    for h in histList:
        h.SetMaximum(yMax*factor)
        out.append(h)
    return out
