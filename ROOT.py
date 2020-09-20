# Makes sure that --help works for python modules

import imp, sys 
for path in sys.path: 
    try: 
        modtuple = imp.find_module( 'ROOT', [path] ) 
    except ImportError: 
        pass 
    else: 
        if modtuple and modtuple[1] not in __file__: 
            sys.modules[ 'ROOT' ] = imp.load_module( 'ROOT', *modtuple ) 
            import ROOT 
            ROOT.PyConfig.IgnoreCommandLineOptions = True
