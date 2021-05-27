class AlphaWrap:
    '''Wraps various input distributions in a common type so that
    distributions can easily be manipulated and compared.
    '''
    def __init__(self,inObj):
        self.inObj = inObj

    def __add__(self,other):
        pass
    def __sub__(self,other):
        pass
    def __mul__(self,other):
        pass
    def __div__(self,other):
        pass

    def __radd__(self,other):
        return self+other
    def __rsub__(self,other):
        raise NotImplementedError
    def __rmul__(self,other):
        return self*other
    def __rdiv__(self,other):
        raise NotImplementedError

    def RooParametricHist(self):
        pass
    def RooDataHist(self):
        pass
