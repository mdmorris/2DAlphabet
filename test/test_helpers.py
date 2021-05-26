import sys
for p in sys.path: print (p)
from TwoDAlphabet.src.helpers import *

def test__open_json():
    test = {
        'this': 'is',
        'a' : {
            'test': 1,
            'dictionary': False
        }
    }
    assert (open_json("test/testjson.json",False) == test)
    open_json('test/twoDtest.json')


def test__arg_dict_to_list():
    pass
def test__parse_arg_dict():
    pass
