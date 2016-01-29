'''Common utilities and things that are needed everywhere'''

############

import importlib

############



def readConfiguration(configname):

    #can this take a file path?
    configmodule = importlib.import_module(configname)
    
    

    

########################

def buildObject(modulename, classname, *args, **kwds):

    if modulename.lower() == 'none' or classname.lower() == 'none':
        return None

    aModule = importlib.import_module(modulename)
    aClass = getattr(aModule, classname)
    anObject = aClass(*args, **kwds)

    return anObject
    
#####
