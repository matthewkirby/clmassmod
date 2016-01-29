'''Common utilities and things that are needed everywhere'''

############

import imp

############



def readConfiguration(configname):

    configmodule = imp.load_source('config', configname)

    config = configmodule.__dict__

    configure(config)

    return config

#######################

def configure(config):

    for val in config.itervalues():
        if hasattr(val, 'configure'):
            val.configure(config)
    

########################

def buildObject(modulename, classname, *args, **kwds):

    if modulename.lower() == 'none' or classname.lower() == 'none':
        return None

    aModule = importlib.import_module(modulename)
    aClass = getattr(aModule, classname)
    anObject = aClass(*args, **kwds)

    return anObject
    
#####

class Composite(object):

    def __init__(self, *pickers):
        self.pickers = pickers

    def configure(self, config):

        for picker in self.pickers:
            if hasattr(picker, 'configure'):
                picker.configure(config)

    def __call__(self, sim):

        curcat = sim
        for picker in self.pickers:
            curcat = picker(curcat)

        return curcat
