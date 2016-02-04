'''Common utilities and things that are needed everywhere'''

############

import imp

############



def readConfiguration(configname):

    configmodule = imp.load_source('config', configname)

    config = configmodule.__dict__

    runConfigure(config)

    return config

#######################

def runConfigure(config):

    toconfigure = []
    for val in config.itervalues():
        if hasattr(val, 'configure'):
            toconfigure.append(val)

    for val in toconfigure:
        val.configure(config)
    

########################


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
