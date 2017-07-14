import struct
import numpy as np

#See http://docs.python.org/2/library/struct.html  for a description of Format strings

######################

def readArray(input, format, shape, order='C', endian='@'):
    nelem = reduce(lambda x,y: x*y, shape)
    structformat = '{0}{1}{2}'.format(endian, nelem, format)
    rawdat = np.array(struct.unpack(structformat, input.read(struct.calcsize(structformat))))
    return rawdat.reshape(shape, order=order)


#########################


def readVal(input, format, endian='@'):
    actualFormat = '{0}{1}'.format(endian, format)
    return struct.unpack(format, input.read(struct.calcsize(format)))[0]


