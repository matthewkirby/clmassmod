from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("cython libraries", ["*.pyx"],
              include_dirs = [numpy.get_include()],
              libraries = ['python2.7'])]

setup(
    name = 'cython extensions',
    ext_modules = cythonize(extensions)
)

