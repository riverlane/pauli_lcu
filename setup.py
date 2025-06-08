from distutils.core import setup, Extension
import numpy

xoptions = ['-O3',
            '-fopenmp']

include_dirs = []
library_dirs = []
extra_link_args = []

if 0:
    include_dirs += ['/opt/homebrew/opt/libomp/include']
    library_dirs += ['/opt/homebrew/opt/libomp/lib']
    extra_link_args += ['-lomp']
    xoptions += ['-Xpreprocessor']

pauli_lcu_module = Extension('pauli_lcu_module',
                          sources=['src/c/python_bindings.c'],
                          define_macros=[('NPY_NO_DEPRECATED_API', '7')],
                          include_dirs=[numpy.get_include()] + include_dirs,
                          library_dirs=library_dirs,
                          extra_compile_args=xoptions,
                          extra_link_args=extra_link_args
                          )

setup(ext_modules=[pauli_lcu_module])