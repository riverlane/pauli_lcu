from distutils.core import setup, Extension
import numpy

xoptions = ['-O3']

pauli_lcu_module = Extension('pauli_lcu_module',
                          sources=['src/c/python_bindings.c'],
                          define_macros=[('NPY_NO_DEPRECATED_API', '7')],
                          include_dirs=[numpy.get_include()],
                          extra_compile_args = xoptions
                          )

setup(ext_modules=[pauli_lcu_module])
