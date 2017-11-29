from distutils.core import setup, Extension
from Cython.Build import cythonize 

ext = Extension(name='cython_wrapper',
		sources=['fdtd3d.c','grid.c','source.c','sensor.c','pml.c','visualize.c','adhie.c','cython_wrapper.pyx'],
		extra_compile_args=['-O3'])

setup(ext_modules=cythonize(ext))




# "cythonize" converts Cython source code to C source code by calling the Cython compiler. 
# The distutils package uses the setup function to compile C source code into an extension module.

# Terminal command: python cython_setup.py build_ext --inplace 

# The "build_ext" argument instructs distutils to build the extension object that cythonize created.
# The optional "--inplace" flag instructs distutils to place each extension module (so) next to its corresponding pyx-file.  

# The default optimization cflag is '-O2'.
# The extra_compile_args are pasted to the end of the command.
# According to the gcc man page, it overrides the default flag:
# "If you use multiple -O options, with or without level numbers, 
#  the last such option is the one that is effective."  

