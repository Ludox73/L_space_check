from setuptools import setup, Extension
from Cython.Build import cythonize
import os, shutil

regina_root = '/sage/local/'
include_dirs = [regina_root + 'include/regina', '/usr/include/libxml2']
extra_link_args = ['-lregina-engine']

pyregina_ext = Extension('pyregina',
                         sources=['pyregina.pyx'],
                         language='c++',
                         extra_compile_args=['-std=c++11', '-Wno-sign-compare',
                                             '-Wno-reorder'],
                         include_dirs=include_dirs,
                         extra_link_args=extra_link_args)

setup(
    name='pyregina',
    version='0.2',
    ext_modules=cythonize(pyregina_ext)
)

stupid_info = 'pyregina.egg-info'
if os.path.exists(stupid_info):
    shutil.rmtree(stupid_info)
