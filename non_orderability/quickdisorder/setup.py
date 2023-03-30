from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os, glob

c_code = glob.glob('src/c/*.c')
cython_code = glob.glob('src/cython/*.pyx')
headers = ['headers/c', 'headers/cython']


ext_modules = [Extension(name = 'quickdisorder.sl2matrix',
                         sources = c_code + cython_code,
                         include_dirs = headers,
                         extra_compile_args = ['-O3'], 
                     )]

setup(
    name = 'quickdisorder',
    version = '0.3',
    install_requires = ['snappy'],
    packages = ['quickdisorder'],
    package_dir = {'quickdisorder':'src/python'},
    ext_modules = cythonize(ext_modules, include_path=headers)
)
