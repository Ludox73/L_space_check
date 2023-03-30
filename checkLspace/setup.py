from setuptools import setup

setup(
    name = 'checkLspace',
    version = '0.1',
    install_requires = ['pandas', 'ast'],
    packages = ['checkLspace'],
    package_dir = {'checkLspace':'src'},
)
