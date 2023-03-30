from setuptools import setup

setup(
    name = 'foliar',
    version = '0.1',
    install_requires = ['snappy'],
    packages = ['foliar'],
    package_dir = {'foliar':'src'},
)
