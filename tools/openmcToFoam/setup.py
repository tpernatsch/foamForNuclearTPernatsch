# From https://python-packaging-tutorial.readthedocs.io/en/latest/setup_py.html

from setuptools import setup

setup(
    name='MultiGroupXS',
    version='0.1.2',
    author='Thomas Guilbaud',
    author_email='thomas.guilbaud@epfl.ch',
    packages=['MultiGroupXS'],
    # scripts=['bin/MultiGroupXS'],
    # url='http://pypi.python.org/pypi/PackageName/',
    # license='LICENSE.txt',
    description='Package to manage the multi-group XS generation from OpenMC to GeN-Foam',
    # long_description=open('README.txt').read(),
    install_requires=[
        "numpy",
        "openmc",
    ],
)
