#!/usr/bin/env python
from setuptools import setup, find_packages
from glob import glob
import os.path

# Set __version__
exec(open('scriboquant/__init__.py').read())
scripts = ['bin/scRiboQuant']
requires = open("requirements.txt").read().strip().split("\n")

setup(
    name='scRiboQuant',
    version=__version__,
    scripts=scripts,
    packages=find_packages(),
    include_package_data=True,
    install_requires=requires,
    url='https://github.com/vivekbhr/scRiboQuant',
    license='MIT',
    description='A workflow for processing scRibo-seq data',
    zip_safe=False,
    data_files=[("", ["LICENSE"])]
)
