#!/usr/bin/env python
from setuptools import setup
import glob

scripts = glob.glob("*.p*")

setup(
    name='GeneLift',
    version='v1.1',
    description='A tool to order and orient genome assembly contigs via minimap2 alignments to a reference genome.',
    author='Srividya Ramakrishnan',
    author_email='srividya.ramki@gmail.com',
    packages=['geneLift'],
    package_dir={'geneLift': 'geneLift/'},
    install_requires=[
          ],
    scripts=scripts,
    zip_safe=True
)
