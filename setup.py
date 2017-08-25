#!/usr/bin/env python

from setuptools import setup

setup(
    name='methylation_pipeline',
    version='0.0.1',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['methylation_pipeline = src.main:main']
    },
    url='https://github.com/bjpop/methylation_pipeline',
    license='LICENSE',
    description='methylation pipeline: run bismark on bisulphite converted DNA', 
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.6.3",
        "pipeline_base == 1.0.0"
    ],
)
