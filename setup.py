#!/usr/bin/env python 3
from os import path
from setuptools import setup, find_packages

base_dir = path.dirname(path.realpath(__file__))

with open(path.join(base_dir, 'README.md'), encoding='utf-8') as f:
    LONG_DESC = f.read()

# Get package metadata from 'seqtools/__about__.py' file
base_dir = path.dirname(path.realpath(__file__))
about = {}
with open(path.join(base_dir, 'seqtools', '__about__.py'), encoding='utf-8') as f:
    exec(f.read(), about)

setup(
    name=about['__title__'],

    version=about['__version__'],
    
    author=about['__author__'],
    
    author_email=about['__email__'],
    
    description=about['__summary__'],

    long_description=LONG_DESC,

    url=about['__url__'],
    
    license=about['__license__'],

    packages=find_packages(exclude=['test*']),
    
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    
    keywords='seqtools',

    package_data={
        'seqtools': [
            'data/*.csv',
        ]
    },

    install_requires=[
        'pandas',
        ],
    
    entry_points={
        'console_scripts': [
            'seqtools = seqtools.cli:main',
        ],
    },
)
