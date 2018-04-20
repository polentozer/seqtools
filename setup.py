#!/usr/bin/env python 3
import os
from setuptools import setup

# Get package metadata from 'iCount/__about__.py' file
about = {}
with open(os.path.join(os.getcwd, 'seqtools', '__about__.py'), encoding='utf-8') as f:
    exec(f.read(), about)

setup(
    name=about['__title__'],
    version=about['__summary__'],
    author=about['__author__'],
    author_email=about['__email__'],
    description=about['__summary__'],
    license=about['__license__'],
    classifiers=[
        'Development Status :: Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='seqtools',
    #
    # TODO: check this section below, what is needed
    # 
    install_requires=['pandas'],
)


    # # exclude tests from built/installed package
    # packages=find_packages(exclude=['*.tests', '*.tests.*', 'docs/presentations',
    #                                 'docs/presentations/*']),
    # package_data={
    #     'iCount': [
    #         'examples/*.sh',
    #     ]
    # },

    # extras_require={
    #     'docs': [
    #         'docutils',
    #         'releases',
    #         'sphinx_rtd_theme',
    #     ],
    #     'package': [
    #         'pypandoc'
    #         'twine',
    #         'wheel',
    #     ],
    #     'test': [
    #         'check-manifest',
    #         'pylint>=1.6.4',
    #         'pycodestyle>=2.1.0',
    #         'pydocstyle>=1.0.0',
    #         'pytest-cov',
    #         'readme_renderer',
    #         'coverage>=4.2',
    #     ],
    # },

    # entry_points={
    #     'console_scripts': [
    #         'iCount = iCount.cli:main',
    #     ],
    # },
