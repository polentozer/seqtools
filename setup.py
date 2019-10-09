'''Packaging settings.'''
from codecs import open
from os.path import abspath, dirname, join
from setuptools import find_packages, setup
from seqtools import __version__


this_dir = abspath(dirname(__file__))
with open(join(this_dir, 'README.rtf'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='seqtools',
    version=__version__,
    description='A command line tool for manipulating `.fasta` files',
    long_description=long_description,
    url='https://github.com/polentozer/seqtools',
    author='Tadej Markus',
    author_email='markus.tadej@gmail.com',
    license='MIT',
    classifiers=[
        'Intended Audience :: Developers',
        'Topic :: Utilities',
        'License :: MIT',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='seqtools',
    packages=find_packages(exclude=['docs', 'tests*']),
    package_data={
        'seqtools': [
            'data/*.csv',
            'data/*.spsum',
        ]
    },
    install_requires=['pandas', 'pyperclip'],
    extras_require={
        'test': ['coverage', 'pytest', 'pytest-cov'],
    },
    entry_points={
        'console_scripts': [
            'seqtools=seqtools.cli:main',
        ],
    },
)
