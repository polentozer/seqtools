from setuptools import setup

setup(
    name="seqtools",
    version="0.1",
    author="Tadej Markus",
    author_email="markus.tadej@gmail.com",
    description=("Command line tool for manipulating `.fasta` files."),
    license="MIT",
    keywords="",
    url="https://github.com/polentozer/seqtools",
    zip_safe=False,  # the package can run out of an .egg file
    classifiers=['Intended Audience :: Science/Research',
                 'Intended Audience :: Developers',
                 'License :: OSI Approved',
                 'Programming Language :: Python',
                 'Topic :: Software Development',
                 'Topic :: Scientific/Engineering',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: Unix',
                 'Operating System :: MacOS'],
    platforms='any',
    scripts=['seqtools'],
    include_package_data=False,
    install_requires=['pandas'],
)