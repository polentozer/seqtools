"""Central place for package metadata."""

# NOTE: We use __title__ instead of simply __name__ since the latter would
#       interfere with a global variable __name__ denoting object's name.
#       Semantic versioning is used. For more information see:
#       https://packaging.python.org/en/latest/distributing/#semantic-versioning-preferred

__title__ = 'seqtools'
__summary__ = 'Command line tools for manipulating `.fasta` files'
__url__ = 'https://github.com/polentozer/seqtools'

__version__ = '0.1.dev'

__author__ = 'Tadej Markus'
__email__ = 'markus.tadej@gmail.com'

__license__ = 'MIT'
__copyright__ = '2018, ' + __author__

__all__ = (
    '__title__',
    '__summary__',
    '__url__',
    '__version__',
    '__author__',
    '__email__',
    '__license__',
    '__copyright__',
)
