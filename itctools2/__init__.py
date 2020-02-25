"""
itctools2
An itctools repo based on the molssi cookiecutter
"""

# Add imports here
from .itctools import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
