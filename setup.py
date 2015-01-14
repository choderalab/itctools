try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import subprocess
import os

VERSION = '0.0.1'
__version__ = VERSION
ISRELEASED = False


def git_version():
    # Modified from numpy/setup.py
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        out = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out
    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = 'Unknown'

    return GIT_REVISION


def get_version_info():
    # Modified from numpy/setup.py
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('itctools/version.py'):
        # must be a source distribution, use existing version file
        try:
            from itctools.version import git_revision as GIT_REVISION
        except ImportError:
            raise ImportError("Unable to import git_revision. Try removing "
                              "itctools/version.py and the build directory "
                              "before building.")
    else:
        GIT_REVISION = "Unknown"

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='itctools/version.py'):
    # Modified from numpy/setup.py
    content = """
# This file is generated from itctools setup.py
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()

    versionpy = open(filename, 'w')
    try:
        versionpy.write(content % {'version': VERSION,
                                   'full_version': FULLVERSION,
                                   'git_revision': GIT_REVISION,
                                   'isrelease': str(ISRELEASED)})
    finally:
        versionpy.close()

write_version_py()
setup(
    name='itctools',
    version=__version__,
    package_dir={'itctools': 'itctools'},
    packages=['itctools'],
    url='https://github.com/choderalab/itctools',
    license='',
    author='Bas Rustenburg, John Chodera',
    author_email='bas.rustenburg@choderalab.org',
    description='Tools for setting up ITC experiments in an automated fashion using the Tecan EVO and Auto-iTC 200',
    test_suite='nose.collector',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
    ],
    install_requires=['nose>=1.0', 'pint', 'rednose', 'behave', 'mock', 'coverage', 'openpyxl'],
    setup_requires=['nose>=1.0', 'pint'],
    tests_require=['nose', 'behave', 'rednose', 'pint', 'mock', 'coverage']
)
