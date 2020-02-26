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
    tests_require=['nose', 'behave', 'rednose', 'pint', 'mock', 'coverage']
)
"""
itctools2
An itctools repo based on the molssi cookiecutter

import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name='itctools2',
    author='Chodera Lab',
    author_email='john.chodera@choderalab.org',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # python_requires=">=3.5",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,
    )
"""