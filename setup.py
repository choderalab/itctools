from distutils.core import setup

setup(
    name='itctools',
    version='0.0.1',
    packages=['itctools'],
    url='https://github.com/choderalab/itctools',
    license='',
    author='Bas Rustenburg, John Chodera',
    author_email='bas.rustenburg@choderalab.org',
    description='Tools for setting up ITC experiments in an automated fashion using the Tecan EVO and Auto-iTC 200',
    test_suite='tests',
    classifiers=[
    'Development Status :: 2 - Pre-Alpha',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.3',
    ],
)
