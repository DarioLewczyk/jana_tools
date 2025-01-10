# Imports: {{{
from setuptools import setup, find_packages
#}}}

setup(
    name = 'jana_tools',
    version = '1.0.0',
    packages = find_packages(where='src'),
    package_dir = {'': 'src'},
    install_requires = [
        # List package dependencies here: 
        'pymatgen',
        'tqdm',
        'plotly',
        'texttable',
        'numpy',
        'scipy',

    ],
    entry_points = {
        'console_scripts': [
            # Define command line scripts if necessary
        ],
    },
    author = 'Dario C. Lewczyk',
    author_email = 'darlewczyk@gmail.com',
    description = 'A suite of tools for assisting in the analysis of data from JANA.',
    long_description = open('README.md').read(),
    url = 'https://github.com/DarioLewczyk/jana_tools.git',
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires = '>=3.6',
)
