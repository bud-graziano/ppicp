#!/usr/bin/env python

from setuptools import setup

description = 'PPICP - Protein-Protein Interaction Calculation Pipeline'

with open('README.rst', 'r') as f:
    readme = f.read()

setup(
    name='ppicp',
    description=description,
    long_description=readme,
    version='1.0',
    author='Andreas Scheck',
    url='https://github.com/bud-graziano/ppicp',
    platforms=['Linux'],
    license='GNU General Public License v2',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Topic :: Documentation :: Sphinx',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    packages=['ppicp'],
    entry_points={
        "console_scripts": ['ppicp = ppicp.cli:main']
        },
    install_requires=[
        'configparser>=3.5.0',
        'pygal>=2.2.3',
        'requests>=2.10.0'
    ]
)


