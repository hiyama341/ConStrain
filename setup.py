#!/usr/bin/env python

"""The setup script."""
import versioneer
from setuptools import setup, find_packages
import re

# Readme & History
with open('README.rst') as readme_file:
    readme = readme_file.read()
with open('HISTORY.rst') as history_file:
    history = history_file.read()

# Requirements
install_requires = ['pydna>=4.0.7',
'pandas>=1.3.0', 
'benchlingapi>=2.1.12',
'numpy>=1.21.0',
'biopython>=1.79',
'python-dotenv>=0.20.0',
'openpyxl>=3.0.9',
'wheel>=0.37.1',
'synbiopython>=0.2.0']

# Test requirements
test_requirements = []
with open("requirements_test.txt", encoding="utf-8") as f:
    for line in f.readlines():
        test_requirements.append(re.split(r"(<|>|=)=", line)[0])

setup(
    author="Lucas Levassor",
    author_email='lucaslevassor@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A Python package for  constructing microbial strains",
    entry_points={
        'console_scripts': [
            'constrain=constrain.cli:main',
        ],
    },
    install_requires=install_requires,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='ConStrain',
    name='ConStrain',
    packages=find_packages(include=['constrain', 'constrain.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/hiyama341/ConStrain',

    ### Change version and put tag to release on PYPI
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
)
