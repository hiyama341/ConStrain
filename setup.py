#!/usr/bin/env python
import re


"""The setup script."""
import versioneer
from setuptools import setup, find_packages


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

install_requires = []
with open("requirements.txt", encoding="utf-8") as f:
    for line in f.readlines():
        install_requires.append(re.split(r"(<|>|=)=", line)[0])

test_requirements = ['pytest>=3', ]

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
