=========
ConStrain
=========


.. image:: https://img.shields.io/pypi/v/constrain.svg
        :target: https://pypi.python.org/pypi/constrain

.. image:: https://img.shields.io/travis/hiyama341/constrain.svg
        :target: https://travis-ci.com/hiyama341/constrain

.. image:: https://readthedocs.org/projects/constrain/badge/?version=latest
        :target: https://constrain.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status



ConStrain is an easy-to-use python package with functions that can be used in literate programming to simulate steps of a strain construction cycle from generating genetic parts, to designing a combinatorial library along with instructions for the assembly. A fully integrated LIMS system is presented to keep track of samples and allocation through both a commercial Benchling API and a low-level CSV file database. 

Here, we demonstrate the use of ConStrain in a complex machine learning-guided metabolic engineering task. We envision that literate programming for biology can be adapted for any experimental workflow and be mixed and matched for the benefit of the user. As this tool is built to be flexible through its open-source Python platform, future repetitive tasks can be automated and thus increase the speed at which we engineer biology. 



* Free software: MIT license
* Documentation: https://constrain.readthedocs.io.


Installation
--------

Stable release
To install ConStrain, run this command in your terminal:

    pip install constrain
    
This is the preferred method to install ConStrain, as it will always install the most recent stable release.

If you dont have pip installed, this Python installation guide can guide you through the process.

From sources

The sources for ConStrain can be downloaded from the Github repo.

You can either clone the public repository:

    git clone git://github.com/hiyama341/constrain

Or download the tarball:

    curl -OJL https://github.com/hiyama341/constrain/tarball/master
    
Once you have a copy of the source, you can install it with:

    python setup.py install

Modules
--------

ConStrain consist of three modules to aid in the strain construction process.

1. Design: Helps with planning and excecuting a cloning workflow
2. Lab: Simulation and verification of lab experiments
3. LIMS: A system to keep track of your samples. 

More features are currently being developed to further help with strain construction. 


Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
