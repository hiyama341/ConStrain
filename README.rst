
.. image:: https://raw.githubusercontent.com/hiyama341/ConStrain/main/pictures/constrain1.svg?token=GHSAT0AAAAAABTYCY2KQITX4BYBOV6ISV3UYZRQFDQ
  :width: 400
  :alt: ConStrain logo 

ConStrain: Literate programming can streamline bioengineering workflows
-----------------------------------------------------------------------

.. summary-start

.. image:: https://badge.fury.io/py/ConStrain.svg
        :target: https://badge.fury.io/py/ConStrain

.. image:: https://github.com/hiyama341/ConStrain/actions/workflows/main.yml/badge.svg
        :target: https://github.com/hiyama341/ConStrain/actions

.. image:: https://readthedocs.org/projects/constrain/badge/?version=latest
        :target: https://constrain.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/github/license/hiyama341/ConStrain
        :target: https://github.com/hiyama341/ConStrain/blob/main/LICENSE

.. image:: https://img.shields.io/pypi/pyversions/ConStrain.svg
        :target: https://pypi.org/project/ConStrain/
        :alt: Supported Python Versions

.. image:: https://codecov.io/gh/hiyama341/ConStrain/branch/main/graph/badge.svg?token=P4457QACUY 
        :target: https://codecov.io/gh/hiyama341/ConStrain

.. image:: https://img.shields.io/badge/code%20style-black-black
        :target: https://black.readthedocs.io/en/stable/


What is ConStrain?
~~~~~~~~~~~~~~~~~~

**ConStrain** is an easy-to-use python package with functions that
can be used in literate programming to simulate steps of a strain 
construction cycle from generating genetic parts, to designing a 
combinatorial library along with instructions for the assembly. 
A fully integrated LIMS system is presented to keep track of samples 
and allocation through both a commercial Benchling API and a low-level CSV file database. 

Here, we demonstrate the use of ConStrain in a complex machine learning-guided
metabolic engineering task. We envision that literate programming for biology 
can be adapted for any experimental workflow and be mixed and matched for the 
benefit of the user. As this tool is built to be flexible through its open-source
Python platform, future repetitive tasks can be automated and thus increase 
the speed at which we engineer biology. 

Curious about how you can build strains easier and faster? Head over to our `Google Colab notebooks <https://github.com/hiyama341/ConStrain/tree/main/colab_notebooks>`__
and give it a try.

Please cite our paper (link tba) if you've used ConStrain in a scientific publication.

.. summary-end


Features
--------

* Combinatorial library generation
* HT cloning and transformation workflows
* Flowbot One instructions
* CSV-based LIMS system as well as integration to Benchling
* Genotyping of microbial strains
* Advanced Machine Learning of biological datasets with the AutoML `H2O <https://docs.h2o.ai/h2o/latest-stable/h2o-docs/automl.html>`__
* Workflows for selecting enzyme homologs
* Promoter selection workflows from RNA-seq datasets
* Data analysis of large LC-MS datasets along with workflows for analysis


.. image:: https://github.com/hiyama341/ConStrain/blob/main/pictures/Overview_of_ConStrain.png
  :width: 800
  :alt: Overview of ConStrain's features throughout the DBTL cycle. 

Getting started
~~~~~~~~~~~~~~~
To get started with making microbial strains in an HT manner please follow the steps below: 

1. Install ConStrain. You will find the necessary information below for installation.

2. Check out our `notebooks <https://github.com/hiyama341/ConStrain/tree/main/colab_notebooks>`__ for inspiration to make HT strain construction with ConStrain.

3. You can start making your own workflows by importing ConStrain into either Google colab or Jupyter lab/notebooks.


Installation
~~~~~~~~~~~~

.. installation-start

Use pip to install ConStrain from `PyPI <https://pypi.org/project/ConStrain/>`__.

::

    $ pip install constrain


If you want to develop or if you cloned the repository from our `GitHub <https://github.com/hiyama341/ConStrain/>`__
you can install ConStrain in the following way.

::

    $ pip install -e <path-to-constrain-repo>  


You might need to run these commands with administrative
privileges if you're not using a virtual environment (using ``sudo`` for example).
Please check the `documentation <https://constrain.readthedocs.io/en/latest/installation.html#>`__
for further details.

.. installation-end

Documentation and Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~

Documentation is available on through numerous Google Colab notebooks with
examples on how to use ConStrain and how we use these notebooks for strain
construnction. The Colab notebooks can be found here 
`constrain.notebooks <https://github.com/hiyama341/ConStrain/tree/main/colab_notebooks>`__. 

* Documentation: https://constrain.readthedocs.io.


Contributions
~~~~~~~~~~~~~

Contributions are very welcome! Check our `guidelines <https://constrain.readthedocs.io/en/latest/contributing.html>`__ for instructions how to contribute.


License
~~~~~~~
* Free software: MIT license

Credits
-------
- This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter

.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

- ConStrains logo was made by Jonas Krogh Fischer. Check out his `website <https://github.com/hiyama341/ConStrain/>`__. 