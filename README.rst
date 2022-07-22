
=========
ConStrain
=========

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

Installation
~~~~~~~~~~~~

.. installation-start

Use pip to install ConStrain from `PyPI <https://pypi.org/project/ConStrain/>`__.

::

    $ pip install constrain


If you want to develop or you cloned the repository from our `GitHub <https://github.com/hiyama341/ConStrain/>`__
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
~~~~~~~~~~~~~
* Free software: MIT license

Credits
-------
This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter

.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
