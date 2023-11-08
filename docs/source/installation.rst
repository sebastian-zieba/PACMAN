.. _installation:

Installation
=============================

Options
____________________________________________________

There are several options to download and install the most recent, stable release of ``PACMAN`` that is v0.3.1:

Using a conda environment (recommended)
----------------------------------------

We generally recommend installing ``PACMAN`` in a ``conda`` environment.  To install conda, follow these `instructions <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.

Once you have Anaconda installed, navigate in a terminal to a directory where you want to store the cloned ``PACMAN`` directory and type the following lines:

.. code-block:: console

    conda create -n pacman python==3.9
    conda activate pacman
    git clone -b v0.3.1 https://github.com/sebastian-zieba/PACMAN
    cd PACMAN
    pip install -e .

With these steps, we created a conda environment called ``pacman`` and installed ``PACMAN`` and its dependencies there.
We also installed the code in the edit mode with ``-e``, giving you the option to edit the source code.

.. note:: In this example, the conda environment was initialized with python v3.9.0. PACMAN was tested with python v3.8 to v3.11, so it is up to the user which version to use. Python 3.7 reached end of life in `June 2023 <https://devguide.python.org/versions/>`_. We therefore stopped supporting it.

.. note:: ``pip install -e .`` will only install the necessary dependencies. If you also want to install the dependencies to run the tests (using pytest) or work on the docs, you have to use the ``[test]`` or ``[docs]`` arguments respectively after ``pip install -e .``. See `Test Your Installation <https://pacmandocs.readthedocs.io/en/latest/installation.html#test-your-installation>`_ for an example.


If you want to update your local installation with a new release, navigate to the ``PACMAN`` directory again and type:

.. code-block:: console

    git pull
    pip install --upgrade --force-reinstall .

This will reinstall ``PACMAN`` and its dependencies.


With pip & GitHub
---------------------------------

You can install ``PACMAN`` using ``pip`` by entering the following line into a terminal:

.. code-block:: console

    conda create -n pacman python==3.9
    conda activate pacman
    pip install 'pacman@git+https://github.com/sebastian-zieba/PACMAN.git@v0.3.1'


Directly from GitHub
---------------------------------

1. You can also download from source on `GitHub <https://github.com/sebastian-zieba/PACMAN>`_.
There are two ways to do that:

* On GitHub, you can either click on **Code** and **Download ZIP** followed by unpacking the distribution by opening up a terminal and typing:

.. code-block:: console

    unzip PACMAN-master.zip

* Or clone the repository using ``git`` by typing:

.. code-block:: console

    git clone https://github.com/sebastian-zieba/PACMAN

2. To install ``PACMAN`` as a package, go into the downloaded PACMAN directory (where setup.cfg is located) and type:

.. code-block:: console

    pip install -e .


Using pip (PyPI)
---------------------------------

Not implemented yet. Might be added in the future...


Test your installation
____________________________________________________

To test if your installation was successful, navigate to your cloned PACMAN directory and type:

.. code-block:: console

    pip install -e .[test]

The ``[test]`` argument will also install the necessary dependencies to run pytest.

You can now run pytest by typing the following line while (from inside the PACMAN directory):

.. code-block:: console

    pytest tests/tests_all.py -s

The optional ``-s`` flag will also output all print statements during the tests so that you can see what happens.
The tests might take a few minutes (depending on your internet connection speed).

You have passed all tests if you get a message like this in the end:

.. code-block:: console

    =========== 12 passed, 197 warnings in 157.00s (0:02:37) ===========
