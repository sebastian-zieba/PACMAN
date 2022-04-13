.. _installation:

Installation
=============================

There are different options to download and install ``PACMAN``:

Using a conda environment
---------------------------------

It is generally recommended to install PACMAN in a separate new environment.
If you do not have Anaconda install yet, you can find instructions here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

If you have that, navigate in a terminal to a directory where you want to store the cloned PACMAN directory and type the following lines:

.. code-block:: console

    conda deactivate
    conda create -n pacman python==3.9
    git clone https://github.com/sebastian-zieba/PACMAN
    cd PACMAN
    conda activate pacman
    pip install .

With these steps we created a conda environment called pacman and installed PACMAN and its dependencies there.

.. note:: In this example, the conda environment was initialized with python v3.9.0. PACMAN was tested with python v3.7 to v3.10, so it is up to the user which version they want to use.

If you want to update your local installation with a new release, navigate to the PACMAN directory again and type:

.. code-block:: console

    git pull
    pip install --upgrade --force-reinstall .

This will reinstall PACMAN and its dependencies.


With pip & GitHub
---------------------------------

You can install ``PACMAN`` using ``pip`` by entering the following line into a terminal:

.. code-block:: console

    pip install git+https://github.com/sebastian-zieba/PACMAN


Directly on GitHub
---------------------------------

1. You can also download from source on `GitHub <https://github.com/sebastian-zieba/PACMAN>`_.
There are two ways to do that:

* On GitHub, you can either click on **Code** and **Download ZIP** followed by unpacking the distribution by opening up a terminal and typing:

.. code-block:: console

    unzip PACMAN-master.zip

* Or Clone the repository using ``git`` by typing:

.. code-block:: console

    git clone https://github.com/sebastian-zieba/PACMAN

2. To install ``PACMAN`` as a package, go into the downloaded PACMAN directory (where setup.cfg is located) and type:

.. code-block:: console

    pip install .


Using pip (PyPI)
---------------------------------

Will be added ...


Tests
=============================

If you want to test if your installation was successful, go back to your cloned PACMAN directory and type:

.. code-block:: console

    pip install .[test]

This will install PACMAN again and the dependencies (including the testing dependencies!).

You can now run pytest by typing:

.. code-block:: console

    pytest tests/tests_all.py -s

The -s flag will also output all print statements during the tests so that you can see what happens.
The tests might need a few minutes (depending on your internet connection speed).

You have passed all tests if you get a message like this in the end:

.. code-block:: console

    =========== 12 passed, 197 warnings in 157.00s (0:02:37) =============
