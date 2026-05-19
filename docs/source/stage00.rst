.. _stage00:

Stage 00
========

.. raw:: html

    <style> .blue {color:blue} </style>

.. topic:: Summary

    - Make sure you followed the steps in :ref:`Before Running <before_running>`.
    - You should now have a run directory containing ``pacman_run_files``.
    - ``pacman_run_files`` should contain ``run_pacman.py``, ``fit_par.txt``, and ``obs_par.pcf``.
    - Make sure the paths to the run directory and data directory are set in ``obs_par.pcf``.
	- Make sure in ``run_pacman.py`` only the line meta = s00.run00(pcf_path=pcf_path) is uncommented and that the other lines (so the other stages) are commented out.
    - Navigate to ``pacman_run_files`` and run ``python run_pacman.py``.
    - Continue with :ref:`Stage 01 <stage01>`.


1) **Set up pcf**

    In :ref:`Before Running <before_running>`, we set up the location of the data directory (``datadir``) and run directory (``rundir``) in ``obs_par.pcf``.

    As a reminder, for this example these directories are:

    - ``datadir = /home/zieba/Desktop/Data/GJ1214_Hubble13021``
    - ``rundir = /home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021``

    The ``pacman_run_files`` directory should contain:

    - ``run_pacman.py``
    - ``fit_par.txt``
    - ``obs_par.pcf``

    These template files can be found in the package directory, on `GitHub <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/pacman_run_files>`_,
    or downloaded directly `here <https://downgit.github.io/#/home?url=https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/pacman_run_files>`_.

    As mentioned in the Introduction, this example analyzes two visits taken in the middle of the GO 13021 program: 2013-03-13 and 2013-03-15.

    **Set up which_visits in the pcf**

    If your ``datadir`` contains many visits, but you only want to analyze a subset, use the ``which_visits`` parameter in ``obs_par.pcf``.

    If you downloaded all 15 visits in GO 13021 (that is 1145 ima files with about 12.5 GB), use:

    .. code-block:: text

        which_visits   [5,6]

    If you only downloaded the two visits as shown in the Download Data instructions, leave the default setting:

    .. code-block:: text

        which_visits   everything


2) **Run PACMAN**

    Navigate to the ``pacman_run_files`` directory in your terminal:

    .. code-block:: console

        cd /home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021/pacman_run_files

    Open ``run_pacman.py``. It should look similar to this (after the imports):

    .. code-block:: python

        if __name__ == "__main__":
            meta = s00.run00(pcf_path=pcf_path)

            # meta = s01.run01(pcf_path=pcf_path)

            # meta = s02.run02(pcf_path=pcf_path)

            # meta = s03.run03(pcf_path=pcf_path)

            # meta = s10.run10(pcf_path=pcf_path)

            # meta = s20.run20(pcf_path=pcf_path)

            # meta = s21.run21(pcf_path=pcf_path)

            # meta = s30.run30(pcf_path=pcf_path)

    To run Stage 00, make sure only the ``s00`` line is uncommented.
    All other stages should remain commented out. In the future you can have multiple stages running at the same time, but for now we will run them one by one.

    Then run:

    .. code-block:: console

        python run_pacman.py

    When Stage 00 is executed, PACMAN creates the ``stage00`` directory and a new timestamped directory inside ``stage00``.
    The directory name has the form:

    .. code-block:: console

        s00_run_YYYY-MM-DD_HH-MM-SS

    For example:

    .. code-block:: console

        /home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021/stage00/s00_run_2022-03-04_15-10-29

    The ``fit_par.txt`` and ``obs_par.pcf`` files used for this run are copied into this stage run directory. This copy is really just a copy and just for future reference. The files in the ``pacman_run_files`` directory are the ones that you should edit to change settings for future runs.

    You should now have a directory structure similar to this:

::

    user
    └── Desktop
        └── Projects
            └── Observations
                └── Hubble
                    └── GJ1214_13021
                        ├── pacman_run_files
                        │   ├── run_pacman.py
                        │   ├── fit_par.txt
                        │   └── obs_par.pcf
                        └── stage00
                            └── s00_run_2022-03-04_15-10-29
                                ├── fit_par.txt
                                ├── obs_par.pcf
                                └── filelist.txt

3) **Results**

    After running Stage 00, you should get terminal output similar to this (in this example I had all 15 visits in the data directory, but only wanted to analyze two of them, so the amount of files analyzed was reduced from 1145 to 158):

    .. code-block:: console

        Starting s00
        Found 1145 data file(s) ending in ima.fits
        Reading in files and their headers: 100%|##########| 1145/1145 [00:03<00:00, 303.42it/s]
        Determining orbit(s) and visit(s): 100%|##########| 1145/1145 [00:00<00:00, 261786.76it/s]
        The user does not want to analyse every visit (which_visits != everything). The amount of files analyzed therefore reduced from 1145 to 158.
        Writing table into filelist.txt
        Saving Metadata
        Finished s00

    Stage 00 creates a file called ``filelist.txt``. It should look like this:

    .. include:: media/s00/filelist.txt
       :literal:

    It has the following columns:

    * ``filenames``

    * ``instr``: The specific filter or grism used in the observation.

    * ``ivisit``: Current visit of the observation.

    * ``iorbit``: Current orbit of the observation.

    * ``t_mjd``: Time in Modified Julian Date (MJD).

    * ``t_visit``: Time since the first exposure in the visit in minutes.

    * ``t_orbit``: Time since the first exposure in the orbit in minutes.

    * ``scan``: Scan direction:

      * ``0``: forward scan

      * ``1``: reverse scan

      * ``-1``: not a spectrum but a direct image

    * ``exp``: exposure time in seconds.