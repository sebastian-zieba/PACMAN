.. _example_introduction:

Introduction
============

Let's apply ``PACMAN`` to real observations of GJ 1214!

Fifteen transits of the planet GJ 1214 b were observed in `HST Program GO 13021 <https://archive.stsci.edu/proposal_search.php?mission=hst&id=13021>`_.
In this example, we will reduce and analyze two of these visits to highlight ``PACMAN``'s ability to fit multiple transits simultaneously.

The observation dates of these visits were (YYYY-MM-DD) 2013-03-13 and 2013-03-15.

When the visits are sorted by time, these visits have indices 5 and 6, using zero-based indexing.

.. note::
    STScI uses visit numbers that are not sorted by time
    (see, for example, the `visit information page for HST program 13021 <https://www.stsci.edu/hst-program-info/visits/?program=13021>`_).
    ``PACMAN`` uses a temporally sorted indexing scheme. This program has 15 visits in total. Within ``PACMAN``, the first visit is referred to with index 0, and the last visit has index 14.

Before we download the observations, we will create the directories needed for the ``PACMAN`` analysis: the data directory, the run directory, and the ``pacman_run_files`` directory. Next: :ref:`Before Running <before_running>`.