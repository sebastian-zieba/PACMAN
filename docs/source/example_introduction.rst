.. _example_introduction:

Introduction
========================

Let's apply ``PACMAN`` to real observations of GJ 1214!

15 transits of the planet GJ 1214 b were observed in `HST Program GO 13021  <https://archive.stsci.edu/proposal_search.php?mission=hst&id=13021>`_.
In this example, we will reduce and analyze two of these visits to highlight ``PACMAN``'s ability to fit multiple transits simultaneously.

The observation dates of these visits were (YYYY-MM-DD): 2013-03-13 and 2013-03-15.

When you sort the visits by time, these visits have an index 5 and 6 (if you start indexing with 0).

.. note::
    STSci uses different visit numbers which are not sorted by time
    (see for example the `visit information page for program 13021 <https://www.stsci.edu/cgi-bin/get-visit-status?id=13021&markupFormat=html&observatory=HST>`_).
    ``PACMAN`` uses an temporally sorted indexing scheme. This program has 15 visits in total. So within ``PACMAN``, the first visit is referred to with index 0. Similarly, the last visit has index 14.

Before we download the observations, we will have create some important directories needed for the PACMAN analysis: the "data directory" and the "run directory". So next: :ref:`Before Running <before_running>`.
