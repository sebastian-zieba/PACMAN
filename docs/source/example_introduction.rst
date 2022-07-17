.. _example_introduction:

Introduction
========================

Let's apply ``PACMAN`` to real observations of GJ 1214.
15 transits of the planet GJ 1214 b were observed in `HST Program GO 13021  <https://archive.stsci.edu/proposal_search.php?mission=hst&id=13021>`_.
In this example, we will reduce and analyze two of these visits to highlight ``PACMAN``'s ability to fit multiple transits simultaneously.

The observation dates of these visits were (YYYY-MM-DD): 2013-03-13,  2013-03-15

If you index visits by time and start with 0, these visits correspond to visit 5 and 6.

.. note::
    STSci uses different visit numbers which are not sorted by time
    (see for example the `visit information page for program 13021 <https://www.stsci.edu/cgi-bin/get-visit-status?id=13021&markupFormat=html&observatory=HST>`_).
    ``PACMAN`` uses an temporally sorted indexing scheme. So in this example with 15 visits in total, the first visit one is refered to with index 0 and the last with index 14.

`In the next step <https://pacmandocs.readthedocs.io/en/latest/astroquery_visits.html>`_, let's download the data using ``astroquery``.
