.. _download_data:

Download data (alternatives)
=============================

There are different ways to download the data. Let's have a look at them by downloading the data taken in HST Program GO 13021.

    .. note:: PACMAN can currently just work with files with an ``ima`` extension, so you want to select these.
	``ima`` is an intermediate data product standing for calibrated intermediate IR multiaccum image.
	From the `WFC3 data handbook (Types of WFC3 Files) <https://hst-docs.stsci.edu/wfc3dhb/chapter-2-wfc3-data-structure/2-1-types-of-wfc3-files>`_: "For the IR detector, an intermediate MultiAccum (ima) file is the result after all calibrations are applied (dark subtraction, linearity correction, flat fielding, etc.) to all of the individual readouts of the IR exposure."


Using the program's site
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Go to the `new HST MAST search tool <https://mast.stsci.edu/search/hst/ui/#/>`_.

.. image:: media/download/mast1.png

Enter 13021 for Proposal ID and click on search.
All files associated with the program will be then shown in the following table.
One of the columns is called Proposal ID.
You can see that all of them say 13021 because that's what we wanted. Click on any of the numbers saying ``13021``.
You should get to `program site <https://archive.stsci.edu/proposal_search.php?mission=hst&id=13021>`_.
It lists interesting information to the program like papers associated with the program and also (further down on the website) all files taken in this program.

.. image:: media/download/mast2.png

Let's click on ``Mark all`` right above the table. And then on ``Submit marked data for retrieval from STDADS``.
A new page will load where you can finally download your data.

.. image:: media/download/mast3.png

1. Enter your email address or use your STScI log in. (If the data you want to download is still proprietary, you will have to use the log-in)
2. File type: Calibrated data (should be marked already)
3. Extension (bottom right): ima
4. Click Send retrieval request to ST-DADS

You should then get more instructions per email.

You can find more information how to download the files with ftp here:
`Retrieval Options <https://archive.stsci.edu/hst/help/retrieval_help.html>`_,
`Retrieving IUE Data via FTP <https://archive.stsci.edu/iue/ftp_retrieve.html>`_
and `MAST FTP Service <https://archive.stsci.edu/ftp.html>`_.

You can access the data using the ``ftplib`` module in python. Here's an example script which you can also find on GitHub. It will save download the data into the current directory.

.. literalinclude:: media/download/data_downloader_ftps.py
   :language: python

If necessary, unpack everything and move these fits files into a new directory.
This data directory should then contain all the downloaded .fits files. You will have to specify the location of this "data directory" then in the pcf file.

If the user has all files of the 15 visits in the data directory, they can use for example ``which_visits = [13,14]`` in the pcf to only analyze the last two visits.


Using the new MAST search tool
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

An alternative is to download it directly after searching for the program files with the `new HST MAST search tool <https://mast.stsci.edu/search/hst/ui/#/>`_.
Note that the website only allows you to mark a maximum of 100 observations to download.
You can also download the files by using the `old HST MAST search tool <https://archive.stsci.edu/hst/search.php>`_ which does not have this limitation but might be retired soon.

Because the program has more than 1000 files, this would take some time.
So let's only download the last two visits which were taken in August of 2013 (2013-08-12 and 2013-08-20).
Set OBS START DATE = 2013-08-12 and as before, enter 13021 for Proposal ID and click on search.

.. image:: media/download/mast5.png

All files associated with the program taken during the last two visits will be then shown in the following table.
Go to the buttom right and seclect Rows per Page = 100, then select the 100 files and click on DOWNLOAD DATA (100 DATASETS).
A small window opens on the website. Go to the EXTENSIONS drop down list and only select ima. Click on Start DOWNLAÒAD

.. image:: media/download/mast6.png

It will download the ima files associated with these 100 datasets as a zip file.
After that, do the same for the remaining 52 DATASETS.

Finally, unpack everything and move these fits files into a new directory. This data directory should then just contain the downloaded .fits files. You will have to specify the location of this "data directory" then in the pcf file.


