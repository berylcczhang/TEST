.. highlight:: rest

.. _sec-filters:

Filters and Filter Data
=======================

SLUG comes with a fairly extensive list of filters, adapted from the list maintained by Charlie Conroy as part of `fsps <https://code.google.com/p/fsps/>`_. However, users may wish to add additional filters, and so the format of the filter list is documented here for convenience.

Filter data is stored in two ASCII text files, ``FILTER_LIST`` and ``allfilters.dat``, which are stored in the ``lib/filters`` directory. The ``FILTER_LIST`` file is an index listing the available filters. In consists of five whitespace-separated columns. The first column is just an numerical index. The second is the name of the filter; this is the name that should be entered in the ``phot_bands`` keyword (see :ref:`ssec-phot-keywords`) to request photometry in that filter. The third and fourth columns the value of :math:`\beta` and :math:`\lambda_c` (the central wavelength) for that filter -- see :ref:`ssec-spec-phot` for definitions. Anything after the fourth column is regarded as a comment, and can be used freely for a description of that filter.

The ``allfilters.dat`` file contains the filter responses. The file contains a series of entires for different filters, each delineated by a header line that begins with ``#``. The order in which filters appear in this file matches that in which they appear in the ``FILTER_LIST``. After the header line, are a series of lines each containing two numbers. The first is the wavelength in Angstrom, and the second is the filter response function at that wavelength.
