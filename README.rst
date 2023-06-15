=========
breizorro
=========
|Pypi Version|
|Python Versions|
|Project License|

A tool for creating a binary mask given a FITS image

==============
Installation
==============

Installation from source_,
working directory where source is checked out

.. code-block:: bash
  
    $ pip install .

This package is available on *PYPI*, allowing

.. code-block:: bash
  
    $ pip install breizorro

To show help message and exit

.. code-block:: bash
   
   $ breizorro --help

     breizorro.breizorro - 2022-08-24 11:07:39,311 INFO - Welcome to breizorro
     breizorro.breizorro - 2022-08-24 11:07:39,375 INFO - Version: 0.1.1
     breizorro.breizorro - 2022-08-24 11:07:39,375 INFO - Usage: breizorro --help
     usage: breizorro [-h] [-r IMAGE] [-m MASK] [-t THRESHOLD] [-b BOXSIZE]
                      [--savenoise] [--merge MASKs|REGs) [MASK(s|REGs) ...]]
                      [--subtract MASK(s|REGs) [MASK(s|REGs ...]]
                      [--number-islands] [--remove-islands N|COORD [N|COORD ...]]
                      [--ignore-missing-islands]
                      [--extract-islands N|COORD [N|COORD ...]]
                      [--minimum-size MINSIZE] [--make-binary] [--invert]
                      [--dilate R] [--erode N] [--fill-holes] [--sum-peak SUM_PEAK]
                      [-o OUTFILE] [--gui]

     breizorro [options] --restored-image restored_image

     optional arguments:
          -h, --help            show this help message and exit
          -r IMAGE, --restored-image IMAGE
                                Restored image file from which to build mask
          -m MASK, --mask-image MASK
                                Input mask file(s). Either --restored-image or --mask-
                                image must be specfied.
          -t THRESHOLD, --threshold THRESHOLD
                                Sigma threshold for masking (default = 6.5)
          -b BOXSIZE, --boxsize BOXSIZE
                                Box size over which to compute stats (default = 50)
          --savenoise           Enable to export noise image as FITS file (default=do
                                not save noise image)
          --merge MASK(s)|REG(s) [MASK(s)|REG(s) ...]
                                Merge in one or more masks or region files
          --subtract MASK(s)|REG(s) [MASK(s)|REG(s) ...]
                                Subract one or more masks or region files
          --number-islands      Number the islands detected (default=do not number
                                islands)
          --remove-islands N|COORD [N|COORD ...]
                                List of islands to remove from input mask. e.g.
                                --remove-islands 1 18 20 20h10m13s,14d15m20s
          --ignore-missing-islands
                                If an island specified by coordinates does not exist,
                                do not throw an error
          --extract-islands N|COORD [N|COORD ...]
                                List of islands to extract from input mask. e.g.
                                --extract-islands 1 18 20 20h10m13s,14d15m20s
          --minimum-size MINSIZE
                                Remove islands that have areas fewer than or equal to
                                the specified number of pixels
          --make-binary         Replace all island numbers with 1
          --invert              Invert the mask
          --dilate R            Apply dilation with a radius of R pixels
          --erode N             Apply N iterations of erosion
          --fill-holes          Fill holes (i.e. entirely closed regions) in mask
          --sum-peak SUM_PEAK   Sum to peak ratio of flux islands to mask in original
                                image.e.g. --sum-peak 100 will mask everything with a
                                ratio above 100
          -o OUTFILE, --outfile OUTFILE
                                Suffix for mask image (default based on input name
          --gui                 Open mask in gui.

=======
License
=======

This project is licensed under the GNU General Public License v3.0 - see license_ for details.

=============
Contribute
=============

Contributions are always welcome! Please ensure that you adhere to our coding
standards pep8_.

.. |Project License| image:: https://img.shields.io/badge/license-GPL-blue.svg
                     :target: https://github.com/ratt-ru/breizorro/blob/main/LICENSE
                     :alt:

.. |Python Versions| image:: https://img.shields.io/pypi/pyversions/breizorro.svg
                     :target: https://pypi.python.org/pypi/breizorro/
                     :alt:

.. |Pypi Version| image:: https://img.shields.io/pypi/v/breizorro.svg
                  :target: https://pypi.python.org/pypi/breizorro
                  :alt:

.. _source: https://github.com/ratt-ru/breizorro
.. _license: https://github.com/ratt-ru/breizorro/blob/main/LICENSE
.. _pep8: https://www.python.org/dev/peps/pep-0008
