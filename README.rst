Write Chimera map (CMAP) files
==============================

Chimera map (CMAP) files are HDF5 files containing series of 3D XYZ datasets.

CMAP files can be visualized using UCSF Chimera [2], a highly extensible
program for interactive visualization and analysis of molecular structures
and related data.

The cmapfile module can create CMAP files from numpy arrays and various file
formats containing volume data, e.g. BIN, TIFF, LSM, OIF, and OIB.

For command line usage run ``python -m cmapfile --help``

:Author:
  `Christoph Gohlke <https://www.lfd.uci.edu/~gohlke/>`_

:Organization:
  Laboratory for Fluorescence Dynamics. University of California, Irvine

:Version: 2018.8.30

Requirements
------------
* `CPython 2.7 or 3.5+ <https://www.python.org>`_
* `Numpy 1.13 <https://www.numpy.org>`_
* `Scipy 1.0 <https://www.scipy.org>`_
* `H5py 2.7 <https://www.h5py.org/>`_
* `Tiffile 2018.8.29 <https://www.lfd.uci.edu/~gohlke/>`_
* `Oiffile 2018.8.29 <https://www.lfd.uci.edu/~gohlke/>`_

References
----------
(1) Thomas Goddard. [Chimera-users] reading in hdf5 files in chimera.
    https://www.cgl.ucsf.edu/pipermail/chimera-users/2008-September/003052.html
(2) UCSF Chimera, an extensible molecular modeling system.
    https://www.cgl.ucsf.edu/chimera/
(3) Globals for Images - SimFCS. https://www.lfd.uci.edu/globals/

Examples
--------
Convert a 5D LSM file to CMAP file::

    python -m cmapfile "/my data directory/large.lsm"

Convert all BIN files in the current directory to test.cmap. The BIN files
are known to contain 128x128x64 samples of 16 bit integers. The CMAP file
will store float32 maps using subsampling up to 16::

    python -m cmapfile --shape 128,128,64 --step 1,1,2 --dtype i2
                       --cmap test.cmap --subsample 16 --astype float32 *.bin

Change the step size in the CMAP file::

    python -m cmapfile --step 1,1,1.5 test.cmap

Notes
-----
The CMAP file format according to [1]::

    Example of HDF format written by Chimera (Chimera map format) follows.
    The Chimera map file reader will allow all fields to be missing (except
    the 3D data).

    /image (group, any name allowed)
     name "centriole" (attribute)
     step (1.2, 1.2, 1.2) (attribute)
     origin (-123.4, -522, 34.5) (attribute)
     cell_angles (90.0, 90.0, 90.0) (attribute)
     rotation_axis (0.0, 0.0, 1.0) (attribute)
     rotation_angle 45.0 (attribute, degrees)
     /data (3d array of uint8 (123,542,82)) (dataset, any name allowed)
     /data_x (3d array of uint8 (123,542,82), alternate chunk shape) (dataset)
     /data_2 (3d array of uint8 (61,271,41)) (dataset, any name allowed)
        subsample_spacing (2, 2, 2) (attribute)
     (more subsampled or alternate chunkshape versions of same data)


Revisions
---------
2018.8.30
    Move cmapfile.py into cmapfile package.
2014.10.10
    Initial release.
