Write Chimera Map (CMAP) files
==============================

Cmapfile is a Python library and console script to write Chimera Map (CMAP)
files, HDF5 files containing series of 3D XYZ datasets.

CMAP files can be created from numpy arrays and various file formats
containing volume data, e.g., BIN, TIFF, LSM, OIF, and OIB.

CMAP files can be visualized using UCSF Chimera [2], a program for interactive
visualization and analysis of molecular structures and related data.

:Author: `Christoph Gohlke <https://www.cgohlke.com>`_
:License: BSD 3-Clause
:Version: 2023.8.30

Quickstart
----------

Install the cmapfile package and all dependencies from the
`Python Package Index <https://pypi.org/project/cmapfile/>`_::

    python -m pip install -U cmapfile[all]

Print the command line usage::

    python -m cmapfile --help

See `Examples`_ for usage cases.

Source code and support are available on
`GitHub <https://github.com/cgohlke/cmapfile>`_.

Requirements
------------

This revision was tested with the following requirements and dependencies
(other versions may work):

- `CPython <https://www.python.org>`_ 3.9.13, 3.10.11, 3.11.5, 3.12rc
- `NumPy <https://pypi.org/project/numpy/>`_ 1.25.2
- `Scipy <https://pypi.org/project/scipy/>`_ 1.11.2
- `H5py <https://pypi.org/project/h5py/>`_ 3.9.0
- `Tifffile <https://pypi.org/project/tifffile/>`_ 2023.8.30 (optional)
- `Oiffile <https://pypi.org/project/oiffile/>`_ 2023.8.30 (optional)

Revisions
---------

2023.8.30

- Drop support for Python 3.8 and numpy < 1.22 (NEP29).

2022.9.29

- Make subsampling compatible with ChimeraX (breaking).
- Fix deprecated import of scipy.ndimage.interpolation.zoom.
- Switch to Google style docstrings.

2022.2.2

- Add type hints.
- Drop support for Python 3.7 and numpy < 1.19 (NEP29).

2021.2.26

- Fix LSM conversion with tifffile >= 2021.2.26.
- Remove support for Python 3.6 (NEP 29).

2020.1.1

- Do not write name attribute.
- Remove support for Python 2.7 and 3.5.
- Update copyright.

2018.8.30

- Move cmapfile.py into cmapfile package.

2014.10.10

- Initial release.

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
     color (1.0, 1.0, 0, 1.0) (attribute, rgba 0-1 float)
     time 5 (attribute, time series frame number)
     channel 0 (attribute, integer for multichannel data)
     /data (3d array of uint8 (123,542,82)) (dataset, any name allowed)
     /data_x (3d array of uint8 (123,542,82), alternate chunk shape) (dataset)
     /data_2 (3d array of uint8 (61,271,41)) (dataset, any name allowed)
        subsample_spacing (2, 2, 2) (attribute)
     (more subsampled or alternate chunkshape versions of same data)

References
----------

1. Thomas Goddard. [Chimera-users] reading in hdf5 files in chimera.
   https://www.cgl.ucsf.edu/pipermail/chimera-users/2008-September/003052.html
2. UCSF Chimera, an extensible molecular modeling system.
   https://www.cgl.ucsf.edu/chimera/
3. Globals for Images - SimFCS. https://www.lfd.uci.edu/globals/

Examples
--------

Convert a 5D LSM file to CMAP file::

    $ python -m cmapfile "/my data directory/large.lsm"

Convert all BIN files in the current directory to test.cmap. The BIN files
are known to contain 128x128x64 samples of 16-bit integers. The CMAP file
will store float32 maps using subsampling up to 16::

    $ python -m cmapfile --shape 128,128,64 --step 1,1,2 --dtype i2 \
                         --cmap test.cmap --subsample 16 --astype float32 *.bin

Change the step size in the CMAP file::

    $ python -m cmapfile --step 1,1,1.5 test.cmap

Print the cmapfile script usage::

    $ python -m cmapfile --help

    Usage: cmapfile [options] files

    Convert volume data files to Chimera MAP files.

    Options:
    --version             show program's version number and exit
    -h, --help            show this help message and exit
    -q, --quiet
    --filetype=FILETYPE   type of input file(s), e.g. BIN, LSM, OIF, TIF
    --dtype=DTYPE         type of data in BIN files. e.g. uint16
    --shape=SHAPE         shape of data in BIN files in F order, e.g 256,256,32
    --offset=OFFSET       number of bytes to skip at beginning of BIN files
    --step=STEP           stepsize of data in files in F order, e.g 1.0,1.0,8.0
    --cmap=CMAP           name of output CMAP file
    --astype=ASTYPE       type of data in CMAP file. e.g. float32
    --subsample=SUBSAMPLE
                            write subsampled datasets to CMAP file
