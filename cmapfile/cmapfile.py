# cmapfile.py

# Copyright (c) 2014-2020, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Write Chimera Map (CMAP) files.

Cmapfile is a Python library and console script to write Chimera Map (CMAP)
files, HDF5 files containing series of 3D XYZ datasets.

CMAP files can be created from numpy arrays and various file formats
containing volume data, e.g. BIN, TIFF, LSM, OIF, and OIB.

CMAP files can be visualized using UCSF Chimera [2], a highly extensible
program for interactive visualization and analysis of molecular structures
and related data.

For command line usage run ``python -m cmapfile --help``

:Author:
  `Christoph Gohlke <https://www.lfd.uci.edu/~gohlke/>`_

:Organization:
  Laboratory for Fluorescence Dynamics. University of California, Irvine

:License: BSD 3-Clause

:Version: 2020.1.1

Requirements
------------
* `CPython >= 3.6 <https://www.python.org>`_
* `Numpy 1.14 <https://www.numpy.org>`_
* `Scipy 1.1 <https://www.scipy.org>`_
* `H5py 2.10 <https://www.h5py.org/>`_
* `Tifffile 2019.1.1 <https://pypi.org/project/tifffile/>`_
* `Oiffile 2020.1.1 <https://pypi.org/project/oiffile/>`_

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
2020.1.1
    Do not write name attribute.
    Remove support for Python 2.7 and 3.5.
    Update copyright.
2018.8.30
    Move cmapfile.py into cmapfile package.
2014.10.10
    Initial release.

"""

__version__ = '2020.1.1'

__all__ = (
    'CmapFile', 'bin2cmap', 'tif2cmap', 'lsm2cmap', 'oif2cmap', 'array2cmap'
)

import sys
import os
import glob
import warnings

import numpy
import h5py

try:
    from ndimage.interpolation import zoom
except ImportError:
    from scipy.ndimage.interpolation import zoom

from tifffile import TiffFile, transpose_axes, natural_sorted, product
from oiffile import OifFile


class CmapFile(h5py.File):
    """Write Chimera MAP formatted HDF5 file."""

    def __init__(self, filename, mode='w', **kwargs):
        """Create new HDF5 file object.

        See h5py.File for parameters.

        """
        h5py.File.__init__(self, name=filename, mode=mode, **kwargs)
        self.mapcounter = 0

    def addmap(self, data, name=None, step=None, origin=None,
               cell_angles=None, rotation_axis=None, rotation_angle=None,
               symmetries=None, astype=None, subsample=16, chunks=True,
               compression=None, verbose=False):
        """Create HDF5 group and datasets according to CMAP format.

        The order of axes is XYZ for 'step', 'origin', 'cell_angles', and
        'rotation_axis'. Data 'shape' and 'chunks' sizes are in ZYX order.

        Parameters
        ----------
        data : array_like
            Map data to store. Must be three dimensional.
        name : str, optional
            Name of map.
        step : sequence of 3 float, optional
            Spacing between samples in XYZ dimensions.
            Chimera defaults to (1.0, 1.0, 1.0)
        origin : sequence of 3 float, optional
            Chimera defaults to (0.0, 0.0, 0.0)
        cell_angles : sequence of 3 float, optional
            Chimera defaults to (90.0, 90.0, 90.0)
        rotation_axis : sequence of 3 float, optional
            Axis to rotate around. Chimera defaults to (0.0, 0.0, 1.0)
        rotation_angle : float, optional
            Extent to rotate around rotation_axis. Chimera defaults to 0.0.
        symmetries : None, optional
            Undocumented.
        astype : numpy dtype, optional
            Datatype of HDF dataset, e.g. 'float32'.
            By default this is data.dtype.
        subsample : int, optional
            Store subsampled datasets up to the specified number (default: 16).
        chunks : bool or sequence of 3 int, optional
            Size of chunks to store in datasets in ZYX order.
            By default HDF5 determines this.
        compression : str, optional
            Type of HDF5 data compression, e.g. None (default) or 'gzip'.
        verbose : bool, optional
            If False (default), do not print messages to stdout.

        """
        data = numpy.atleast_3d(data)
        if data.ndim != 3:
            raise ValueError('map data must be 3 dimensional')
        if astype:
            data = numpy.ascontiguousarray(data, astype)
        else:
            data = numpy.ascontiguousarray(data)
        # create group and write attributes
        group = self.create_group(f'map{self.mapcounter:05d}')

        # do not write name attribute to work around UnicodeDecodeError:
        # File "Chimera\share\VolumeData\cmap\cmap_grid.py", line 21
        # name += ' ' + image_name
        # UnicodeDecodeError: 'utf8' codec can't decode byte 0xf0 in position 1
        # if name:
        #     group.attrs['name'] = name

        if step:
            group.attrs['step'] = step
        if origin:
            group.attrs['origin'] = origin
        if cell_angles:
            group.attrs['cell_angles'] = cell_angles
        if rotation_axis:
            group.attrs['rotation_axis'] = rotation_axis
        if rotation_angle:
            group.attrs['rotation_angle'] = rotation_angle
        if symmetries:
            group.attrs['symmetries'] = symmetries
        # create main dataset
        if verbose:
            print('1 ', end='', flush=True)
        dset = group.create_dataset(f'data{self.mapcounter:05d}',
                                    data=data, chunks=chunks,
                                    compression=compression)
        # create subsampled datasets
        for i, data in enumerate(subsamples(data, int(subsample))):
            sample = 2**(i + 1)
            if verbose:
                print(f'{sample} ', end='', flush=True)
            dset = group.create_dataset(
                f'data{self.mapcounter:05d}_{i + 2}',
                data=data, chunks=chunks, compression=compression)
            dset.attrs['subsample_spacing'] = sample, sample, sample
        self.mapcounter += 1

    def setstep(self, step):
        """Set 'step' attribute on all datasets.

        Parameters
        ----------
        step : sequence of 3 float
            Spacing between samples in XYZ dimensions.

        """
        for name, group in self.items():
            if name.startswith('map'):
                group.attrs.modify('step', step)


def bin2cmap(binfiles, shape, dtype, offset=0, cmapfile=None, fail=True,
             **kwargs):
    r"""Convert series of SimFCS BIN files to Chimera MAP file.

    SimFCS BIN files contain homogeneous data of any type and shape,
    stored C-contiguously in little endian order.
    A common format is: shape=(-1, 256, 256), dtype='uint16'.

    TODO: Support generic strides, storage order, and byteorder.

    Parameters
    ----------
    binfiles : str or sequence of str
        List of BIN file names or file pattern, e.g. '\*.bin'
    shape : tuple of 3 int
        Shape of data in BIN files in ZYX order, e.g. (32, 256, 256).
    dtype : numpy dtype
        Type of data in the BIN files, e.g. 'uint16'.
    offset : int, optional
        Number of bytes to skip at beginning of BIN file (default: 0).
    cmapfile : str, optional
        Name of the output CMAP file. If None (default), the name is
        derived from the first BIN file.
    fail : bool, optional
        If True (default), raise error when reading invalid BIN files.
    kwargs : dict, optional
        Additional parameters passed to the CmapFile.addmap function,
        e.g. verbose, step, origin, cell_angles, rotation_axis,
        rotation_angle, subsample, chunks, and compression.

    """
    binfiles = parse_files(binfiles)
    validate_shape(shape, 3)
    shape = tuple(shape)
    dtype = numpy.dtype(dtype)
    count = product(shape)
    if count < 0:
        count = -1
    if not cmapfile:
        cmapfile = binfiles[0] + '.cmap'
    verbose = kwargs.get('verbose', False)
    if verbose:
        print(f"Creating '{cmapfile}'", flush=True)
    with CmapFile(cmapfile, 'w') as cmap:
        for binfile in binfiles:
            if verbose:
                print('+', os.path.basename(binfile), end=' ', flush=True)
            try:
                with open(binfile, 'rb') as fh:
                    fh.seek(offset)
                    data = numpy.fromfile(fh, dtype=dtype, count=count)
                    data.shape = shape
                    shape = data.shape
            except Exception:
                if fail:
                    raise
                if verbose:
                    print('failed!', end='', flush=True)
                continue
            cmap.addmap(data, name=os.path.basename(binfile), **kwargs)
            if verbose:
                print(flush=True)


def tif2cmap(tiffiles, cmapfile=None, fail=True, **kwargs):
    r"""Convert series of 3D TIFF files to Chimera MAP file.

    Parameters
    ----------
    tiffiles : str or sequence of str
        List of TIFF file names or file pattern, e.g. '\*.tif'.
        Files must contain 3D data of matching shape and dtype.
    cmapfile : str, optional
        Name of the output CMAP file. If None (default), the name is
        derived from the first TIFF file.
    fail : bool, optional
        If True (default), raise error when processing incompatible TIFF files.
    kwargs : dict, optional
        Additional parameters passed to the CmapFile.addmap function,
        e.g. verbose, step, origin, cell_angles, rotation_axis,
        rotation_angle, subsample, chunks, and compression.

    """
    tiffiles = parse_files(tiffiles)
    if not cmapfile:
        cmapfile = tiffiles[0] + '.cmap'
    verbose = kwargs.get('verbose', False)
    if verbose:
        print(f"Creating '{cmapfile}'", flush=True)
    shape = dtype = None
    with CmapFile(cmapfile, 'w') as cmap:
        for tiffile in tiffiles:
            if verbose:
                print('+', os.path.basename(tiffile), end=' ', flush=True)
            try:
                with TiffFile(tiffile) as tif:
                    data = tif.asarray()
                    data = numpy.atleast_3d(numpy.squeeze(data))
                    if not shape:
                        shape = data.shape
                        dtype = data.dtype
                        if len(shape) != 3 or any(i <= 4 for i in shape):
                            raise ValueError('not a 3D map')
                    elif shape != data.shape or dtype != data.dtype:
                        raise ValueError('shape or dtype mismatch')
            except Exception as exc:
                if fail:
                    raise
                if verbose:
                    print(exc, end='', flush=True)
                continue
            cmap.addmap(data, name=os.path.basename(tiffile), **kwargs)
            if verbose:
                print(flush=True)


def lsm2cmap(lsmfile, cmapfile=None, **kwargs):
    """Convert 5D TZCYX LSM file to Chimera MAP files, one per channel.

    Parameters
    ----------
    lsmfile : str
        Name of the LSM file to convert.
    cmapfile : str, optional
        Name of the output CMAP file. If None (default), the name is
        derived from lsmfile.
    kwargs : dict, optional
        Additional parameters passed to the CmapFile.addmap function,
        e.g. verbose, step, origin, cell_angles, rotation_axis,
        rotation_angle, subsample, chunks, and compression.

    """
    verbose = kwargs.get('verbose', False)
    try:
        cmaps = []
        lsm = None
        # open LSM file
        lsm = TiffFile(lsmfile)
        series = lsm.series[0]  # first series contains the image data
        if series.axes != 'TZCYX':
            raise ValueError(
                f'not a 5D LSM file (expected TZCYX, got {series.axes})'
            )
        if verbose:
            print(lsm)
            print(series.shape, series.axes, flush=True)
        # create one CMAP file per channel
        if cmapfile:
            cmapfile = '{}.ch%04d{}'.format(*os.path.splitext(cmapfile))
        else:
            cmapfile = f'{lsmfile}.ch%04d.cmap'
        cmaps = [CmapFile(cmapfile % i) for i in range(series.shape[2])]
        # voxel/step sizes
        if not kwargs.get('step', None):
            try:
                attrs = lsm[0].cz_lsm_info
                kwargs['step'] = (
                    attrs['voxel_size_x'] / attrs['voxel_size_x'],
                    attrs['voxel_size_y'] / attrs['voxel_size_x'],
                    attrs['voxel_size_z'] / attrs['voxel_size_x'])
            except Exception:
                pass
        # iterate over Tiff pages containing data
        pages = iter(series.pages)
        for _ in range(series.shape[0]):  # iterate over time axis
            data = []
            for _ in range(series.shape[1]):  # iterate over z slices
                data.append(next(pages).asarray())
            data = numpy.vstack(data).reshape(series.shape[1:])
            for c in range(series.shape[2]):  # iterate over channels
                # write datasets and attributes
                cmaps[c].addmap(data=data[:, c], **kwargs)
    finally:
        if lsm:
            lsm.close()
        for f in cmaps:
            f.close()


def array2cmap(data, axes, cmapfile, **kwargs):
    """Save numpy ndarray to Chimera MAP files, one per channel.

    Parameters
    ----------
    data : ndarray
        Three to 5 dimensional array.
    axes : str
        Specifies type and order of axes in data array.
        May contain only 'CTZYX'.
    cmapfile : str
        Name of the output CMAP file.
    kwargs : dict, optional
        Additional parameters passed to the CmapFile.addmap function,
        e.g. verbose, step, origin, cell_angles, rotation_axis,
        rotation_angle, subsample, chunks, and compression.

    """
    if len(data.shape) != len(axes):
        raise ValueError('Number of axes do not match data shape')
    data = transpose_axes(data, axes, 'CTZYX')
    try:
        # create one CMAP file per channel
        cmaps = []
        if cmapfile.lower().endswith('.cmap'):
            cmapfile = cmapfile[:-5]
        if data.shape[0] > 1:
            cmaps = [CmapFile(f'{cmapfile}.ch{i:04d}.cmap')
                     for i in range(data.shape[0])]
        else:
            cmaps = [CmapFile(f'{cmapfile}.cmap')]
        # iterate over data and write cmaps
        for c in range(data.shape[0]):  # channels
            for t in range(data.shape[1]):  # times
                cmaps[c].addmap(data=data[c, t], **kwargs)
    finally:
        for f in cmaps:
            f.close()


def oif2cmap(oiffile, cmapfile=None, **kwargs):
    """Convert OIF or OIB files to Chimera MAP files, one per channel.

    Parameters
    ----------
    oiffile : str
        Name of the OIF or OIB file to convert.
    cmapfile : str, optional
        Name of the output CMAP file. If None (default), the name is
        derived from oiffile.
    kwargs : dict, optional
        Additional parameters passed to the CmapFile.addmap function,
        e.g. verbose, step, origin, cell_angles, rotation_axis,
        rotation_angle, subsample, chunks, and compression.

    """
    verbose = kwargs.get('verbose', False)
    with OifFile(oiffile) as oif:
        if verbose:
            print(oif)
        try:
            tiffs = oif.series[0]
        except Exception:
            # oiffile < 2020.1.1
            tiffs = oif.tiffs
        data = tiffs.asarray()
        axes = tiffs.axes + 'YX'
        if verbose:
            print(data.shape, axes, flush=True)
        # voxel/step sizes
        if not kwargs.get('step', None):
            try:
                size = oif_axis_size(oif.mainfile)
                shape = data.shape
                xsize = size['X'] / (shape[-1] - 1)
                kwargs['step'] = (
                    1.0,
                    (size['Y'] / (shape[-2] - 1)) / xsize,
                    (size['Z'] / (shape[axes.index('Z')] - 1)) / xsize)
            except Exception:
                pass
    if cmapfile is None:
        cmapfile = oiffile
    array2cmap(data, axes, cmapfile, **kwargs)


def oif_axis_size(oifsettings):
    """Return dict of axes sizes from OIF main settings."""
    scale = {'nm': 1000.0, 'ms': 1000.0}
    result = {}
    i = 0
    while True:
        try:
            axis = oifsettings[f'Axis {i} Parameters Common']
        except KeyError:
            break
        size = abs(axis['EndPosition'] - axis['StartPosition'])
        size /= scale.get(axis['UnitName'], 1.0)
        result[axis['AxisCode']] = size
        i += 1
    return result


def subsamples(data, maxsample=16, minshape=4):
    """Return iterator over data zoomed by 0.5."""
    # TODO: use faster mipmap or gaussian pyramid generator
    sample = 2
    while sample <= maxsample and all(i >= minshape for i in data.shape):
        data = zoom(data, 0.5, prefilter=False)
        sample *= 2
        yield data


def validate_shape(shape, length=None):
    """Raise ValueError if shape is not a sequence of positive integers."""
    try:
        if length is not None and len(shape) != length:
            raise ValueError()
        if any(i < 1 and i != -1 for i in shape):
            raise ValueError()
    except Exception as exc:
        raise ValueError('invalid shape') from exc


def parse_numbers(numbers, dtype=float, sep=','):
    """Return list of numbers from string of separated numbers."""
    if not numbers:
        return []
    try:
        return [dtype(i) for i in numbers.split(sep)]
    except Exception as exc:
        raise ValueError(f"not a '{sep}' separated list of numbers") from exc


def parse_files(files):
    """Return list of file names from pattern or list of file names.

    Raise ValueError if no files are found.

    """
    #    # list of files as string
    #    if isinstance(files, str):
    #        files = natural_sorted(
    #            match.group(1) or match.group(2)
    #            for match in re.finditer(r'(?:"([^"\t\n\r\f\v]+))"|(\S+)',
    #                                     files))
    try:  # list of files
        if os.path.isfile(files[0]):
            return files
    except Exception:
        pass
    try:  # glob pattern
        files = natural_sorted(glob.glob(files[0]))
        files[0]  # noqa: validation
        return files
    except Exception as exc:
        raise ValueError('no files found') from exc


def main(argv=None):
    """Command line usage main function."""
    if argv is None:
        argv = sys.argv

    import optparse  # TODO: use argparse

    parser = optparse.OptionParser(
        usage='usage: %prog [options] files',
        description='Convert volume data files to Chimera MAP files.',
        version=f'%prog {__version__}', prog='cmapfile')

    opt = parser.add_option
    opt('-q', '--quiet', dest='verbose', action='store_false', default=True)
    opt('--filetype', dest='filetype', default=None,
        help='type of input file(s), e.g. BIN, LSM, OIF, TIF')
    opt('--dtype', dest='dtype', default=None,
        help='type of data in BIN files. e.g. uint16')
    opt('--shape', dest='shape', default=None,
        help='shape of data in BIN files in F order, e.g. 256,256,32')
    opt('--offset', dest='offset', type='int', default=0,
        help='number of bytes to skip at beginning of BIN files')
    opt('--step', dest='step', default=None,
        help='stepsize of data in files in F order, e.g. 1.0,1.0,8.0')
    opt('--cmap', dest='cmap', default=None,
        help='name of output CMAP file')
    opt('--astype', dest='astype', default=None,
        help='type of data in CMAP file. e.g. float32')
    opt('--subsample', dest='subsample', type='int', default=16,
        help='write subsampled datasets to CMAP file')

    options, files = parser.parse_args()
    if not files:
        parser.error('no input files specified')
    try:
        files = parse_files(files)
    except ValueError:
        parser.error('input file not found')
    shape = parse_numbers(options.shape, int)
    if shape and len(shape) != 3:
        parser.error('invalid shape: expected 3 integers')
    shape = tuple(reversed(shape))  # C order
    step = parse_numbers(options.step, float)
    if step and len(step) != 3:
        parser.error('invalid step: expected 3 numbers')
    if options.filetype:
        filetype = options.filetype.upper()
    else:
        filetype = os.path.splitext(files[0])[-1][1:].upper()

    if filetype == 'LSM':
        if len(files) > 1:
            warnings.warn('too many input files')
        lsm2cmap(
            files[0],
            step=step,
            cmapfile=options.cmap,
            astype=options.astype,
            subsample=options.subsample,
            verbose=options.verbose)
    elif filetype in ('OIB', 'OIF'):
        if len(files) > 1:
            warnings.warn('too many input files')
        oif2cmap(
            files[0],
            step=step,
            cmapfile=options.cmap,
            astype=options.astype,
            subsample=options.subsample,
            verbose=options.verbose)
    elif filetype in ('TIF', 'TIFF'):
        tif2cmap(
            files,
            step=step,
            cmapfile=options.cmap,
            astype=options.astype,
            subsample=options.subsample,
            verbose=options.verbose)
    elif filetype == 'CMAP':
        if not step:
            parser.error('no step size specified for CMAP file')
        if options.verbose:
            print(f"Changing step size in '{os.path.basename(files[0])}'",
                  flush=True)
        with CmapFile(files[0], mode='r+') as cmap:
            cmap.setstep(step)
    elif options.dtype and options.shape:
        bin2cmap(
            files,
            dtype=options.dtype,
            shape=shape,
            offset=options.offset,
            step=step,
            cmapfile=options.cmap,
            astype=options.astype,
            subsample=options.subsample,
            verbose=options.verbose)
    else:
        if not options.shape:
            parser.error('no data shape specified')
        if not options.dtype:
            parser.error('no data type specified')
        parser.error(f'do not know how to convert {filetype} to CMAP')
    if options.verbose:
        print('Done.', flush=True)


if __name__ == '__main__':
    sys.exit(main())
