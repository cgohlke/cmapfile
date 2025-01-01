# cmapfile.py

# Copyright (c) 2014-2025, Christoph Gohlke
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
containing volume data, for example, BIN, TIFF, LSM, OIF, and OIB.

CMAP files can be visualized using UCSF Chimera [2], a program for interactive
visualization and analysis of molecular structures and related data.

:Author: `Christoph Gohlke <https://www.cgohlke.com>`_
:License: BSD 3-Clause
:Version: 2025.1.1

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

- `CPython <https://www.python.org>`_ 3.10.11, 3.11.9, 3.12.8, 3.13.1 64-bit
- `NumPy <https://pypi.org/project/numpy/>`_ 2.1.3
- `Scipy <https://pypi.org/project/scipy/>`_ 1.14.1
- `H5py <https://pypi.org/project/h5py/>`_ 3.12.1
- `Tifffile <https://pypi.org/project/tifffile/>`_ 2024.12.12
- `Oiffile <https://pypi.org/project/oiffile/>`_ 2025.1.1

Revisions
---------

2025.1.1

- Improve type hints.
- Support Python 3.13.

2024.8.28

- Fix lsm2cmap with tifffile > 2024.8.24.
- Drop support for Python 3.9 and numpy < 1.24 (NEP29).

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
- Drop support for Python 3.6 (NEP 29).

2020.1.1

- Do not write name attribute.
- Drop support for Python 2.7 and 3.5.
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
    --filetype=FILETYPE   type of input file(s).
                          For example, BIN, LSM, OIF, TIF
    --dtype=DTYPE         type of data in BIN files. For example, uint16
    --shape=SHAPE         shape of data in BIN files in F order.
                          For example, 256,256,32
    --offset=OFFSET       number of bytes to skip at beginning of BIN files
    --step=STEP           stepsize of data in files in F order.
                          For example, 1.0,1.0,8.0
    --cmap=CMAP           name of output CMAP file
    --astype=ASTYPE       type of data in CMAP file. For example, float32
    --subsample=SUBSAMPLE
                          write subsampled datasets to CMAP file

"""

from __future__ import annotations

__version__ = '2025.1.1'

__all__ = [
    '__version__',
    'CmapFile',
    'bin2cmap',
    'tif2cmap',
    'lsm2cmap',
    'oif2cmap',
    'array2cmap',
]

import glob
import os
import sys
import warnings
from typing import TYPE_CHECKING

import h5py
import numpy
from oiffile import OifFile
from scipy.ndimage import zoom
from tifffile import TiffFile, natural_sorted, product, transpose_axes

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Any

    from numpy.typing import ArrayLike, DTypeLike


class CmapFile(h5py.File):  # type: ignore[misc]
    """Write Chimera MAP formatted HDF5 file.

    Parameters:
        filename:
            Name of file to write.
        mode:
            File open mode.
        **kwargs:
            Arguments passed to `h5py.File`.

    """

    mapcounter: int

    def __init__(
        self,
        filename: str | os.PathLike[Any],
        /,
        mode: str = 'w',
        **kwargs: Any,
    ) -> None:
        h5py.File.__init__(self, name=filename, mode=mode, **kwargs)
        self.mapcounter = 0

    def addmap(
        self,
        data: ArrayLike,
        /,
        *,
        name: str | None = None,
        step: Sequence[float] | None = None,
        origin: Sequence[float] | None = None,
        cell_angles: Sequence[float] | None = None,
        rotation_axis: Sequence[float] | None = None,
        rotation_angle: float | None = None,
        color: Sequence[float] | None = None,
        time: int | None = None,
        channel: int | None = None,
        symmetries: Any | None = None,
        astype: DTypeLike | None = None,
        subsample: int = 16,
        chunks: bool = True,
        compression: str | None = None,
        verbose: bool = False,
    ) -> None:
        """Create HDF5 group and datasets according to CMAP format.

        The order of axes is XYZ for `step`, `origin`, `cell_angles`, and
        `rotation_axis`. Data shape and `chunks` sizes are in ZYX order.

        Parameters:
            data:
                Map data to store. Must be three dimensional.
            name:
                Name of map.
            step:
                Spacing between samples in XYZ dimensions.
                Chimera defaults to (1.0, 1.0, 1.0).
            origin:
                Origin in XYZ dimensions.
                Chimera defaults to (0.0, 0.0, 0.0).
            cell_angles:
                Chimera defaults to (90.0, 90.0, 90.0).
            rotation_axis:
                Axis to rotate around. Chimera defaults to (0.0, 0.0, 1.0).
            rotation_angle:
                Extent to rotate around rotation_axis. Chimera defaults to 0.0.
            color:
                RGBA color values.
                For example, (1.0, 1.0, 0, 1.0).
            time:
                Time series frame number.
            channel:
                Multi-channel index.
            symmetries:
                Undocumented.
            astype:
                Datatype of HDF dataset.
                For example, 'float32'.
                By default, this is data.dtype.
            subsample:
                Store subsampled datasets up to specified number.
            chunks:
                Size of chunks to store in datasets in ZYX order.
                By default, the chunk size is determined by HDF5.
            compression:
                Type of HDF5 data compression.
                For example, None or 'gzip'.
            verbose:
                Print messages to stdout.

        """
        if astype:
            data = numpy.ascontiguousarray(data, astype)
        else:
            data = numpy.ascontiguousarray(data)
        if data.ndim != 3:
            raise ValueError('map data must be 3 dimensional')

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
        if color:
            group.attrs['color'] = color
        if time:
            group.attrs['time'] = time
        if channel:
            group.attrs['channel'] = channel
        if symmetries:
            group.attrs['symmetries'] = symmetries
        # create main dataset
        if verbose:
            print('1 ', end='', flush=True)
        dset = group.create_dataset(
            f'data{self.mapcounter:05d}',
            data=data,
            chunks=chunks,
            compression=compression,
        )
        # create subsampled datasets
        for i, subsampled in enumerate(subsamples(data, int(subsample))):
            sample = 2 ** (i + 1)
            if verbose:
                print(f'{sample} ', end='', flush=True)
            dset = group.create_dataset(
                f'data{self.mapcounter:05d}_{i + 2}',
                data=subsampled,
                chunks=chunks,
                compression=compression,
            )
            dset.attrs['subsample_spacing'] = sample, sample, sample
        self.mapcounter += 1

    def setstep(self, step: Sequence[float], /) -> None:
        """Set `step` attribute on all datasets.

        Parameters:
            step: Spacing between samples in XYZ dimensions.

        """
        for name, group in self.items():
            if name.startswith('map'):
                group.attrs.modify('step', step)


def bin2cmap(
    binfiles: Sequence[str | os.PathLike[Any]] | str,
    /,
    shape: tuple[int, ...],
    dtype: DTypeLike,
    offset: int = 0,
    cmapfile: str | os.PathLike[Any] | None = None,
    fail: bool = True,
    **kwargs: Any,
) -> None:
    r"""Convert series of SimFCS BIN files to Chimera MAP file.

    SimFCS BIN files contain homogeneous data of any type and shape,
    stored C-contiguously in little endian order.
    A common format is `shape=(-1, 256, 256), dtype='uint16'`.

    TODO: Support generic strides, storage order, and byteorder.

    Parameters:
        binfiles:
            List of BIN file names or file pattern.
            For example, '\*.bin'
        shape:
            Shape of data in BIN files in ZYX order,
            For example, (32, 256, 256).
        dtype:
            Type of data in BIN files.
            For example, 'uint16'.
        offset:
            Number of bytes to skip at beginning of BIN file.
        cmapfile:
            Name of output CMAP file.
            If None, the name is derived from the first BIN file.
        fail:
            Raise error when reading invalid BIN files.
        **kwargs:
            Arguments passed to :py:meth:`CmapFile.addmap`.

    """
    binfiles_list = parse_files(binfiles)
    validate_shape(shape, 3)
    shape = tuple(shape)
    dtype = numpy.dtype(dtype)
    count = product(shape)
    if count < 0:
        count = -1
    if not cmapfile:
        cmapfile = os.fspath(binfiles_list[0]) + '.cmap'
    verbose = kwargs.get('verbose', False)
    if verbose:
        print(f"Creating '{cmapfile}'", flush=True)
    with CmapFile(cmapfile, 'w') as cmap:
        for binfile in binfiles_list:
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


def tif2cmap(
    tiffiles: Sequence[str | os.PathLike[Any]],
    /,
    cmapfile: str | os.PathLike[Any] | None = None,
    fail: bool = True,
    **kwargs: Any,
) -> None:
    r"""Convert series of 3D TIFF files to Chimera MAP file.

    Parameters:
        tiffiles:
            TIFF file names or file pattern.
            For example, '\*.tif'.
            Files must contain 3D data of matching shape and dtype.
        cmapfile:
            Name of output CMAP file.
            By default, the name is derived from the first TIFF file.
        fail:
            Raise error when processing incompatible TIFF files.
        **kwargs:
            Arguments passed to :py:meth:`CmapFile.addmap`.

    """
    tiffiles_list = parse_files(tiffiles)
    if not cmapfile:
        cmapfile = os.fspath(tiffiles_list[0]) + '.cmap'
    verbose = kwargs.get('verbose', False)
    if verbose:
        print(f"Creating '{cmapfile}'", flush=True)
    shape: tuple[int, ...] | None = None
    dtype = None
    with CmapFile(cmapfile, 'w') as cmap:
        for tiffile in tiffiles_list:
            if verbose:
                print('+', os.path.basename(tiffile), end=' ', flush=True)
            try:
                with TiffFile(tiffile) as tif:
                    data = tif.asarray()
                    data = numpy.atleast_3d(numpy.squeeze(data))
                    if shape is None:
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


def lsm2cmap(
    lsmfile: str | os.PathLike[Any],
    /,
    cmapfile: str | os.PathLike[Any] | None = None,
    **kwargs: Any,
) -> None:
    """Convert 5D TZCYX LSM file to Chimera MAP files, one per channel.

    Parameters:
        lsmfile:
            Name of LSM file to convert.
        cmapfile:
            Name of output CMAP file.
            If None, the name is derived from `lsmfile`.
        **kwargs:
            Arguments passed to :py:meth:`CmapFile.addmap`.

    """
    verbose = kwargs.get('verbose', False)
    try:
        cmaps = []
        lsm = None
        # open LSM file
        lsm = TiffFile(lsmfile)
        series = lsm.series[0]  # first series contains the image data
        if hasattr(series, 'get_shape'):
            # tifffile > 2020.2.25 return squeezed shape and axes
            shape = series.get_shape(False)
            axes = series.get_axes(False)
            if axes[:2] == 'MP' and shape[:2] == (1, 1):
                axes = axes[2:]
                shape = shape[2:]
        else:
            shape = series.shape
            axes = series.axes
        if axes != 'TZCYX':
            if axes == 'ZCYX':
                axes = 'TZCYX'
                shape = (1,) + shape
            else:
                raise ValueError(f'not a 5D LSM file ({axes=} != TZCYX)')
        if verbose:
            print(lsm)
            print(shape, axes, flush=True)
        # create one CMAP file per channel
        if cmapfile:
            cmapfile = '{}.ch%04d{}'.format(*os.path.splitext(cmapfile))
        else:
            cmapfile = f'{lsmfile}.ch%04d.cmap'
        cmaps = [CmapFile(cmapfile % i) for i in range(shape[2])]
        # voxel/step sizes
        if not kwargs.get('step', None):
            try:
                attrs = lsm.lsm_metadata
                if attrs is not None:
                    kwargs['step'] = (
                        attrs['voxel_size_x'] / attrs['voxel_size_x'],
                        attrs['voxel_size_y'] / attrs['voxel_size_x'],
                        attrs['voxel_size_z'] / attrs['voxel_size_x'],
                    )
            except Exception:
                pass
        # iterate over Tiff pages containing data
        pages = iter(series.pages)
        for t in range(shape[0]):  # iterate over time axis
            datalist = []
            for _ in range(shape[1]):  # iterate over z slices
                page = next(pages)
                assert page is not None
                datalist.append(page.asarray())
            data = numpy.vstack(datalist).reshape(shape[1:])
            for c in range(shape[2]):  # iterate over channels
                # write datasets and attributes
                cmaps[c].addmap(data[:, c], time=t, **kwargs)
    finally:
        if lsm is not None:
            lsm.close()
        for f in cmaps:
            f.close()


def array2cmap(
    data: numpy.ndarray[Any, Any],
    /,
    axes: str,
    cmapfile: str | os.PathLike[Any],
    **kwargs: Any,
) -> None:
    """Save numpy ndarray to Chimera MAP files, one per channel.

    Parameters:
        data:
            Three to 5 dimensional array.
        axes:
            Specifies type and order of axes in data array.
            May contain only 'CTZYX'.
        cmapfile:
            Name of output CMAP file.
        **kwargs:
            Arguments passed to :py:meth:`CmapFile.addmap`.

    """
    if len(data.shape) != len(axes):
        raise ValueError('Number of axes do not match data shape')
    data = transpose_axes(data, axes, 'CTZYX')
    try:
        # create one CMAP file per channel
        cmaps = []
        cmapfile = os.fspath(cmapfile)
        if cmapfile.lower().endswith('.cmap'):
            cmapfile = cmapfile[:-5]
        if data.shape[0] > 1:
            cmaps = [
                CmapFile(f'{cmapfile}.ch{i:04d}.cmap')
                for i in range(data.shape[0])
            ]
        else:
            cmaps = [CmapFile(f'{cmapfile}.cmap')]
        # iterate over data and write cmaps
        for c in range(data.shape[0]):  # channels
            for t in range(data.shape[1]):  # times
                cmaps[c].addmap(data[c, t], time=t, **kwargs)
    finally:
        for f in cmaps:
            f.close()


def oif2cmap(
    oiffile: str | os.PathLike[Any],
    /,
    cmapfile: str | os.PathLike[Any] | None = None,
    **kwargs: Any,
) -> None:
    """Convert OIF or OIB files to Chimera MAP files, one per channel.

    Parameters:
        oiffile:
            Name of OIF or OIB file to convert.
        cmapfile:
            Name of output CMAP file.
            By default, the name is derived from `oiffile`.
        **kwargs:
            Arguments passed to :py:meth:`CmapFile.addmap`.

    """
    verbose = kwargs.get('verbose', False)
    with OifFile(oiffile) as oif:
        if verbose:
            print(oif)
        tiffs = oif.series[0]
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
                    (size['Z'] / (shape[axes.index('Z')] - 1)) / xsize,
                )
            except Exception:
                pass
    if cmapfile is None:
        cmapfile = oiffile
    array2cmap(data, axes, cmapfile, **kwargs)


def oif_axis_size(oifsettings: dict[str, Any], /) -> dict[str, Any]:
    """Return mapping of axes sizes from OIF main settings.

    Parameters:
        oifsettings: OIF main settings.

    """
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


def subsamples(
    data: numpy.ndarray[Any, Any], /, maxsample: int = 16, minshape: int = 4
) -> Iterator[numpy.ndarray[Any, Any]]:
    """Return iterator over data zoomed by ~0.5.

    Parameters:
        data:
            Data to be resampled.
        maxsample:
            Inverse of maximum zoom factor.
        minshape:
            Minimum size of any dimension.

    """
    # TODO: use faster mipmap or gaussian pyramid generator
    zoomed = data
    zooms = [1.0 for size in data.shape]
    sample = 2
    while sample <= maxsample and all(i >= minshape for i in zoomed.shape):
        # this formula is used by ChimeraX to calculate subsample zoom factors
        zooms = [((i + sample - 1) // sample) / i for i in data.shape]
        zoomed = zoom(data, zooms, prefilter=False)
        yield zoomed
        sample *= 2


def validate_shape(
    shape: tuple[int, ...], /, length: int | None = None
) -> None:
    """Raise ValueError if shape is not a sequence of positive integers.

    Parameters:
        shape:
            Shape to validate.
        length:
            Expected length of `shape`.

    """
    try:
        if length is not None and len(shape) != length:
            raise ValueError()
        if any(i < 1 and i != -1 for i in shape):
            raise ValueError()
    except Exception as exc:
        raise ValueError('invalid shape') from exc


def parse_numbers(
    numbers: str, /, dtype: type = float, sep: str = ','
) -> list[Any]:
    """Return list of numbers from string of separated numbers.

    Parameters:
        numbers:
            Numbers of type `dtype` separated by `sep.`
        dtype:
            Type of numbers.
        sep:
            Separator used to split numbers.

    """
    if not numbers:
        return []
    try:
        return [dtype(i) for i in numbers.split(sep)]
    except Exception as exc:
        raise ValueError(f'not a {sep!r} separated list of numbers') from exc


def parse_files(
    files: Sequence[str | os.PathLike[Any]], /
) -> Sequence[str | os.PathLike[Any]]:
    """Return list of file names from pattern or list of file names.

    Parameters:
        files: Sequence of file names.

    Raises:
        ValueError: No files are found.

    """
    #    # list of files as string
    #    if isinstance(files, str):
    #        files = natural_sorted(
    #            match.group(1) or match.group(2)
    #            for match in re.finditer(r'(?:"([^"\t\n\r\f\v]+))"|(\S+)',
    try:
        # list of files
        if isinstance(files[0], os.PathLike) or os.path.isfile(files[0]):
            return files
    except Exception:
        pass
    try:
        # glob pattern
        assert isinstance(files[0], str)
        files = natural_sorted(glob.glob(files[0]))
        files[0]  # noqa: validation
        return files
    except Exception as exc:
        raise ValueError('no files found') from exc


def main(argv: list[str] | None = None) -> int:
    """Command line usage main function."""
    if argv is None:
        argv = sys.argv

    import optparse  # TODO: use argparse

    parser = optparse.OptionParser(
        usage='usage: %prog [options] files',
        description='Convert volume data files to Chimera MAP files.',
        version=f'%prog {__version__}',
        prog='cmapfile',
    )

    opt = parser.add_option
    opt('-q', '--quiet', dest='verbose', action='store_false', default=True)
    opt(
        '--filetype',
        dest='filetype',
        default=None,
        help='type of input file(s). For example, BIN, LSM, OIF, TIF',
    )
    opt(
        '--dtype',
        dest='dtype',
        default=None,
        help='type of data in BIN files. For example, uint16',
    )
    opt(
        '--shape',
        dest='shape',
        default=None,
        help='shape of data in BIN files in F order. For example, 256,256,32',
    )
    opt(
        '--offset',
        dest='offset',
        type='int',
        default=0,
        help='number of bytes to skip at beginning of BIN files',
    )
    opt(
        '--step',
        dest='step',
        default=None,
        help='stepsize of data in files in F order. For example, 1.0,1.0,8.0',
    )
    opt('--cmap', dest='cmap', default=None, help='name of output CMAP file')
    opt(
        '--astype',
        dest='astype',
        default=None,
        help='type of data in CMAP file. For example, float32',
    )
    opt(
        '--subsample',
        dest='subsample',
        type='int',
        default=16,
        help='write subsampled datasets to CMAP file',
    )

    options, filesarg = parser.parse_args()
    if not filesarg:
        parser.error('no input files specified')
    try:
        files = parse_files(filesarg)
        if len(files) == 0:
            raise ValueError
    except ValueError:
        parser.error('input file not found')
    shape = tuple(parse_numbers(options.shape, int))
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
            verbose=options.verbose,
        )
    elif filetype in ('OIB', 'OIF'):
        if len(files) > 1:
            warnings.warn('too many input files')
        oif2cmap(
            files[0],
            step=step,
            cmapfile=options.cmap,
            astype=options.astype,
            subsample=options.subsample,
            verbose=options.verbose,
        )
    elif filetype in ('TIF', 'TIFF'):
        tif2cmap(
            files,
            step=step,
            cmapfile=options.cmap,
            astype=options.astype,
            subsample=options.subsample,
            verbose=options.verbose,
        )
    elif filetype == 'CMAP':
        if not step:
            parser.error('no step size specified for CMAP file')
        if options.verbose:
            print(
                f"Changing step size in '{os.path.basename(files[0])}'",
                flush=True,
            )
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
            verbose=options.verbose,
        )
    else:
        if not options.shape:
            parser.error('no data shape specified')
        if not options.dtype:
            parser.error('no data type specified')
        parser.error(f'do not know how to convert {filetype} to CMAP')
    if options.verbose:
        print('Done.', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
