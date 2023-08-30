# cmapfile/setup.py

"""Cmapfile package Setuptools script."""

import re
import sys

from setuptools import setup


def search(pattern, code, flags=0):
    # return first match for pattern in code
    match = re.search(pattern, code, flags)
    if match is None:
        raise ValueError(f'{pattern!r} not found')
    return match.groups()[0]


with open('cmapfile/cmapfile.py', encoding='utf-8') as fh:
    code = fh.read()

version = search(r"__version__ = '(.*?)'", code).replace('.x.x', '.dev')

description = search(r'"""(.*)\.(?:\r\n|\r|\n)', code)

readme = search(
    r'(?:\r\n|\r|\n){2}"""(.*)"""(?:\r\n|\r|\n){2}from __future__',
    code,
    re.MULTILINE | re.DOTALL,
)

readme = '\n'.join(
    [description, '=' * len(description)] + readme.splitlines()[1:]
)

license = search(
    r'(# Copyright.*?(?:\r\n|\r|\n))(?:\r\n|\r|\n)+""',
    code,
    re.MULTILINE | re.DOTALL,
)

license = license.replace('# ', '').replace('#', '')

if 'sdist' in sys.argv:
    with open('LICENSE', 'w', encoding='utf-8') as fh:
        fh.write('BSD 3-Clause License\n\n')
        fh.write(license)
    with open('README.rst', 'w', encoding='utf-8') as fh:
        fh.write(readme)

setup(
    name='cmapfile',
    version=version,
    license='BSD',
    description=description,
    long_description=readme,
    long_description_content_type='text/x-rst',
    author='Christoph Gohlke',
    author_email='cgohlke@cgohlke.com',
    url='https://www.cgohlke.com',
    project_urls={
        'Bug Tracker': 'https://github.com/cgohlke/cmapfile/issues',
        'Source Code': 'https://github.com/cgohlke/cmapfile',
        # 'Documentation': 'https://',
    },
    packages=['cmapfile'],
    package_data={'cmapfile': ['py.typed']},
    python_requires='>=3.9',
    install_requires=[
        'numpy',
        'scipy',
        'h5py',
        'tifffile',
        'oiffile',
    ],
    extras_require={'all': ['imagecodecs']},
    entry_points={'console_scripts': ['cmapfile = cmapfile:main']},
    platforms=['any'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)
