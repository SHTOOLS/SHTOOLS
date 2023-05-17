#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Setup script for pyshtools"""

import sys

min_version = (3, 6)

if sys.version_info < min_version:
    error = '\n*** Beginning with pyshtools 4.6, Python {:} or above \n' \
            '*** is required. This error may be a result of using a \n' \
            '*** python-2.7 version of pip. \n' \
            '*** {:}'.format(('.'.join(str(n) for n in min_version)),
                             sys.version_info)
    raise SystemError(error)

import os  # noqa: E402
from pathlib import Path  # noqa: E402
import pkg_resources  # noqa: E402
import sysconfig  # noqa: E402

import setuptools  # noqa: E402
import numpy  # noqa: E402
import versioneer  # noqa: E402
from numpy.distutils.core import setup  # noqa: E402
from numpy.distutils.core import numpy_cmdclass  # noqa: E402
from numpy.distutils.fcompiler import FCompiler  # noqa: E402
from numpy.distutils.fcompiler import get_default_fcompiler  # noqa: E402
from numpy.distutils.misc_util import Configuration  # noqa: E402
from numpy.distutils.system_info import get_info, dict_append  # noqa: E402


CLASSIFIERS = [
    'Development Status :: 5 - Production/Stable',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Fortran',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: GIS',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics'
]

KEYWORDS = ['Spherical Harmonics', 'Spectral Estimation', 'Slepian Functions',
            'Legendre Functions', 'Gravity Field', 'Magnetic Field']

PYTHON_REQUIRES = '>={}'.format('.'.join(str(n) for n in min_version))

INSTALL_REQUIRES = [
    'numpy>=' + str(numpy.__version__),
    'scipy>=0.14.0',
    'matplotlib>=3.3',
    'astropy>=4.0',
    'xarray',
    'requests',
    'pooch>=1.1',
    'tqdm'
]

EXTRAS_REQUIRE = {
    'cartopy': ['cython', 'pyshp', 'six', 'shapely', 'cartopy>=0.18.0'],
    'pygmt': ['pygmt>=0.3'],
    'palettable': ['palettable>=3.3'],
    'ducc': ['ducc0>=0.15']
}

VERSION = versioneer.get_version()
print('INSTALLING SHTOOLS {}'.format(VERSION))


def distutils_dir_name(dname):
    """Returns the name of a distutils build directory"""
    parse_version = pkg_resources.packaging.version.parse
    if parse_version(setuptools.__version__) < parse_version('62.1.0'):
        f = "{dirname}.{platform}-{version[0]}.{version[1]}"
        return f.format(dirname=dname,
                        platform=sysconfig.get_platform(),
                        version=sys.version_info)
    else:
        f = "{dirname}.{platform}-{cache_tag}"
        return f.format(dirname=dname,
                        platform=sysconfig.get_platform(),
                        cache_tag=sys.implementation.cache_tag)


# configure python extension to be compiled with f2py
def configuration(parent_package='', top_path=None):
    """Configure all packages that need to be built."""
    config = Configuration('', parent_package, top_path)

    kwargs = {
        'libraries': [],
        'include_dirs': [],
        'library_dirs': []
    }

    # numpy.distutils.fcompiler.FCompiler doesn't support .F95 extension
    compiler = FCompiler(get_default_fcompiler())
    compiler.src_extensions.append('.F95')
    compiler.language_map['.F95'] = 'f90'

    # collect all Fortran sources
    files = os.listdir('src')
    exclude_sources = ['PlanetsConstants.f95', 'PythonWrapper.f95',
                       'cWrapper.f95']
    sources = [os.path.join('src', file) for file in files if
               file.lower().endswith(('.f95', '.c')) and file not in
               exclude_sources]

    libdir = os.path.join('build', distutils_dir_name('temp'))
    print('searching SHTOOLS in:', libdir)

    # Fortran compilation
    config.add_library('SHTOOLS', sources=sources)

    # SHTOOLS
    kwargs['libraries'].extend(['SHTOOLS'])
    kwargs['include_dirs'].extend([libdir])
    kwargs['library_dirs'].extend([libdir])
    kwargs['f2py_options'] = ['--quiet']

    # FFTW info
    fftw_info = get_info('fftw', notfound_action=2)
    dict_append(kwargs, **fftw_info)

    if sys.platform != 'win32':
        kwargs['libraries'].extend(['m'])

    # BLAS / Lapack info
    lapack_info = get_info('lapack_opt', notfound_action=2)
    blas_info = get_info('blas_opt', notfound_action=2)
    dict_append(kwargs, **blas_info)
    dict_append(kwargs, **lapack_info)

    if sys.platform == 'win32':
        kwargs['runtime_library_dirs'] = []

    config.add_extension('pyshtools._SHTOOLS',
                         sources=['src/pyshtools.pyf',
                                  'src/PythonWrapper.f95'],
                         **kwargs)

    return config


CMDCLASS = numpy_cmdclass
CMDCLASS = versioneer.get_cmdclass(CMDCLASS)

metadata = dict(
    name='pyshtools',
    version=VERSION,
    description='SHTOOLS - Spherical Harmonic Tools',
    long_description=Path('README.md').read_text(encoding='utf-8'),
    long_description_content_type='text/markdown',
    url='https://shtools.github.io/SHTOOLS/',
    download_url='https://github.com/SHTOOLS/SHTOOLS/zipball/master',
    author='The SHTOOLS developers',
    author_email="mark.wieczorek@ipgp.fr",
    license='BSD-3-Clause',
    keywords=KEYWORDS,
    python_requires=PYTHON_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    platforms='OS Independent',
    packages=setuptools.find_packages(),
    package_data={'pyshtools': ['doc/*.doc']},
    classifiers=CLASSIFIERS,
    configuration=configuration,
    cmdclass=CMDCLASS
)


setup(**metadata)
