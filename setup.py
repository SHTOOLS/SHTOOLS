#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Setup script for SHTOOLS."""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function
import os
import re
import sys
import sysconfig
# the setuptools import dummy patches the distutil commands such that
# python setup.py develop works
import setuptools
import numpy

from numpy.distutils.core import setup
from numpy.distutils.command.build import build as _build
from numpy.distutils.command.install import install as _install
from numpy.distutils.command.develop import develop as _develop
from numpy.distutils.fcompiler import FCompiler, get_default_fcompiler
from numpy.distutils.misc_util import Configuration
from numpy.distutils.system_info import get_info, dict_append
from subprocess import CalledProcessError, check_output, check_call


# Convert markdown README.md to restructured text (.rst) for PyPi, and
# remove the first 5 lines that contain a reference to the shtools LOGO.
# pandoc can be installed either by conda or pip:
# conda install -c conda-forge pandoc pypandoc
# pip install pypandoc
try:
    import pypandoc
    rst = pypandoc.convert_file('README.md', 'rst')
    long_description = rst.split('\n', 5)[5]
except(IOError, ImportError):
    print('*** pypandoc is not installed. PYPI description will not be '
          'formatted correctly. ***')
    long_description = open('README.md').read()


# This flag has to be True if the version number indicated in the file
# VERSION has already been released and to False if this is a development
# version of a future release.
ISRELEASED = True


def get_version():
    """Get version from git and VERSION file.

    In the case where the version is not tagged in git, this function appends
    .post0+commit if the version has been released and .dev0+commit if the
    version has not yet been released.

    Derived from: https://github.com/Changaco/version.py
    """
    d = os.path.dirname(__file__)
    # get release number from VERSION
    with open(os.path.join(d, 'VERSION')) as f:
        vre = re.compile('.Version: (.+)$', re.M)
        version = vre.search(f.read()).group(1)

    if os.path.isdir(os.path.join(d, '.git')):
        # Get the version using "git describe".
        cmd = 'git describe --tags'
        try:
            git_version = check_output(cmd.split()).decode().strip()[1:]
        except CalledProcessError:
            print('Unable to get version number from git tags\n'
                  'Setting to x.x')
            git_version = 'x.x'

        # PEP440 compatibility
        if '-' in git_version:
            git_revision = check_output(['git', 'rev-parse', 'HEAD'])
            git_revision = git_revision.strip().decode('ascii')
            # add post0 if the version is released
            # otherwise add dev0 if the version is not yet released
            if ISRELEASED:
                version += '.post0+' + git_revision[:7]
            else:
                version += '.dev0+' + git_revision[:7]

    return version


VERSION = get_version()
print('INSTALLING SHTOOLS {}'.format(VERSION))


# Custom Builder
class build(_build):
    """This overrides the standard build class to include the doc build."""

    description = "builds python documentation"

    def run(self):
        """Build the Fortran library, all python extensions and the docs."""
        print('---- BUILDING ----')
        _build.run(self)

        # build documentation
        print('---- BUILDING DOCS ----')
        docdir = os.path.join(self.build_lib, 'pyshtools', 'doc')
        self.mkpath(docdir)
        doc_builder = os.path.join(self.build_lib, 'pyshtools', 'make_docs.py')
        doc_source = '.'
        check_call([sys.executable, doc_builder, doc_source, self.build_lib])

        print('---- ALL DONE ----')


# Custom Installer
class install(_install):
    """This overrides the standard build class to include the doc build."""

    description = "builds python documentation"

    def run(self):
        """Build the Fortran library, all python extensions and the docs."""
        print('---- CUSTOM INSTALL ----')
        _install.run(self)


# Custom Installer
class develop(_develop):
    """This overrides the standard build class to include the doc build."""

    description = "builds python documentation"

    def run(self):
        """Build the Fortran library, all python extensions and the docs."""
        print('---- CUSTOM DEVELOP ----')
        _develop.run(self)

        # build documentation
        print('---- BUILDING DOCS ----')
        docdir = os.path.join(self.setup_path, 'pyshtools', 'doc')
        self.mkpath(docdir)
        doc_builder = os.path.join(self.setup_path, 'pyshtools',
                                   'make_docs.py')
        doc_source = '.'
        check_call([sys.executable, doc_builder, doc_source, self.setup_path])

        print('---- ALL DONE ----')


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
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: GIS',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics'
]


KEYWORDS = ['Spherical Harmonics', 'Spectral Estimation', 'Slepian Functions',
            'Wigner Symbols', 'Legendre Functions', 'Gravity Field',
            'Magnetic Field']


INSTALL_REQUIRES = [
    'numpy>=' + str(numpy.__version__),
    'scipy>=0.14.0',
    'matplotlib',
    'astropy'
]

# configure python extension to be compiled with f2py


def get_compiler_flags():
    """Set fortran flags depending on the compiler."""
    compiler = get_default_fcompiler()
    if compiler == 'absoft':
        flags = ['-m64', '-O3', '-YEXT_NAMES=LCS', '-YEXT_SFX=_',
                 '-fpic', '-speed_math=10']
    elif compiler == 'gnu95':
        flags = ['-m64', '-fPIC', '-O3', '-ffast-math']
    elif compiler == 'intel':
        flags = ['-m64', '-fpp', '-free', '-O3', '-Tf']
    elif compiler == 'g95':
        flags = ['-O3', '-fno-second-underscore']
    elif compiler == 'pg':
        flags = ['-fast']
    else:
        flags = ['-m64', '-O3']
    return flags


def configuration(parent_package='', top_path=None):
    """Configure all packages that need to be built."""
    config = Configuration('', parent_package, top_path)

    F95FLAGS = get_compiler_flags()

    kwargs = {
        'libraries': [],
        'include_dirs': [],
        'library_dirs': [],
    }
    kwargs['extra_compile_args'] = F95FLAGS
    kwargs['f2py_options'] = ['--quiet']

    # numpy.distutils.fcompiler.FCompiler doesn't support .F95 extension
    compiler = FCompiler(get_default_fcompiler())
    compiler.src_extensions.append('.F95')
    compiler.language_map['.F95'] = 'f90'

    # collect all Fortran sources
    files = os.listdir('src')
    exclude_sources = ['PlanetsConstants.f95', 'PythonWrapper.f95']
    sources = [os.path.join('src', file) for file in files if
               file.lower().endswith(('.f95', '.c')) and file not in
               exclude_sources]

    # (from http://stackoverflow.com/questions/14320220/
    #              testing-python-c-libraries-get-build-path)):
    build_lib_dir = "{dirname}.{platform}-{version[0]}.{version[1]}"
    dirparams = {'dirname': 'temp',
                 'platform': sysconfig.get_platform(),
                 'version': sys.version_info}
    libdir = os.path.join('build', build_lib_dir.format(**dirparams))
    print('searching SHTOOLS in:', libdir)

    # Fortran compilation
    config.add_library('SHTOOLS',
                       sources=sources,
                       **kwargs)

    # SHTOOLS
    kwargs['libraries'].extend(['SHTOOLS'])
    kwargs['include_dirs'].extend([libdir])
    kwargs['library_dirs'].extend([libdir])

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

    config.add_extension('pyshtools._SHTOOLS',
                         sources=['src/pyshtools.pyf',
                                  'src/PythonWrapper.f95'],
                         **kwargs)

    return config


metadata = dict(
    name='pyshtools',
    version=VERSION,
    description='SHTOOLS - Spherical Harmonic Tools',
    long_description=long_description,
    url='https://shtools.github.io/SHTOOLS/',
    download_url='https://github.com/SHTOOLS/SHTOOLS/zipball/master',
    author='The SHTOOLS developers',
    author_email="mark.a.wieczorek@gmail.com",
    license='BSD',
    keywords=KEYWORDS,
    install_requires=INSTALL_REQUIRES,
    platforms='OS Independent',
    packages=setuptools.find_packages(),
    classifiers=CLASSIFIERS,
    configuration=configuration,
    cmdclass={'build': build, 'install': install, 'develop': develop}
)


setup(**metadata)
