#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function
import os
import re

try:
    import numpy  # @UnusedImport # NOQA
except:
    msg = ("No module named numpy. Install numpy before SHTOOLS")
    raise ImportError(msg)

from setuptools import find_packages
from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build import build
from subprocess import CalledProcessError, check_output, check_call
from multiprocessing import cpu_count


def get_version():
    """Get version from git and VERSION file

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
            raise RuntimeError('Unable to get version number from git tags')

        # PEP440 compatibility
        if '-' in git_version:
            # increase version by 0.1 if any new revision exists in repo
            version = '{:.1f}'.format(float(version) + 0.1)
            git_revision = check_output(['git', 'rev-parse', 'HEAD'])
            git_revision = git_revision.strip().decode('ascii')
            version += '.dev0+' + git_revision[:7]

    return version


# Custom Builder
class SHTOOLS_build(build):
    description = "builds python documentation"

    def run(self):
        # build Fortran library using the makefile
        print('---- BUILDING FORTRAN ----')
        make_fortran = ['make', 'fortran']

        try:
            make_fortran.append('-j%d' % cpu_count())
        except:
            pass

        # build python module
        print('---- BUILDING PYTHON ----')
        check_call(make_fortran)
        build.run(self)

        # build documentation
        print('---- BUILDING DOCS ----')
        docdir = os.path.join(self.build_lib, 'pyshtools', 'doc')
        self.mkpath(docdir)
        doc_builder = os.path.join(self.build_lib, 'pyshtools', 'make_docs.py')
        doc_source = '.'
        check_call(['python', doc_builder, doc_source, self.build_lib])

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


KEYWORDS = ['Spherical Harmonics', 'Wigner Symbols']


INSTALL_REQUIRES = [
    'future>=0.12.4',
    'numpy>=1.0.0',
    'setuptools']

# configure python extension to be compiled with f2py

# Absoft f95 flags:
# F95FLAGS = ['m64', 'O3', 'YEXT_NAMES=LCS', 'YEXT_SFX=_', 'fpic',
#             'speed_math=10']
# gfortran flags:
F95FLAGS = ['-m64', '-fPIC', '-O3', '-ffast-math']


ext1 = Extension(name='pyshtools._SHTOOLS',
                 include_dirs=['modules'],
                 libraries=['SHTOOLS', 'fftw3', 'm', 'lapack', 'blas'],
                 library_dirs=['/usr/local/lib', 'lib'],
                 extra_link_args=F95FLAGS,
                 extra_compile_args=F95FLAGS,
                 sources=['src/pyshtools.pyf', 'src/PythonWrapper.f95'])


ext2 = Extension(name='pyshtools._constant',
                 extra_link_args=F95FLAGS,
                 extra_compile_args=F95FLAGS,
                 sources=['src/PlanetsConstants.f95'])


metadata = dict(
    name='pyshtools',
    version=get_version(),
    description='SHTOOLS - Tools for working with spherical harmonics',
    url='http://shtools.ipgp.fr',
    download_url='https://github.com/SHTOOLS/SHTOOLS/zipball/master',
    author='Mark Wieczorek, Matthias Meschede et al.',
    license='BSD',
    keywords=KEYWORDS,
    install_requires=INSTALL_REQUIRES,
    platforms='OS Independent',
    packages=find_packages(),
    package_data={'': ['doc/*.doc', '*.so']},
    include_package_data=True,
    classifiers=CLASSIFIERS,
    ext_modules=[ext1, ext2],
    cmdclass={'build': SHTOOLS_build}
)


setup(**metadata)
