
import sys
import os
import re
from setuptools import setup, find_packages
from distutils.command.build import build as _build
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


class build(_build):
    description = "run make command from setup.py and than usual build"

    def run(self):
        make_fortran = ['make', 'fortran']

        try:
            make_fortran.append('-j%d' % cpu_count())
        except:
            pass

        check_call(make_fortran)

        if sys.version_info.major == 3:
            make_python = ['make', 'python3']
        else:
            make_python = ['make', 'python2']

        check_call(make_python)

        _build.run(self)


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


metadata = dict(
    name='pyshtools',
    version=get_version(),
    description='SHTOOLS - Tools for working with spherical harmonics',
    url='http://shtools.ipgp.fr',
    download_url='https://github.com/SHTOOLS/SHTOOLS/releases',
    author='Mark Wieczorek, Matthias Meschede et al.',
    license='BSD',
    platforms='OS Independent',
    packages=find_packages(),
    package_data={'': ['doc/*.doc', '*.so']},
    include_package_data=True,
    cmdclass={
        'build': build,
    },
    classifiers=CLASSIFIERS,
)

setup(**metadata)
