python_sources = [
  '__init__.py',
  '_version.py',
]

pymod = import('python')
python = pymod.find_installation(pure: false)

# NumPy include directory
# This section was copied from SciPy, except for the a) part about `signal`
# since that directory does not exist here, and b) the absolute paths, since we
# don't link to NumPy core libraries.
# The try-except is needed because when things are split across drives on
# Windows, there is no relative path and an exception gets raised. There may be
# other such cases, so add a catch-all and switch to an absolute path. Relative
# paths are needed when for example a virtualenv is placed inside the source
# tree; Meson rejects absolute paths to places inside the source tree.
# For cross-compilation it is often not possible to run the Python interpreter
# in order to retrieve numpy's include directory. It can be specified in the
# cross file instead:
#
#   [properties]
#   numpy-include-dir = /abspath/to/host-pythons/site-packages/numpy/core/include
#
# This uses the path as is, and avoids running the interpreter.
incdir_numpy = meson.get_external_property('numpy-include-dir', 'not-given')
if incdir_numpy == 'not-given'
  incdir_numpy = run_command(python,
    [
      '-c',
      '''import os
import numpy as np
try:
    incdir = os.path.relpath(np.get_include())
except Exception:
    incdir = np.get_include()
print(incdir)'''
    ],
    check: true
  ).stdout().strip()
endif
inc_np = include_directories(incdir_numpy)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'


pyshtools_pyf = custom_target('pyshtools_c',
  input: '../src/pyshtools.pyf',
  output: ['_SHTOOLS-f2pywrappers.f', '_SHTOOLSmodule.c'],
  command: [python, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@'])

pyshtool = python.extension_module('_SHTOOLS',
  pyshtools_pyf, '../src/PythonWrapper.f95', fortranobject_c,
  c_args: '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',
  include_directories: [inc_np, inc_f2py],
  link_with: shtools,
  dependencies: python.dependency(),
  install: true,
  subdir: 'pyshtools')

python.install_sources(python_sources, subdir: 'pyshtools')

subdir('backends')
subdir('constants')
subdir('datasets')
subdir('expand')
subdir('gravmag')
subdir('legendre')
subdir('rotate')
subdir('shclasses')
subdir('shio')
subdir('spectralanalysis')
subdir('utils')
