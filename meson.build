project('shtools', 'c', 'fortran',
  version: run_command(find_program('python3'), '-c', 'import os; os.chdir("pyshtools"); import _version; print(_version.get_versions()["version"])').stdout().strip(),
  license: 'BSD-3-Clause',
  meson_version: '>=1.1.0',
  default_options: [
    'buildtype=release',
  ])

cc = meson.get_compiler('c')

main_incdir = include_directories('include')

blas = get_option('blas')
if blas == ''
  blas = dependency('openblas', 'OpenBLAS', 'BLAS')
else
  blas = dependency(blas)
endif
lapack = get_option('lapack')
if lapack == ''
  lapack = dependency('openblas', 'OpenBLAS', 'lapack')
else
  lapack = dependency(lapack)
endif

subdir('src')
subdir('pyshtools')
subdir('examples')
