project('shtools', 'c', 'fortran',
  version: run_command(find_program('python3'), '-m', 'setuptools_scm').stdout().strip(),
  license: 'BSD-3-Clause',
  meson_version: '>=1.0.0')

cc = meson.get_compiler('c')

main_incdir = include_directories('include')

blas = get_option('blas')
if blas == ''
  blas = dependency('openblas', 'BLAS')
else
  blas = dependency(blas)
endif
lapack = get_option('lapack')
if lapack == ''
  lapack = dependency('lapack')
else
  lapack = dependency(lapack)
endif

subdir('src')
subdir('pyshtools')
subdir('examples')
