jupyter = find_program('jupyter')

ipynb_files = [
  'grids-and-coefficients.ipynb',
  'localized-spectral-analysis.ipynb',
  'gravity-and-magnetic-fields.ipynb',
  'plotting-maps.ipynb',
  'low-level-spherical-harmonic-analyses.ipynb',
  'advanced-localized-spectral-analysis.ipynb',
  'advanced-shcoeffs-and-shgrid-usage.ipynb',
  'spherical-harmonic-normalizations.ipynb',
  'advanced-shwindow-usage.ipynb',
  '3d-plots.ipynb',
]

foreach file : ipynb_files
  test_name = file.split('.')[0]
  test(
    test_name,
    jupyter,
    args: [
      'execute',
      meson.current_source_dir() / file
    ],
    timeout: 0,
    suite: 'python-notebooks')
endforeach
