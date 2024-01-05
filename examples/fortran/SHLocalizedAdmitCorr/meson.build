TestLocalizedAdmitCorr = executable(
  'LocalizedAdmitCorr',
  'LocalizedAdmitCorr.f95',
  link_with: shtools)

test(
  'LocalizedAdmitCorr',
  TestLocalizedAdmitCorr,
  args: [
    meson.current_source_dir() / '../../ExampleDataFiles',
    meson.current_source_dir() / 'input.txt',
  ],
  timeout: 120,
  suite: 'fortran')
