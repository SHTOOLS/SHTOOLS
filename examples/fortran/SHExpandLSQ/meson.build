TestSHExpandLSQ = executable(
  'TestSHExpandLSQ',
  'TestSHExpandLSQ.f95',
  link_with: shtools)

test(
  'SHExpandLSQ',
  TestSHExpandLSQ,
  args: [
    meson.current_source_dir() / '../../ExampleDataFiles',
    meson.current_source_dir() / 'input.txt',
  ],
  suite: 'fortran')
