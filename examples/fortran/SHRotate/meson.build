TestSHRotate = executable(
  'TestSHRotate',
  'TestSHRotate.f95',
  link_with: shtools)

test(
  'SHRotate',
  TestSHRotate,
  args: [
    meson.current_source_dir() / '../../ExampleDataFiles',
    meson.current_source_dir() / 'input.txt',
  ],
  suite: 'fortran')
