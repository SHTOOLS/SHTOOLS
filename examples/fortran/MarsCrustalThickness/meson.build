TestMarsCrustalThickness = executable(
  'MarsCrustalThickness',
  'MarsCrustalThickness.f95',
  link_with: shtools,
  build_by_default: false,
  install: false)

test(
  'MarsCrustalThickness',
  TestMarsCrustalThickness,
  args: [
    meson.current_source_dir() / '../../ExampleDataFiles',
    meson.current_source_dir() / 'input.txt',
  ],
  suite: 'fortran')
