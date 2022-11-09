TestMarsCrustalThickness = executable(
  'MarsCrustalThickness',
  'MarsCrustalThickness.f95',
  link_with: shtools)

test(
  'MarsCrustalThickness',
  TestMarsCrustalThickness,
  args: [
    meson.current_source_dir() / '../../ExampleDataFiles',
    meson.current_source_dir() / 'input.txt',
  ])
