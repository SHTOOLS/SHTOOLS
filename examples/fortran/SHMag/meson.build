TestSHMag = executable(
  'SHMag',
  'SHMag.f95',
  link_with: shtools)

test(
  'SHMag',
  TestSHMag,
  args: [
    meson.current_source_dir() / '../../ExampleDataFiles',
  ])
