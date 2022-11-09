TestExpandDH = executable(
  'TestExpandDH',
  'TestExpandDH.f95',
  link_with: shtools)

test(
  'ExpandDH',
  TestExpandDH,
  args: [
    meson.current_source_dir() / '../../ExampleDataFiles',
  ])
