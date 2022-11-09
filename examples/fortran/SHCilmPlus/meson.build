TestCilmPlus = executable(
  'TestCilmPlus',
  'TestCilmPlus.f95',
  link_with: shtools)

test(
  'CilmPlus',
  TestCilmPlus,
  args: [
    meson.current_source_dir() / '../../ExampleDataFiles',
  ])
