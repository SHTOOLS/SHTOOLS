add_languages('cpp', native: false)

cpptest = executable('cpptest',
  'test.cpp',
  include_directories: main_incdir,
  link_with: shtools,
  build_by_default: false)

# FIXME: Add this when the file actually compiles.
# test('cpptest', cpptest, args: [meson.current_source_dir() / '../ExampleDataFiles'])
