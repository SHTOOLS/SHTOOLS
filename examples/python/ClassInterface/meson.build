test(
  'ClassExample',
  python,
  args: files('ClassExample.py'),
  timeout: 120,
  suite: 'python')

test(
  'ExactPower',
  python,
  args: files('exact_power.py'),
  timeout: 60,
  suite: 'python')

test(
  'WindowExample',
  python,
  args: [files('WindowExample.py'),
         meson.current_source_dir() / '../../ExampleDataFiles'],
  timeout: 120,
  suite: 'python')
