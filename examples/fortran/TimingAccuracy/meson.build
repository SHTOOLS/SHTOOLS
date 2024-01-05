timing_tests = {}

foreach timing: ['TimingAccuracyGLQ', 'TimingAccuracyDH', 'TimingAccuracyDHC', 'TimingAccuracyGLQC']
  e = executable(
    timing,
    f'@timing@.f95',
    link_with: shtools)
  timing_tests += {timing: e}
endforeach

BENCHS = {
  'TimingAccuracyGLQ1': ['TimingAccuracyGLQ', 'input1'],
  'TimingAccuracyGLQ2': ['TimingAccuracyGLQ', 'input2'],
  'TimingAccuracyDH1': ['TimingAccuracyDH', 'input3'],
  'TimingAccuracyDH2': ['TimingAccuracyDH', 'input4'],
  'TimingAccuracyDHC1': ['TimingAccuracyDHC', 'input5'],
  'TimingAccuracyDHC2': ['TimingAccuracyDHC', 'input6'],
  'TimingAccuracyGLQC1': ['TimingAccuracyGLQC', 'input7'],
  'TimingAccuracyGLQC2': ['TimingAccuracyGLQC', 'input8'],
}
foreach name, args : BENCHS
  test_name = args[0]
  input = args[1]
  benchmark(
    name,
    timing_tests[test_name],
    args: [
      meson.current_source_dir() / input,
    ],
    suite: 'fortran')
endforeach
