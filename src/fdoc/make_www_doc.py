"""
This script takes the fortran documentation files, adds a YAML header with a
custom title, strips the first two lines of the file, and outputs a new file
in doc/pages/mydoc/fdoc.
"""
import os
import re

fdocfiles = os.listdir('.')

reh1 = re.compile('\A# (.*?)\n', re.DOTALL)

for fn in fdocfiles:
    name, ext = os.path.splitext(fn)

    if ext == '.md' and name != "planetsconstants":
        with open(fn, 'r') as mdfile:
            line = mdfile.readline()
            title = reh1.search(line).group(1)
            mdfile.readline()
            doc = mdfile.read()
            doc = doc.replace('# ', '## ')

        string = ('---\n' + 'title: ' + title + ' (Fortran)\n' +
                  'keywords: spherical harmonics software package, ' +
                  'spherical harmonic transform, legendre functions, ' +
                  'multitaper spectral analysis, fortran, Python, ' +
                  'gravity, magnetic field\n' + 'sidebar: mydoc_sidebar\n' +
                  'permalink: ' + name + '.html\n' + 'summary:\n' +
                  'tags: [fortran]\n' + 'toc: false\n' + 'editdoc: fdoc\n' +
                  '---\n\n' + doc)

        with open('../../docs/pages/mydoc/fdoc/' + name +
                  '.md', 'w') as www:
            www.write(string)
