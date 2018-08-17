"""
This script takes the python documentation files, adds a YAML header with a
custom title, strips the first two lines of the file, and outputs a new file
in doc/pages/mydoc/pydoc.
"""
import os
import re

pydocfiles = os.listdir('.')

reh1 = re.compile('\A# (.*?)\n', re.DOTALL)

for fn in pydocfiles:
    name, ext = os.path.splitext(fn)

    if ext == '.md':
        with open(fn, 'r') as mdfile:
            line = mdfile.readline()
            title = reh1.search(line).group(1)
            mdfile.readline()
            doc = mdfile.read()
            doc = doc.replace('# ', '## ')

        string = ('---\n' + 'title: ' + title + ' (Python)\n' +
                  'keywords: spherical harmonics software package, ' +
                  'spherical harmonic transform, legendre functions, ' +
                  'multitaper spectral analysis, fortran, Python, ' +
                  'gravity, magnetic field\n' + 'sidebar: mydoc_sidebar\n' +
                  'permalink: ' + name + '.html\n' + 'summary:\n' +
                  'tags: [python]\n' + 'toc: false\n' + 'editdoc: pydoc\n' +
                  '---\n\n' + doc)

        with open('../../docs/pages/mydoc/pydoc/' + name +
                  '.md', 'w') as www:
            www.write(string)
