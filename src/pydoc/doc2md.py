"""
Convert docstrings of pure python functions to a simple markdown format
suitable for web documentation.
"""
import pyshtools as pysh

functions = [
    ['pysh.expand.spharm', 'pyspharm.md'],
    ['pysh.expand.spharm_lm', 'pyspharm_lm.md'],
    ['pysh.gravmag.mag_spectrum', 'mag_spectrum.md'],
    ['pysh.legendre.legendre', 'pylegendre.md'],
    ['pysh.legendre.legendre_lm', 'pylegendre_lm.md'],
    ['pysh.legendre.PlmIndex', 'pyplmindex.md'],
    ['pysh.shio.convert', 'convert.md'],
    ['pysh.shio.read_bshc', 'read_bshc.md'],
    ['pysh.shio.write_bshc', 'write_bshc.md'],
    ['pysh.shio.read_dov', 'read_dov.md'],
    ['pysh.shio.write_dov', 'write_dov.md'],
    ['pysh.shio.read_icgem_gfc', 'read_icgem_gfc.md'],
    ['pysh.shio.write_icgem_gfc', 'write_icgem_gfc.md'],
    ['pysh.shio.read_igrf', 'read_igrf.md'],
    ['pysh.shio.shread', 'pyshread.md'],
    ['pysh.shio.shwrite', 'shwrite.md'],
    ['pysh.shio.YilmIndexVector', 'pyyilmindexvector.md'],
    ['pysh.spectralanalysis.cross_spectrum', 'cross_spectrum.md'],
    ['pysh.spectralanalysis.spectrum', 'spectrum.md'],
    ['pysh.utils.figstyle', 'figstyle.md'],
    ['pysh.backends.select_preferred_backend', 'select_preferred_backend.md'],
    ['pysh.backends.preferred_backend', 'preferred_backend.md'],
    ['pysh.backends.preferred_backend_module', 'preferred_backend_module.md'],
    ['pysh.backends.backend_module', 'backend_module.md']
]

for func, file in functions:
    list = func.split(sep='.')
    with open(file, 'w') as f:
        doc = getattr(getattr(getattr(pysh, list[1]), list[2]), '__doc__')
        lines = doc.splitlines(keepends=True)
        f.write('# {}()\n'.format(list[2]))

        kind = [[] for i in range(len(lines))]
        # Classify each line before processing
        for i in range(len(lines)):
            lines[i] = lines[i].replace('*', '\*')
            if lines[i].isspace():
                kind[i] = 'blank'
            elif lines[i].strip() == 'Usage':
                kind[i] = 'header_usage'
            elif lines[i].strip()[0] == '-':
                kind[i] = 'underline'
                if kind[i-1] != 'header_usage':
                    kind[i-1] = 'header'
            elif len(lines[i].strip().split()) == 1:
                if len(lines[i]) > 8:
                    if lines[i][0:8] == '        ' and kind[i-1] != 'blank':
                        kind[i] = 'continuation'
                    else:
                        kind[i] = 'text'
                else:
                    kind[i] = 'text'
            elif len(lines[i].strip().split()) > 1:
                if lines[i].strip().split()[1] == ":":
                    kind[i] = 'parameter'
                elif lines[i].strip().split()[1] == "=":
                    kind[i] = 'usage'
                elif len(lines[i]) > 8:
                    if lines[i][0:9] == '         ':
                        if kind[i-1] == 'parameter':
                            kind[i] = 'continuation'
                        else:
                            kind[i] = 'continuation'
                    elif lines[i][0:8] == '        ':
                        if kind[i-1] == 'parameter':
                            kind[i] = 'desc'
                        elif kind[i-1] == 'continuation' and kind[i-2] == \
                                'parameter':
                            kind[i] = 'desc'
                        else:
                            kind[i] = 'continuation'
                    else:
                        kind[i] = 'text'
                else:
                    kind[i] = 'text'

            else:
                kind[i] = 'text'

            if kind[i] == 'text':
                if kind[i-2] == 'header_usage':
                    kind[i] = 'usage'

        for i in range(len(lines)):
            if kind[i] == 'header' or kind[i] == 'header_usage':
                lines[i] = '# ' + lines[i].strip() + '\n'
            elif kind[i] == 'underline':
                lines[i] = '\n'
            elif kind[i] == 'parameter':
                lines[i] = lines[i].strip() + '\n'
                if lines[i-1] != '\n':
                    lines[i-1] = lines[i-1] + '\n'
            elif kind[i] == 'usage':
                lines[i] = lines[i].strip() + '\n'
            elif kind[i] == 'continuation' and kind[i-1] == 'parameter':
                lines[i] = lines[i].strip() + '\n'
                lines[i-1] = lines[i-1].strip() + ' '
            elif kind[i] == 'continuation' and kind[i-1] == 'usage':
                lines[i] = '    ' + lines[i].strip() + '\n'
                lines[i-1] = lines[i-1].strip() + '\n'
            elif kind[i] == 'desc' and kind[i-1] == 'parameter':
                lines[i] = ':   ' + lines[i].strip() + '\n'
            elif kind[i] == 'desc' and kind[i-1] == 'continuation':
                 lines[i] = ':   ' + lines[i].strip() + '\n'
            elif kind[i] == 'continuation' and kind[i-1] == 'continuation':
                lines[i] = '    ' + lines[i].strip() + '\n'
            elif kind[i] == 'continuation' and kind[i-1] == 'desc':
                lines[i] = '    ' + lines[i].strip() + '\n'
            elif kind[i] == 'text':
                lines[i] = lines[i].strip() + '\n'
            elif kind[i] == 'blank':
                lines[i] = '\n'

        f.writelines(lines)
