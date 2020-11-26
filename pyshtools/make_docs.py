"""
This script builds the python documentation from the function signature and the
customized markdown files. The processed documentation is saved as ascii text
files which are loaded on runtime and replace the __doc__ string of the f2py
wrapped functions.
"""
import sys
import os
import re
import textwrap

import _SHTOOLS


def main():
    # ---- input/output folders ----
    docfolder = os.path.abspath(sys.argv[1])
    libfolder = os.path.abspath(sys.argv[2])
    mddocfolder = os.path.join(docfolder, 'src', 'pydoc')
    pydocfolder = os.path.join(libfolder, 'pyshtools', 'doc')
    print('-- searching documentation in folder: {} --'.format(mddocfolder))

    # ---- loop through the f2py _SHTOOLS functions and make docstrings ----
    for name, func in _SHTOOLS.__dict__.items():
        if callable(func):
            try:
                # ---- process and load documentation ----
                # ---- read md file documentation: ----
                fname_mddoc = os.path.join(mddocfolder, 'py' +
                                           name.lower() + '.md')
                if (os.path.isfile(fname_mddoc)):
                    docstring = process_mddoc(fname_mddoc)
                    # ---- save combined docstring in the pydoc folder ----
                    fname_pydoc = os.path.join(pydocfolder, name.lower() +
                                               '.doc')
                    with open(fname_pydoc, 'w') as pydocfile:
                        pydocfile.write(docstring)

            except IOError as msg:
                print(msg)

    # ---- loop through functions that are defined in python ----
    pyfunctions = ['PlmIndex', 'YilmIndexVector']
    for name in pyfunctions:
        try:
            # ---- process and load documentation
            # read md file documentation:
            fname_mddoc = os.path.join(mddocfolder, 'py' + name.lower() +
                                       '.md')
            docstring = process_mddoc(fname_mddoc)
            # ---- save combined docstring in the pydoc folder--
            fname_pydoc = os.path.join(pydocfolder, name.lower() + '.doc')
            with open(fname_pydoc, 'w') as pydocfile:
                pydocfile.write(docstring)

        except IOError as msg:
            print(msg)


# ===== PROCESS MD DOCUMENTATION FILE ====
def process_mddoc(fname_mddoc):
    # ---- md file search patterns ----
    revalue = re.compile('## Value\n\n', re.DOTALL)
    retail = re.compile('# See (.*)', re.DOTALL)
    reh2 = re.compile('## (.*?)\n', re.DOTALL)
    reh1 = re.compile('\A# (.*?)\n', re.DOTALL)
    reh1b = re.compile('\n# (.*?)\n', re.DOTALL)
    recode = re.compile('`(.*?)`', re.DOTALL)
    restaresc = re.compile(r'(\\\*)', re.DOTALL)
    repycodestart = re.compile(r'```python\n', re.DOTALL)
    repycodeend = re.compile(r'```\n', re.DOTALL)
    # rebold = re.compile('(?![\])[*](.*?)(?![\])[*]',re.DOTALL)

    # ---- open md file and search for patterns ----
    with open(fname_mddoc, 'r') as mdfile:
        # remove the first two lines
        mdstring = mdfile.read().split('\n', 2)[2]

    # First, remove '## Value\n\n' from constant documentation
    match = revalue.search(mdstring)
    if match is not None:
        mdstring = re.sub(match.group(0), '', mdstring)

    # Remove python code block around usage StringEnd
    match = repycodestart.search(mdstring)
    if match is not None:
        mdstring = mdstring.replace(match.group(0), '')
    match = repycodeend.search(mdstring)
    if match is not None:
        mdstring = mdstring.replace(match.group(0), '')

    match = retail.search(mdstring)
    if match is not None:
        mdstring = mdstring.replace(match.group(0), '')

    match = reh1.search(mdstring)
    while match is not None:
        mdstring = re.sub(match.group(0), match.group(1) + '\n' +
                          len(match.group(1)) * '-', mdstring)
        match = reh1.search(mdstring)

    match = reh1b.search(mdstring)
    while match is not None:
        mdstring = re.sub(match.group(0), '\n' + match.group(1) + '\n' +
                          len(match.group(1)) * '-', mdstring)
        match = reh1b.search(mdstring)

    match = reh2.search(mdstring)
    while match is not None:
        mdstring = re.sub(match.group(0), match.group(1) + '\n' +
                          len(match.group(1)) * '-', mdstring)
        match = reh2.search(mdstring)

    match = recode.search(mdstring)
    while match is not None:
        mdstring = mdstring.replace(match.group(0), match.group(1))
        match = recode.search(mdstring)

    match = restaresc.search(mdstring)
    while match is not None:
        mdstring = mdstring.replace(match.group(0), '*')
        match = recode.search(mdstring)

    # ---- combine into docstring ----
    docstring = ''
    tmp = mdstring.splitlines(True)

    # --- remove line breaks between parameters ---
    for i in range(0, len(tmp)-3):
        if tmp[i][0:4] == ':   ' and tmp[i+3][0:4] == ':   ':
            tmp[i+1] = ''

    for i in range(0, len(tmp)):
        if tmp[i][0:4] == ':   ':
            docstring += textwrap.fill(tmp[i][4:], width=80,
                                       replace_whitespace=False,
                                       initial_indent='    ',
                                       subsequent_indent='    ') + '\n'
        elif tmp[i] == '':
            pass
        else:
            docstring += textwrap.fill(tmp[i], width=80,
                                       replace_whitespace=False) + '\n'

    return docstring


# ===== PROCESS F2PY DOCUMENTATION ====
def process_f2pydoc(f2pydoc):
    """
    this function replace all optional _d0 arguments with their default values
    in the function signature. These arguments are not intended to be used and
    signify merely the array dimensions of the associated argument.
    """
    # ---- split f2py document in its parts
    # 0=Call Signature
    # 1=Parameters
    # 2=Other (optional) Parameters (only if present)
    # 3=Returns
    docparts = re.split('\n--', f2pydoc)

    if len(docparts) == 4:
        doc_has_optionals = True
    elif len(docparts) == 3:
        doc_has_optionals = False
    else:
        print('-- uninterpretable f2py documentation --')
        return f2pydoc

    # ---- replace arguments with _d suffix with empty string in ----
    # ---- function signature (remove them): ----
    docparts[0] = re.sub('[\[(,]\w+_d\d', '', docparts[0])

    # ---- replace _d arguments of the return arrays with their default value:
    if doc_has_optionals:

        returnarray_dims = re.findall('[\[(,](\w+_d\d)', docparts[3])
        for arg in returnarray_dims:
            searchpattern = arg + ' : input.*\n.*Default: (.*)\n'
            match = re.search(searchpattern, docparts[2])
            if match:
                default = match.group(1)
                docparts[3] = re.sub(arg, default, docparts[3])
                docparts[2] = re.sub(searchpattern, '', docparts[2])

    # ---- remove all optional _d# from optional argument list:
    if doc_has_optionals:
        searchpattern = '\w+_d\d : input.*\n.*Default: (.*)\n'
        docparts[2] = re.sub(searchpattern, '', docparts[2])

    # ---- combine doc parts to a single string
    processed_signature = '\n--'.join(docparts)

    return processed_signature


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
