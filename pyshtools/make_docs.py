#!/usr/bin/env python
"""
This script builds the python documentation from the function signature and the
".pod" files that document the Fortran library. The processed documentation is
saved as ascii text files which are loaded on runtime and replace the __doc__
string of the f2py wrapped functions.
"""

import sys
import os
import re
import textwrap
import string
import _SHTOOLS
import _constant


def main():
    #---- input/output folders ----
    libfolder = os.path.abspath(sys.argv[1])
    poddocfolder = os.path.join(libfolder, 'src/pydoc')
    mddocfolder = os.path.join(libfolder, 'src/pydoc')
    pydocfolder = os.path.join(libfolder, 'pyshtools/doc')
    print '---- searching documentation in folder: {} ----'.format(poddocfolder)

    #---- pod file search patterns ----
    redescr = re.compile('DESCRIPTION(.*?)=head1 ARGUMENTS', re.DOTALL)
    renotes = re.compile('NOTES(.*?)=head1', re.DOTALL)

    #---- md file search patterns ----
    retail = re.compile('# See (.*)', re.DOTALL)
    reh2 = re.compile('## (.*?)\n', re.DOTALL)
    reh1 = re.compile('\A# (.*?)\n', re.DOTALL)
    reh1b = re.compile('\n# (.*?)\n', re.DOTALL)
    recode = re.compile('`(.*?)`', re.DOTALL)
    restaresc = re.compile(r'(\\\*)', re.DOTALL)
    #rebold = re.compile('(?![\])[*](.*?)(?![\])[*]',re.DOTALL)

    #---- loop through the f2py _SHTOOLS functions and make docstrings ----
    for name, func in _SHTOOLS.__dict__.items():
        if callable(func):
            try:
                #---- process and load documentation
                # read md file documentation:
                fname_mddoc = os.path.join(mddocfolder, 'py' + name.lower() + '.md')
                if (os.path.isfile(fname_mddoc)):
                    mdfile = open(fname_mddoc, 'r')
                    mdstring = mdfile.read()

                    match = retail.search(mdstring)
                    if match != None:
                        #    mdstring = re.sub(match.group(0),'',mdstring) doesn't work. don't know why
                        mdstring = string.replace(mdstring, match.group(0), '')

                    match = reh1.search(mdstring)
                    while match != None:
                        mdstring = re.sub(match.group(0), match.group(1) + '\n' + len(match.group(1)) * '-' + '\n', mdstring)
                        match = reh1.search(mdstring)

                    match = reh1b.search(mdstring)
                    while match != None:
                        mdstring = re.sub(match.group(0), '\n' + match.group(1) + '\n' + len(match.group(1)) * '-' + '\n', mdstring)
                        match = reh1b.search(mdstring)

                    match = recode.search(mdstring)
                    while match != None:
                        #    mdstring = re.sub(match.group(0),match.group(1),mdstring) doesn't work. don't know why
                        mdstring = string.replace(mdstring, match.group(0), match.group(1))
                        match = recode.search(mdstring)

                    match = restaresc.search(mdstring)
                    while match != None:
                        mdstring = string.replace(mdstring, match.group(0), '*')
                        match = recode.search(mdstring)

                    # wrap each line of the string individually.
                    docstring = ''
                    tmp = mdstring.splitlines(True)
                    for i in range(0, len(tmp)):
                        if tmp[i][0:4] == ':   ':
                            docstring += textwrap.fill(tmp[i][4:], width=80,
                                                       replace_whitespace=False, initial_indent='    ',
                                                       subsequent_indent='    ') + '\n'
                        else:
                            docstring += textwrap.fill(tmp[i], width=80, replace_whitespace=False) + '\n'

                    #---- save combined docstring in the pydoc folder--
                    fname_pydoc = os.path.join(pydocfolder, name.lower() + '.doc')
                    pydocfile = open(fname_pydoc, 'w')
                    pydocfile.write(docstring)
                    pydocfile.close()

                else:
                    # delete this after converted all doc files from pod to md!!!!
                    #-- process and load documentation from different sources --
                    # get and process f2py generated call signature:
                    callsign = process_f2pydoc(func.__doc__)
                    signature = 'SIGNATURE\n---------\n' + callsign
                    # read pod file documentation:
                    fname_poddoc = os.path.join(poddocfolder, name.lower() + '.pod')
                    podfile = open(fname_poddoc, 'r')
                    podstring = podfile.read()
                    match = redescr.search(podstring)
                    if match is None:
                        description = '\nDESCRIPTION\n-----------\nDOCUMENTATION NOT AVAILABLE'
                    else:
                        poddoc = process_pod(match.group(1))
                        description = '\nDESCRIPTION\n-----------' + poddoc
                    match = renotes.search(podstring)
                    if match is None:
                        notes = ''
                    else:
                        poddoc = process_pod(match.group(1))
                        notes = '\nNOTES\n-----' + poddoc
                    podfile.close()
                    # combine docstring:
                    docstring = signature + description + notes
                    #-- save combined docstring in the pydoc folder--
                    fname_pydoc = os.path.join(pydocfolder, name.lower() + '.doc')
                    pydocfile = open(fname_pydoc, 'w')
                    pydocfile.write(docstring)
                    pydocfile.close()

            except IOError, msg:
                print msg

    #---- loop through the f2py constants and make docstrings ----
    for name, value in _constant.planetsconstants.__dict__.items():
        try:
            #---- read md file documentation:
            fname_mddoc = os.path.join(mddocfolder, 'constant_' + name.lower() + '.md')
            mdfile = open(fname_mddoc, 'r')
            mdstring = mdfile.read()

            match = reh1.search(mdstring)
            while match != None:
                mdstring = re.sub(match.group(0), match.group(1) + '\n' + len(match.group(1)) * '=' + '\n', mdstring)
                match = reh1.search(mdstring)

            match = reh2.search(mdstring)
            while match != None:
                mdstring = re.sub(match.group(0), match.group(1) + '\n' + len(match.group(1)) * '-' + '\n', mdstring)
                match = reh2.search(mdstring)

            # wrap each line of the string individually.
            docstring = ''
            tmp = mdstring.splitlines(True)
            for i in range(0, len(tmp)):
                docstring += textwrap.fill(tmp[i], width=80, replace_whitespace=False) + '\n'

            #-- save docstring in the pydoc folder--
            fname_pydoc = os.path.join(pydocfolder, 'constant_' + name.lower() + '.doc')
            pydocfile = open(fname_pydoc, 'w')
            pydocfile.write(docstring)
            pydocfile.close()

        except IOError, msg:
            print msg

#===== PROCESS F2PY DOCUMENTATION ====


def process_f2pydoc(f2pydoc):
    """
    this function replace all optional _d0 arguments with their default values
    in the function signature. These arguments are not intended to be used and
    signify merely the array dimensions of the associated argument.
    """
    #---- split f2py document in its parts
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
        print '-- uninterpretable f2py documentation --'
        return f2pydoc

    #---- replace arguments with _d suffix with empty string in function signature (remove them):
    docparts[0] = re.sub('[\[(,]\w+_d\d', '', docparts[0])

    #---- replace _d arguments of the return arrays with their default value:
    if doc_has_optionals:

        returnarray_dims = re.findall('[\[(,](\w+_d\d)', docparts[3])
        for arg in returnarray_dims:
            searchpattern = arg + ' : input.*\n.*Default: (.*)\n'
            match = re.search(searchpattern, docparts[2])
            if match:
                default = match.group(1)  # this returns the value in brackets in the search pattern
                docparts[3] = re.sub(arg, default, docparts[3])
                docparts[2] = re.sub(searchpattern, '', docparts[2])

    #---- remove all optional _d# from optional argument list:
    if doc_has_optionals:
        searchpattern = '\w+_d\d : input.*\n.*Default: (.*)\n'
        docparts[2] = re.sub(searchpattern, '', docparts[2])

    #---- combine doc parts to a single string
    processed_signature = '\n--'.join(docparts)

    return processed_signature

#===== PROCESS POD DOCUMENTATION ====


def process_pod(string):
    """This function removes some of the pod language commands to get, for example, bold text."""
    while True:
        match = re.search('I<(\w+)>', string)
        if match:
            string = re.sub(match.group(0), match.group(1), string)
        else:
            break
    return string

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
