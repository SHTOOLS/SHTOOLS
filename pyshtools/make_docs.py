#!/usr/bin/env python
"""
This script builds the python documentation from the function signature and the
".pod" files that document the Fortran library. The processed documentation is
saved as ascii text files which are loaded on runtime and replace the __doc__
string of the f2py wrapped functions.
"""

import sys, os, re
import _SHTOOLS

def main():
    #---- input/output folders ----
    libfolder = os.path.abspath(sys.argv[1])
    poddocfolder = os.path.join(libfolder,'src/doc')
    pydocfolder = os.path.join(libfolder,'pyshtools/doc')
    print '---- searching documentation in folder: {} ----'.format(poddocfolder)
    
    #---- pod file search patterns ----
    redescr = re.compile('DESCRIPTION(.*?)=head1 ARGUMENTS',re.DOTALL)
    renotes = re.compile('NOTES(.*?)=head1',re.DOTALL)

    #---- loop through the f2py _SHTOOLS functions and make docstrings ----
    for name,func in _SHTOOLS.__dict__.items():
        if callable(func):
            try:
                #-- process and load documentation from different sources --
                #get and process f2py generated call signature:
                callsign = process_f2pydoc(func.__doc__)
                signature = 'SIGNATURE\n---------\n' + callsign
    
                #read pod file documentation:
                fname_poddoc = os.path.join(poddocfolder,name.lower()+'.pod')
                podfile   = open(fname_poddoc,'r')
                podstring = podfile.read()
    
                match = redescr.search(podstring)
                if match is None:
                    description = '\nDESCRIPTION\n-----------\nDOCUMENTATION NOT AVAILABLE'
                else:
                    description = '\nDESCRIPTION\n-----------' + match.group(1)
    
                match = renotes.search(podstring)
                if match is None:
                    notes       = ''
                else:
                    notes       = '\nNOTES\n-----'+ match.group(1)
    
                podfile.close()
    
                #combine docstring:
                docstring       =  signature + description + notes
    
                #-- save combined docstring in the pydoc folder--
                fname_pydoc  = os.path.join( pydocfolder,name.lower()+'.doc')
                pydocfile = open(fname_pydoc,'w')
                pydocfile.write(docstring)
                pydocfile.close()
    
            except IOError,msg:
                print msg

#===== HELPER FUNCTIONS ====
def process_f2pydoc(f2pydoc):
    """
    this function replace all optional _d0 arguments with their default values
    in the function signature. These arguments are not intended to be used and
    signify merely the array dimensions of the associated argument.
    """
    #I.  split f2py document in its parts
    docparts    = re.split('\n--',f2pydoc)

    #II. replace arguments with _d suffix with empty string in function signature (remove them):
    docparts[0] = re.sub('[\[(,]\w+_d\d','',docparts[0])

    # combine doc parts to a single string
    processed_signature = '\n--'.join(docparts)

    return processed_signature

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
