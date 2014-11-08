"""
Python Wrapper for the SHTOOLS library by Mark Wieczorek

Authors: Matthias Meschede, Mark Wieczorek, 2014
"""

def load_documentation():
    import _SHTOOLS
    import os
    import re
    print 'loading shtools documentation'
    redescr = re.compile('DESCRIPTION(.*?)=head1 ARGUMENTS',re.DOTALL)
    renotes = re.compile('NOTES(.*?)=head1',re.DOTALL)
    DOC_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'doc'))
    for name,func in _SHTOOLS.__dict__.items():
        if callable(func):
            try:
                path = os.path.join(DOC_DIR,name.lower()+'.pod')
                docfile = open(path,'r').read()
                description = '\nDESCRIPTION\n-----------' + redescr.search(docfile).group(1) + '\n'
                try:
                    notes       = 'NOTES\n-----'         + renotes.search(docfile).group(1) + '\n'
                except:
                    notes = ''
                func.__doc__ = 'SIGNATURE\n--------\n\n'+func.__doc__ +description + notes
            except Exception,msg:
                print msg

#try to load documentation if interactive
#import __main__ as main
#if not hasattr(main, '__file__'):
#    load_documentation()

#import into main namespace
from _constants import planetsconstants as constants

from _SHTOOLS import *
