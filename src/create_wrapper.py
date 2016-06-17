#!/usr/bin/env python
"""
Automatically creates python-wrapper subroutines from the interface file
SHTOOLS.f95. Unfortunately all assumed array shapes have to be changed because
their structure is only known by the Fortran compiler and can not be directly
exposed to C. It is possible that newer f2py versions can handle assumed array
shapes using a similar procedure.
"""
from __future__ import absolute_import, division, print_function

#==== IMPORTS ====
from numpy.f2py import crackfortran
import re
from copy import deepcopy

#==== MAIN FUNCTION ====


def main():
    fname_fortran = 'SHTOOLS.f95'
    fname_wrapper = 'PythonWrapper.f95'

    print('now cracking Fortran file SHTOOLS.f95 using f2py function...')
    crackfortran.verbose = False
    crackfortran.dolowercase = False
    cracked_shtools = crackfortran.crackfortran(fname_fortran)

    print('decending through shtools module tree...')
    module = cracked_shtools[0]
    interface_old = module['body'][0]
    interface_new = deepcopy(interface_old)
    for subroutine in interface_new['body']:
        modify_subroutine(subroutine)

    print('create interface string...')
    wrapper = crackfortran.crack2fortran(interface_new)
    wrapperlines = wrapper.split('\n')

    print('add implicit none statements')
    # search for the indices of 'use shtools,' to insert 'implicit none' after
    iusestatement = [iline for iline, line in enumerate(wrapperlines) if 'use shtools,' in line]
    assert len(iusestatement) == len(interface_new['body']), 'number of subroutines don\'t match'
    for iline in iusestatement[::-1]:
        wrapperlines.insert(iline + 1, 2 * crackfortran.tabchar + 'implicit none')

    print('add shtools subroutine calls...')
    # search for the indices of 'end subroutine'
    iendsubroutine = [iline for iline, line in enumerate(wrapperlines)
                      if 'end subroutine' in line or 'end function' in line]
    assert len(iendsubroutine) == len(interface_new['body']), 'number of subroutines don\'t match'

    # insert call statements before 'end subroutine' line starting from the end such that we
    # don't change the preceding indices
    for sroutine_new, sroutine_old, iline in zip(interface_new['body'],
                                                 interface_old['body'],
                                                 iendsubroutine)[::-1]:
        if sroutine_new['block'] == 'function':
            newline = 2 * crackfortran.tabchar +\
                '%s=%s(' % (sroutine_new['name'], sroutine_old['name']) +\
                ','.join(sroutine_old['args']) + ')'
        elif sroutine_new['block'] == 'subroutine':
            newline = 2 * crackfortran.tabchar +\
                'call %s(' % sroutine_old['name'] +\
                ','.join(sroutine_old['args']) + ')'
        wrapperlines.insert(iline + 1, '')
        wrapperlines.insert(iline, newline)

    print('writing wrapper to file %s' % fname_wrapper)
    for iline, line in enumerate(wrapperlines):
        try:
            firstword = line.split()[0]
            secondword = line.split()[1]
            words = ['real*8', 'integer', 'integer(kind=4)', 'character*(*)', 'complex*16']
            for word in words:
                if firstword == word and not secondword[0] == ':' or secondword[0] == ',':
                    line = line.replace(word, word + ',')
            wrapperlines[iline] = line
        except IndexError:
            pass

    with open(fname_wrapper, 'w') as outfile:
        for line in wrapperlines[4:-5]:
            line = line.replace('! in SHTOOLS.f95:SHTOOLS:unknown_interface',
                                '')
            if len(line) <= 100:
                outfile.write(line + '\n')
            else:
                elems = line.split(',')
                newline = elems[0]
                for elem in elems[1:]:
                    if len(newline) > 100:
                        outfile.write(newline + '&\n')
                        newline = ' ' * len(elems[0])
                    newline += ',' + elem
                outfile.write(newline + '\n')

    print('\n==== ALL DONE ====\n')

#==== FUNCTIONS ====


def modify_subroutine(subroutine):
    """loops through variables of a subroutine and modifies them"""
    # print('\n----',subroutine['name'],'----')

    #-- use original function from shtools:
    subroutine['use'] = {'shtools': {'map': {subroutine['name']: subroutine['name']}, 'only': 1}}

    #-- loop through variables:
    for varname, varattribs in subroutine['vars'].items():
        #-- prefix function return variables with 'py'
        if varname == subroutine['name']:
            subroutine['vars']['py' + varname] = subroutine['vars'].pop(varname)
            varname = 'py' + varname
            # print('prefix added:',varname)
        #-- change assumed to explicit:
        if has_assumed_shape(varattribs):
            make_explicit(subroutine, varname, varattribs)
            # print('assumed shape variable modified to:',varname,varattribs['dimension'])

    #-- add py prefix to subroutine:
    subroutine['name'] = 'py' + subroutine['name']


def make_explicit(subroutine, varname, varattribs):
    dimattribs = {'attrspec': [], 'typespec': 'integer', 'intent': ['in']}
    for idim, size in enumerate(varattribs['dimension']):
        if size == ':':
            # change assumed array to explicit
            dimname = '%s_d%d' % (varname, idim)
            varattribs['dimension'][idim] = dimname
            # declare dimension in subroutine variables
            subroutine['vars'][dimname] = dimattribs
            # add dimension to subroutine arguments
            subroutine['args'].append(dimname)


def has_assumed_shape(varattribs):
    """checks if variable has assumed shape"""
    try:
        if ':' in varattribs['dimension']:
            return True
        else:
            return False
    except KeyError:
        return False

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
