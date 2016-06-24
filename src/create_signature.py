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

#==== MAIN FUNCTION ====


def main():
    fname_wrapper = 'PythonWrapper.f95'
    fname_signature = 'pyshtools.pyf'

    print('now cracking Fortran file SHTOOLS.f95 using f2py function...')
    crackfortran.verbose = False
    crackfortran.dolowercase = False
    cracked_shtools = crackfortran.crackfortran(fname_wrapper)
    for subroutine in cracked_shtools:
        subroutine['f2pyenhancements'] = {'fortranname': subroutine['name'].lower()}
    for subroutine in cracked_shtools:
        subroutine['name'] = subroutine['name'][2:]
    interface = {'block': 'interface', 'name': 'unknown_interface', 'from': '',
                 'body': cracked_shtools, 'externals': [], 'interfaced': [],
                 'vars': {}}
    module = {'block': 'python module', 'name': '_SHTOOLS', 'from': '',
              'body': interface, 'externals': [], 'interfaced': [],
              'vars': {}}
    out = crackfortran.crack2fortran(module)
    with open(fname_signature, 'w') as outfile:
        outfile.write(out)

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
