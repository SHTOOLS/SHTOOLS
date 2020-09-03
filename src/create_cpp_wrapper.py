#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 09:17:37 2020
automatically creats 'cWrapper.h' which can be included by c or cpp programs
"""
from numpy.f2py import crackfortran
import re
import os

TYPEMAP = {'integer':'int',
           'real':'double',
           'character': 'char',
           'complex':'double _Complex'}

BYVALUE = {True:' ',
           False:'* '}

# ==== MAIN FUNCTION ====
def main():
    fname_wrapper = 'cWrapper.f95'
    fname_signature = 'cWrapper.h'

    print('now cracking Fortran file ' + fname_wrapper + ' using f2py function...')
    crackfortran.verbose = False
    crackfortran.dolowercase = False
    cracked_shtools = crackfortran.crackfortran(fname_wrapper)
    
  
    with open(fname_signature, 'w') as outfile:
       
        outfile.write('#pragma once\n')
        outfile.write('namespace shtools{\n')
        outfile.write('extern "C"\n{\n')
        for subroutine in cracked_shtools:
            return_type = 'void'
            if subroutine['block'] == 'function':
                
                return_type = TYPEMAP[subroutine['vars'][subroutine['name']]['typespec']]
                
            newline = return_type + ' ' + subroutine['name'] \
                + '( ' + create_signature(subroutine) + ')' + ';'
            outfile.write(newline+'\n')
        outfile.write('}\n}\n')
    
    os.system('clang-format -i -style=Mozilla ' + fname_signature)

    
def create_signature(subroutine):
    args = []
    for arg in subroutine['args']:
        var = subroutine['vars'][arg]
                
        is_called_by_value = 'attrspec' in var and 'value' in var['attrspec']
        carg = TYPEMAP[var['typespec']] + BYVALUE[is_called_by_value] + arg
        if 'intent' in var and 'in' in var['intent']:
            carg = 'const ' + carg
        args.append(carg)
    return ', '.join(args) 
        

def camel_to_snake(name):
  name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
  return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
