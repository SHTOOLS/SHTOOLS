#!/usr/bin/env python3
"""
Automatically creates python-wrapper subroutines from the interface file
SHTOOLS.f95. Unfortunately all assumed array shapes have to be changed because
their structure is only known by the Fortran compiler and can not be directly
exposed to C. It is possible that newer f2py versions can handle assumed array
shapes using a similar procedure.
"""
from numpy.f2py import crackfortran
from copy import deepcopy
import re

# ==== MAIN FUNCTION ====
def main():
    fname_fortran = 'SHTOOLS.f95'
    fname_wrapper = 'cWrapper.f95'

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

   
    # search for the indices of 'use shtools,' to insert 'implicit none' after
    iusestatement = [iline for iline, line in enumerate(wrapperlines)
                     if 'use shtools,' in line]
    assert len(iusestatement) == len(interface_new['body']), \
        'number of subroutines don\'t match'
        
    # search for the indices of 'end subroutine'
    iendsubroutine = [iline for iline, line in enumerate(wrapperlines)
                      if 'end subroutine' in line or 'end function' in line]
    assert len(iendsubroutine) == len(interface_new['body']), \
        'number of subroutines don\'t match'
    
    print('sort variables...')
    for i in range(len(iusestatement)):
        declaration = wrapperlines[iusestatement[i]+1:iendsubroutine[i]]
        current_order =  [ (line.split('::')[1]).strip() for line in declaration ]
        order = interface_new['body'][i]['sortvars']
        
        idx = []
        for j in range(len(order)):
            idx.append( current_order.index(order[j]) )
            
        declaration[:] = [declaration[k] for k in idx]
        wrapperlines[iusestatement[i]+1:iendsubroutine[i]] =  declaration[:] 
       
    print('add implicit none statements')
    for iline in iusestatement[::-1]:
        wrapperlines.insert(iline + 1,
                            2 * crackfortran.tabchar + 'implicit none')

    print('add shtools subroutine calls...')
    
     # search for the indices of 'end subroutine'
    iendsubroutine = [iline for iline, line in enumerate(wrapperlines)
                      if 'end subroutine' in line or 'end function' in line]
    assert len(iendsubroutine) == len(interface_new['body']), \
        'number of subroutines don\'t match'
    
    # insert call statements before 'end subroutine' line starting from the
    # end such that we don't change the preceding indices
    for sroutine_new, sroutine_old, iline in list(zip(interface_new['body'],
                                                 interface_old['body'],
                                                 iendsubroutine))[::-1]:
        args = create_arg_list(sroutine_old)
        if sroutine_new['block'] == 'function':
            newline = 2 * crackfortran.tabchar +\
                '%s=%s(' % (sroutine_new['name'], sroutine_old['name']) +\
                args+ ')'
        elif sroutine_new['block'] == 'subroutine':
            newline = 2 * crackfortran.tabchar +\
                'call %s(' % sroutine_old['name'] +\
                args+ ')'
        wrapperlines.insert(iline + 1, '')
        wrapperlines.insert(iline, newline)
        
    print('add bind statment...')
    p = re.compile('\s*(subroutine|function)')
    # search for the indices of 'subroutine'
    isubroutine = [iline for iline, line in enumerate(wrapperlines)
                      if p.match(line) is not None]
    assert len(isubroutine) == len(interface_new['body']), \
        'number of subroutines don\'t match'
    
    for sroutine_new, sroutine_old, iline in list(zip(interface_new['body'],
                                                interface_old['body'],
                                                isubroutine))[::-1]:
       wrapperlines[iline] = wrapperlines[iline] + ' bind(c, name=\"' \
         + sroutine_new['name'] + '\")'
       newline = 2 * crackfortran.tabchar + 'use, intrinsic :: iso_c_binding'
       wrapperlines.insert(iline+1, newline)
              
    print('writing wrapper to file %s' % fname_wrapper)
    for iline, line in enumerate(wrapperlines):
        try:
            firstword = line.split()[0]
            secondword = line.split()[1]
            words = ['real(kind=c_double)', 'complex(kind=c_double_complex)', 
                     'integer(kind=c_int)', 'character(kind=c_char)',
                     'real(kind=8)','real*8', 'integer', 'integer(kind=4)', 'character*(*)',
                     'complex*16']
            for word in words:
                if (firstword == word and not secondword[0] == ':' or
                        secondword[0] == ','):
                    line = line.replace(word, word + ',')
            wrapperlines[iline] = line
        except IndexError:
            pass

    with open(fname_wrapper, 'w') as outfile:
        for line in wrapperlines[4:-5]:
            line = line.replace('! in SHTOOLS.f95:SHTOOLS:unknown_interface',
                                '')
            if len(line) <= 80:
                outfile.write(line + '\n')
            else:
                elems = line.split(',')
                newline = elems[0]
                for elem in elems[1:]:
                    if len(newline) > 80:
                        outfile.write(newline + '&\n')
                        newline = ' ' * len(elems[0])
                    newline += ',' + elem
                outfile.write(newline + '\n')

    print('\n==== ALL DONE ====\n')


# ==== FUNCTIONS ====


def create_arg_list(subroutine):
   
    sig = []
    for arg in subroutine['args']:
         if 'attrspec' in subroutine['vars'][arg] and 'optional' in subroutine['vars'][arg]['attrspec']:
                 sig.append(arg + '=' + arg)
         else: 
             sig.append(arg)
             
    return ','.join(sig) 
                 


def modify_subroutine(subroutine):
    """loops through variables of a subroutine and modifies them"""
    # print('\n----',subroutine['name'],'----')

    # use original function from shtools:
        
    prepend = 'c'
        
    subroutine['use'] = {'shtools':
                         {'map': {subroutine['name']: subroutine['name']},
                          'only': 1}}

    # -- loop through variables:
    for varname, varattribs in list(subroutine['vars'].items()):
        
        # prefix function return variables with prepend
        if varname == subroutine['name']:                
            subroutine['vars'][prepend + varname] = \
                subroutine['vars'].pop(varname)
            subroutine['sortvars'][subroutine['sortvars'].index(varname)] = \
                prepend + varname
            varname = prepend + varname
            # print('prefix added:',varname)
        # -- change assumed to explicit:
        if has_assumed_shape(varattribs):
            make_explicit(subroutine, varname, varattribs)

        if 'dimension' not in varattribs:
            is_intent_in = 'intent' in varattribs and 'in' in varattribs['intent']
            is_optional = 'attrspec' in varattribs and 'optional' in varattribs['attrspec']
            if is_intent_in and not is_optional:
                varattribs['attrspec'].append('value')

         
    # change type
    for varname, varattribs in list(subroutine['vars'].items()):
        if 'kindselector' not in varattribs:
            varattribs['kindselector'] = {}
            
        # modify kind
        if varattribs['typespec'] == 'real':
            varattribs['kindselector']['kind']='c_double'
        elif varattribs['typespec'] == 'integer':
            varattribs['kindselector']['kind']='c_int'
        elif varattribs['typespec'] == 'complex':
            varattribs['kindselector']['kind']='c_double_complex'
        elif varattribs['typespec'] == 'character':
            varattribs['kindselector']['kind']='c_char'



    # add py prefix to subroutine:
    subroutine['name'] = prepend + subroutine['name']      
        

def insert_dim(subroutine, dimname, arg_pos, declartaion_pos):
    dimattribs = {'attrspec': ['value'], 'typespec': 'integer', 'intent': ['in']}
    # declare dimension in subroutine variables
    subroutine['vars'][dimname] = dimattribs
    # add dimension to subroutine arguments
    subroutine['args'].insert(arg_pos, dimname)
    # place indices at the beginning
    subroutine['sortvars'].insert(declartaion_pos, dimname)

def make_explicit(subroutine, varname, varattribs):
    
    argpos = subroutine['args'].index(varname) + 1 
    decpos = subroutine['sortvars'].index(varname) 
    
    if varattribs['typespec'] == 'character':
        dimname = '%s_d%d' % (varname, 1)
        varattribs['dimension'] = [dimname]
        insert_dim(subroutine,dimname,argpos,decpos)
    elif varname=='cilm':
        dimname = 'cilm_d'
        varattribs['dimension'] = ['2', dimname, dimname]
        insert_dim(subroutine,dimname,argpos,decpos)
    else:
        for idim, size in enumerate(varattribs['dimension']):
            if size == ':':
                # change assumed array to explicit
                dimname = '%s_d%d' % (varname, idim)
                varattribs['dimension'][idim] = dimname
                insert_dim(subroutine,dimname,argpos,decpos)
                decpos = decpos+1
                argpos = argpos + 1


def has_assumed_shape(varattribs):
    """checks if variable has assumed shape"""
    if varattribs['typespec'] == 'character' and 'charselector' in varattribs:
        if '*' in varattribs['charselector']['*']:
            return True 
    try:
        if ':' in varattribs['dimension']:
            return True
        else:
            return False
    except KeyError:
        return False


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
