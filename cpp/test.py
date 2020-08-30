#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 19:01:24 2020

@author: corbin
"""
import pyshtools as pysh
import numpy as np

if __name__ == '__main__':
    infile = '../examples/ExampleDataFiles/MarsTopo719.shape'
    cilm = pysh.shio.shread(infile,lmax=15)
    
    value = pysh.expand.MakeGridPoint(cilm[0], 10, 30)

