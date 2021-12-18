# Programed in Python3 but run in python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 15:39:49 2020

@author: carlos
"""


import sys
import subprocess as sp
import os

if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]
                
    print("\nCSS - A CFD code:\n")    

    os.system("rm ./executable %ssolution.csv" % (path))
    os.system("gcc ./readTables.c ./input.c ./mesh.c ./solver.c -o ./executable -lm -fopenmp")
    os.system("./executable %s" % path)
    os.system("python3 %sanalisys.py %s" % (path, path))
