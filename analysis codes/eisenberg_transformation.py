#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 9/21/2021
@author: Sudeep
"""
#Code to perform crystallographic transformation from a single peptide

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import sys

#pdb_in = PDB file input
#pdb_out = PDB file output
#num_sheets = number of sheets in the peptide assembly
#pep_per_sheet = number of peptides in sheet

pdb_in='comp.pdb'
pdb_out='comp_o.pdb'
num_sheets=2
pep_per_sheet=24.0

with open(pdb_in,'r') as f2:
    lines=f2.readlines()
    z=len(lines)
    anum=np.zeros(z,dtype=int)
    atype=["" for aa in range(z)]
    rtype=["" for aa in range(z)]
    rnum=np.zeros(z,dtype=int)
    xcoo=np.zeros(z,dtype=float)
    ycoo=np.zeros(z,dtype=float)
    zcoo=np.zeros(z,dtype=float)
    for a, line in enumerate(lines):
        fields=line.strip().split()
        anum[a]=fields[1]
        atype[a]=str(fields[2])
        rtype[a]=str(fields[3])
        rnum[a]=fields[4]
        xcoo[a]=fields[5]
        ycoo[a]=fields[6]
        zcoo[a]=fields[7]
f2.close
        
#Perform Matrix Operations
#Matrix operation information to create peptide-assembly is generally available in the PDB file from the Protein Data Bank

peplength=rnum[z-1]    
r1=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
r2=np.array([[-1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0]])
t1=np.array([[0.000],[0.000],[0.000]])
t2=np.array([[25.000],[-2.467],[18.770]])
t3=np.array([[0.000],[4.934],[0.000]])

#Output new PDB file

with open(pdb_out,'a') as f2:
    for c in range(2):
        if(c==0):
            r=r1
            t=t1
        elif(c==1):
            r=r2
            t=t2
        for b in range(int(pep_per_sheet)):
            t4=t+int(b)*t3
            print(t4)
            for a in range(z):
                v=np.array([[xcoo[a]],[ycoo[a]],[zcoo[a]]])
                x=np.dot(r,v)
                x=x+t4
                x1_str=np.array2string(x[0],precision=3,formatter={'float_kind':'{0:.3f}'.format})
                x2_str=np.array2string(x[1],precision=3,formatter={'float_kind':'{0:.3f}'.format})
                x3_str=np.array2string(x[2],precision=3,formatter={'float_kind':'{0:.3f}'.format})
                y1=x1_str[1:-1]
                y2=x2_str[1:-1]
                y3=x3_str[1:-1]
                anum_pdb=anum[a]
                rnum_pdb=rnum[a]
                if(len(str(anum_pdb))==1):
                    space1="      "
                elif(len(str(anum_pdb))==2):
                    space1="     "
                elif(len(str(anum_pdb))==3):
                    space1="    "
                elif(len(str(anum_pdb))==4):
                    space1="   "
                if(len(str(atype[a]))==4):
                    space2=" "
                else:
                    space2="  "
                if(len(str(atype[a]))==1):
                    space3="   "
                elif(len(str(atype[a]))==2):
                    space3="  "
                elif(len(str(atype[a]))==3 or len(str(atype[a]))==4):
                    space3=" "
                if(len(str(rnum_pdb))==1):
                    space4="     "
                elif(len(str(rnum_pdb))==2):
                    space4="    "
                elif(len(str(rnum_pdb))==3):
                    space4="   "
                if(len(str(y1))==5):
                    space5="       "
                elif(len(str(y1))==6):
                    space5="      "
                elif(len(str(y1))==7):
                    space5="     "    
                if(len(str(y2))==5):
                    space6="   "
                elif(len(str(y2))==6):
                    space6="  "
                elif(len(str(y2))==7):
                    space6=" "    
                if(len(str(y3))==5):
                    space7="   "
                elif(len(str(y3))==6):
                    space7="  "
                elif(len(str(y3))==7):
                    space7=" "    
                L1=["ATOM", space1, str(anum_pdb), space2, str(atype[a]), space3, str(rtype[a]), space4, \
                    str(rnum_pdb), space5, y1, space6, y2, space7, y3, "  ", "1.00", "  ", "0.00", "\n"]
                f2.writelines(L1)
            f2.writelines(["TER\n"])
            anum=anum+z
            rnum=rnum+peplength
        
f2.close



          