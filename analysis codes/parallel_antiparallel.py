# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:26:49 2020

make arrays instead of files
 
@author: Sudeep

purpose: set of functions/definitions to help analyze results from PRIME20 simulations
    this version is written to be ran on a pc or mac, vs a linux system

bash script to rewrite all pdb files
and also zip them
should zip together pdb and bptnr files
zip catch44results.zip *pdb* *bptnr*

"""
from fortran import FortranFile
import pandas as pd
import numpy as np
#import networkx as nx
import math
import sys, math, array, csv, os, time, glob

def read_bptnr(filename,nchain,noptotal,nb):
    '''
    Parameters
    ----------
    filename : bptnr file
        two rows, first row is bead index, second is the hydrogen bonding partner
    nchain : int
        total number of chains in system
    noptotal : int
        total number of beads in system
    nb: number of beads in a chain

    Returns
    -------
    collision_num : int
        ncoll
    hb_data : 2 by noptotal matrix
        readable version of bptnr file
    hbmatrix : TYPE
        nchain by nchain matrix of hydrogen bonds
    num_of_hb : int
        number of hydrogen bonds

    '''
    # nb - number of beads per chain
    # noptotal - number of total beads in the system
    dtype=[('col','i8'),('hb','i4',(noptotal))]
    collision_num=0 #arbitrary just need to set
    # rel=f.read_ints() #reads all integers?
    # forces it to read through a whole file to make 
    # sure we get the last set of data points    
    with FortranFile(filename,'r') as f:
      while (( collision_num%1e7) != 1.0):
        record=f.read_record(dtype)
        collision_num=record['col'][0]
    hb_data=np.zeros([noptotal,2], dtype=int)
    for j in range(noptotal):
        hb_data[j,:]=[j+1,record['hb'][0][j]]
        
    ''' hydrogen bonding matrix'''
    hbmatrix=np.zeros((nchain,nchain),dtype=int)  
    
    hb_data_chain=hb_data[hb_data[:,1]!=0] # removing lack of hbond beads
    hb_data_chain=np.asarray( (hb_data_chain-1)/nb, dtype=int) # chain indexing
    for ii in range(len(hb_data_chain)):
        hbmatrix[hb_data_chain[ii][0]][hb_data_chain[ii][1]] = hbmatrix[hb_data_chain[ii][0]][hb_data_chain[ii][1]] + 1

    #num_of_hb=sum(hb_data[:,1]!=0)/2 # this works as well, just here to double check
    num_of_hb=sum(sum(hbmatrix))/2
    return collision_num, hb_data, hbmatrix, num_of_hb


def anti_para(bptnr_file,nchain,nb1,nres):
  '''
    returns the number of antiparallel and parallel peptide pairs
    function works for two component systems that have the same number of chains and residues
    
    Parameters
    ----------
    hb_data : TYPE
        is an array of hydrogen bonded pairs
    hbmatrix : TYPE
        DESCRIPTION.
    nchain : TYPE
        number of total chains in the simulation
    nb1 : TYPE
        number of beads in a chain (remember to subtract glycines)
    nres : TYPE, optional
        number of residues per chain. The default is 11.

    Returns
    -------
    anti : integer
        number of antiparallel peptides
    para : integer
        number of parallel peptides
    order : array, matrix?
        gives the neigboring chains (i and j), parallel vs antiparallel, 
        last digit gives register shift???

  ''' 
  noptotal = nchain*nb1
  collision_num, hb_data, hbmatrix, num_of_hb = read_bptnr(bptnr_file,nchain,noptotal,nb1)
  
  anti=0; para=0; order=[]
  for i in range(nchain-1):  
    for j in range(i+1,nchain):
      if ( hbmatrix[i][j]>2):
        #i and j pair are in a beta-sheet, they have at least 3 hydrogen bonds
        nh_i=np.zeros(nres)
        nh_j=np.zeros(nres)
        co_i=np.zeros(nres)
        co_j=np.zeros(nres)
        for ii in range(nres):
          #jj is going through NH beads
          jj=hb_data[i*nb1+nres+ii,1]
          if(int((jj-1)/nb1)==j): # then j is the partner chain in the loop
            nh_i[ii]=i*nb1+nres+ii+1
            co_j[ii]=jj
        nh=nh_i[nh_i!=0]%nres # gives bead position in chain
        co=co_j[co_j!=0]%nres
        nh[nh==0]=nres # need this because python indexes at 0 
        co[co==0]=nres
        #determining if values are increasing or decreasing
        a=sum(np.diff(co_j[co_j!=0]))
        
        # to be more rigorous
        for jj in range(nres):
          ii=hb_data[j*nb1+nres+jj,1]
          if(int((ii-1)/nb1)==i and ii!=0):
            nh_j[jj]=j*nb1+nres+jj+1
            co_i[jj]=ii
        nh2=nh_j[nh_j!=0]%nres
        co2=co_i[co_i!=0]%nres
        nh2[nh2==0]=nres
        co2[co2==0]=nres
        b=sum(np.diff(co_i[co_i!=0]))
        
        if ((a>0 and b>=0) or (a>=0 and b>0)) :
          para=para+1
          #order.append([i,j,1]) #,abs(nh[0]-co[0])])
        elif ((a<0 and b<=0) or (a<=0 and b<0)):
          anti=anti+1
          #order.append([i,j,-1]) #,abs(12-nh[0]-co[0])])
        else:
          print("error in determining direction in pair:", i, j, a, b)  
          print("this can happen when a peptide has a turn")
          print(nh, co, nh2, co2)
  #order=np.asarray(order)
  return anti, para
  
def main():
    
    #nchain = number of peptide chains in the DMD/PRIME20 simulations
    #nres = number of residues on the peptide
    #nb1 = number of beads in a peptide chain
    #.bptnr files are output files are produced in DMD/PRIME20 simulations. These files are available on request
    
    nchain = 48
    nb1 = 27
    nres = 7
    
    os.chdir(path)
    print(os.getcwd())
    bptnr_files = glob.glob('*.bptnr')
    bptnr_files.sort()
    if 'run0000.bptnr' in bptnr_files:
        bptnr_files.remove('run0000.bptnr')
    
    #numFrames = len(bptnr_files)
    numFrames = 508
    #numFilesToRead = int(round(0.05*numFrames))
    numFilesToRead = 25
    #bptnr_files_analyze =  bptnr_files[-numFilesToRead:]
    bptnr_files_analyze =  bptnr_files[458:508]
     
    antiparallel, parallel = [],[]
    for i in range(len(bptnr_files_analyze)):
        anti, para = anti_para(bptnr_files_analyze[i],nchain,nb1,nres)
        antiparallel.append(anti)
        parallel.append(para)
    
    para_percentage = [x / nchain for x in parallel]    
    parallel_mean = float(sum(para_percentage)/len(para_percentage))
    print(parallel_mean)
    parallel_std = np.std(para_percentage)
    print(parallel_std)
    
main()    
  