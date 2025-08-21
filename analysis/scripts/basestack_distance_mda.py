#!/usr/bin/env python

# B Dutagaci, 2022

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.contacts import contact_matrix
from MDAnalysis.analysis.distances import distance_array

import os,sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '--traj_path',type=str, default='md.dcd',
    help='Trajectory path.')
parser.add_argument(
    '--psf_path',type=str, default='protein.psf',
    help='PSF path.')
parser.add_argument(
    '--out_path',type=str, default='basestacking_dist.dat',
    help='Output path.')
arg = parser.parse_args()

def distance_mda(u,filename):
  select = "((segid N01N and resid 1:12) or (segid N01T and resid 32:43))"
  sela = u.select_atoms(select, updating=True)
  selb = u.select_atoms(select, updating=True)
  sela_res = sela.residues
  selb_res = selb.residues

  i=0
  #iteration over the trajectory
  for frame in u.trajectory:
    print(frame)
    distance = ""
    for n in sela_res:
      for m in selb_res:
        if (n.resid-m.resid==1) and (m.resid<12):
           xx = u.select_atoms("((segid N01N) and (resid %s) and (name N1 or name C2 or name N3 or name C4 or name C5 or name C6 or name N7 or name C8 or name N9))"%n.resid)
           yy = u.select_atoms("((segid N01N) and (resid %s) and (name N1 or name C2 or name N3 or name C4 or name C5 or name C6 or name N7 or name C8 or name N9))"%m.resid)
           xx_com = xx.center_of_mass()
           yy_com = yy.center_of_mass()
           distance = distance+str(np.linalg.norm(xx_com-yy_com))+" "
        elif (n.resid-m.resid==1) and (m.resid>=32):
           xx = u.select_atoms("((segid N01T) and (resid %s) and (name N1 or name C2 or name N3 or name C4 or name C5 or name C6 or name N7 or name C8 or name N9))"%n.resid)
           yy = u.select_atoms("((segid N01T) and (resid %s) and (name N1 or name C2 or name N3 or name C4 or name C5 or name C6 or name N7 or name C8 or name N9))"%m.resid)
           xx_com = xx.center_of_mass()
           yy_com = yy.center_of_mass()
           distance = distance+str(np.linalg.norm(xx_com-yy_com))+" "
    print(distance,file=filename)
    i=i+1

dcds = arg.traj_path 

U = mda.Universe(arg.psf_path, dcds, in_memory=True)
outputfile = open(arg.out_path,"w")
distance_mda(U,outputfile)
outputfile.close()

