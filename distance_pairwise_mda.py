#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# B Dutagaci, 2022

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PDB_small
from MDAnalysis.analysis import distances
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
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
    '--sela',type=str, default='protein and resid 1',
    help='First selection.')
parser.add_argument(
    '--selb',type=str, default='protein and resid 1',
    help='Second selection.')
parser.add_argument(
    '--out_path',type=str, default='distance.dat',
    help='Output path.')
arg = parser.parse_args()

def dist_mda(u,selecta,selectb):
  sela = u.select_atoms(selecta)
  selb = u.select_atoms(selectb)

  n_sela = len(sela)
  n_selb = len(selb)

  print('Sela has {} residues and Selb has {} residues'.format(n_sela, n_selb))

  sela_res = sela.residues
  selb_res = selb.residues
  sela_resids = sela_res.resids
  selb_resids = selb_res.resids
  dist_arr = np.zeros(shape=(len(u.trajectory)))
  boxsize = u.trajectory.ts.dimensions
  i=0
  for frame in u.trajectory:
    for n in range(len(sela_res)):
      for m in range(len(selb_res)):
        distance = distances.distance_array(sela_res[n].atoms.positions,selb_res[m].atoms.positions,box=boxsize,backend='OpenMP')
        min_dist = np.min(distance)
        dist_arr[i] = min_dist
    i=i+1
  return dist_arr, sela_resids, selb_resids
outputfile = open(arg.out_path,"w")
U = mda.Universe(arg.psf_path, arg.traj_path, in_memory=True)

dist_arr, sela_resids, selb_resids = dist_mda(U,arg.sela,arg.selb)
dist_arr.shape
for i in range(len(dist_arr)):
    print(dist_arr[i],file=outputfile)
outputfile.close()

