#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# B Dutagaci, B Bogart, July, 2021


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
    '--sela',type=str, default='protein',
    help='First selection.')
parser.add_argument(
    '--selb',type=str, default='protein',
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
  dist_arr = np.zeros(shape=(len(sela_res), len(selb_res)))
  boxsize = u.trajectory.ts.dimensions
  for frame in u.trajectory:
    for n in range(len(sela_res)):
      for m in range(len(selb_res)):
        distance = distances.distance_array(sela_res[n].atoms.positions,selb_res[m].atoms.positions,box=boxsize,backend='OpenMP')
        min_dist = np.min(distance)
        dist_arr[n,m] += min_dist
  dist_arr = np.divide(dist_arr,len(u.trajectory))
  return dist_arr, sela_resids, selb_resids

U = mda.Universe(arg.psf_path, arg.traj_path, in_memory=True)
dist_arr, sela_resids, selb_resids = dist_mda(U,arg.sela,arg.selb)
dist_arr.shape
outputfile = open(arg.out_path,"w")
for i in range(len(dist_arr)):
  for j in range(len(dist_arr[i])):
    print(sela_resids[i],selb_resids[j],dist_arr[i][j],file=outputfile)
outputfile.close()

