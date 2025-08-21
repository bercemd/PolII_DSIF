#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# B Dutagaci, B Bogart, July, 2021 - Modified from the Original script from MDAnalysis

import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms
import MDAnalysis.analysis.align as align
import matplotlib.pyplot as plt
import os,sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '--traj_path',type=str, default='md.kow1.ff14sb.bsc1/',
    help='Trajectory path.')
parser.add_argument(
    '--psf_path',type=str, default='protein.psf',
    help='PSF path.')
parser.add_argument(
    '--pdb_path',type=str, default='protein.pdb',
    help='PDB path.')
parser.add_argument(
    '--fit_selection',type=str, default='protein and (name CA or name CB)',
    help='Selection for alignment.')
parser.add_argument(
    '--rmsd_selection',type=str, default='protein and name CA',
    help='Selection for rmsd calculation.')
parser.add_argument(
    '--out_path',type=str, default='rmsd.dat',
    help='Output path.')
arg = parser.parse_args()

U1 = mda.Universe(arg.psf_path, arg.traj_path)
U2 = mda.Universe(arg.psf_path, arg.pdb_path)

R = mda.analysis.rms.RMSD(U1,U2,select=arg.fit_selection,groupselections=[arg.rmsd_selection])
R.run()

data1 = R.rmsd.T

outputfile = open(arg.out_path,"w")
for i in range(len(data1[3])):
  print(i+1, data1[3][i], file=outputfile)
outputfile.close()
