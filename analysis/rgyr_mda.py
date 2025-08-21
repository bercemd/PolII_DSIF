#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# B Dutagaci, Weththasinghage D. Amith, 2023 - Modified from the Original script from MDAnalysis

import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms
import MDAnalysis.analysis.align as align
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
    '--selection',type=str, default='protein and (name CA or name CB)',
    help='Selection for alignment.')
parser.add_argument(
    '--out_path',type=str, default='rgyr.dat',
    help='Output path.')
arg = parser.parse_args()

u = mda.Universe(arg.psf_path, arg.traj_path)

outputfile=open(arg.out_path,"w")

for ts in u.trajectory:     # iterate through all frames
    sele = u.select_atoms(arg.selection)
    rgyr = sele.radius_of_gyration()  # method of AtomGroup
    print("{0} {1}".format(
          ts.frame, rgyr),file=outputfile)

outputfile.close()

