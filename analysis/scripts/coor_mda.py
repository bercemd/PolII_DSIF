#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# B Dutagaci, 2025

import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms
from MDAnalysis.analysis.rms import RMSF
import MDAnalysis.analysis.align as align
from MDAnalysis.coordinates.memory import MemoryReader
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
    '--pdb_path',type=str, default='protein.pdb',
    help='PDB path.')
parser.add_argument(
    '--ref_selection',type=str, default='protein',
    help='Selection for reference.')
parser.add_argument(
    '--fit_selection',type=str, default='protein and (name CA or name CB)',
    help='Selection for alignment.')
parser.add_argument(
    '--out_path',type=str, default='coor.dat',
    help='Output path.')
arg = parser.parse_args()


def coor_mda(u_in_mem, pdb, selection=arg.ref_selection):
	#Calculate average structure by loading entire trajectory into memory
        protein = u_in_mem.select_atoms(selection)

	#fitting trajectory to initial frame for better average structure
        init_frame = mda.Universe(pdb)
        aligner = align.AlignTraj(u_in_mem, init_frame, select=arg.fit_selection, in_memory=True).run()

        coor = []
        for frame in u_in_mem.trajectory:
            backbone = u_in_mem.select_atoms(arg.fit_selection)
            coor.append(backbone.atoms.positions.tolist())

        return coor

U = mda.Universe(arg.pdb_path, arg.traj_path, in_memory=True)	

outputfile = open(arg.out_path,"w")
coor = coor_mda(U, arg.pdb_path, selection=arg.ref_selection)
for i in range(len(coor)):
    line = ""
    for j in range(len(coor[i])):
        line = "%s %s %s %s"%(line,coor[i][j][0],coor[i][j][1],coor[i][j][2])
    print(line,file=outputfile)


outputfile.close()
