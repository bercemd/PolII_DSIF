#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# B Dutagaci, B Bogart, July, 2021


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
    '--sela',type=str, default='protein',
    help='First selection.')
parser.add_argument(
    '--selb',type=str, default='protein',
    help='Second selection.')
parser.add_argument(
    '--out_path',type=str, default='contact.dat',
    help='Output path.')
arg = parser.parse_args()

def contacts_mda(u,selecta,selectb):
        sel = u.select_atoms(selecta, updating=True)
        ref = u.select_atoms(selectb, updating=True)

        #getting residue names and residue number from selected groups
        sel_res = sel.residues
        sel_seg = sel_res.segids
        sel_resids = sel_res.resids
        sel_nums = sel_res.resnums

        ref_res = ref.residues
        ref_seg = ref_res.segids
        ref_resids = ref_res.resids
        ref_nums = ref_res.resnums

        #empty array to stack all frame's counts of contacts at positions i,j
        count_stack = np.zeros((len(sel_res), len(ref_res)))

        #iteration over the trajectory
        for frame in u.trajectory:
        #for frame in u.trajectory[1:10]:
                print(frame)

                #array of 100's to be replaced at positions i,j by the minimum distance between residue i and residue j
                res_dist = np.full((len(sel_res), len(ref_res)), 100)

                #iterating through each of the contacting residues
                for n in sel_res:
                        for m in ref_res:
                                #selecting the atoms of sel and ref residues
                                sel_atoms = n.atoms.positions
                                ref_atoms = m.atoms.positions

                                #calculating the minimum distance between the residues
                                distance = distance_array(sel_atoms, ref_atoms, backend='OpenMP')
                                min_dist = np.min(distance)

                                #replacing the minimum distance value of residue i and residue j at position [i,j] in res_dist
                                sel_i = (np.where(n == sel_res)[0][0])
                                ref_j = (np.where(m == ref_res)[0][0])
                                res_dist[sel_i, ref_j] = min_dist

                #Calculating where in res_dist is less than or equal to a cutoff distance
                cutoff=5
                contacts = contact_matrix(res_dist, cutoff)
                #Converting True to 1, and False to 0
                contact_int = contacts * 1
                #adding the contacts in the frame to the total contacts in count_stack
                count_stack = count_stack + contact_int

        #Averaging the total number of contacts by the number of frames in the trajectory (needs to be transposed after division)
        count_avg = count_stack / len(u.trajectory)
        #counts = count_avg.transpose()
        counts = count_avg

        return counts, sel_resids, ref_resids, sel_seg, ref_seg

dcds = arg.traj_path
U = mda.Universe(arg.psf_path, dcds, in_memory=True)
contacts_arr, sela_resids, selb_resids, sela_seg, selb_seg = contacts_mda(U,arg.sela,arg.selb)
contacts_arr.shape
outputfile = open(arg.out_path,"w")
for i in range(len(contacts_arr)):
  for j in range(len(contacts_arr[i])):
    print(sela_resids[i],selb_resids[j],contacts_arr[i][j],sela_seg[i],selb_seg[j],file=outputfile)
outputfile.close()

