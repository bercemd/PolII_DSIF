#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# B Dutagaci, July, 2025 

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
    '--fit_selection',type=str, default='protein and name CA',
    help='Selection for fitting.')
parser.add_argument(
    '--sela',type=str, default='protein and name CA',
    help='First selection for cross-correlation.')
parser.add_argument(
    '--selb',type=str, default='protein and name CA',
    help='Second selection for cross-correlation.')
parser.add_argument(
    '--out_path',type=str, default='corr.dat',
    help='Output path.')
arg = parser.parse_args()


def cross_corr(u_in_mem, pdb):
    #Calculate average structure by loading entire trajectory into memory
    protein = u_in_mem.select_atoms(arg.fit_selection)
    #fitting trajectory to initial frame for better average structure
    init_frame = mda.Universe(pdb)
    prealigner = align.AlignTraj(u_in_mem, init_frame, select=arg.fit_selection, in_memory=True).run()
    #building average structure
    reference_coordinates = u_in_mem.trajectory.timeseries(asel=protein).mean(axis=1)
    reference = mda.Merge(protein).load_new(reference_coordinates[:, None, :], order="afc")
    #align whole trajectory to reference structure to minimize RMSD
    aligner = align.AlignTraj(u_in_mem, reference, select=arg.fit_selection, in_memory=True).run()
    #RMSD is stored as aligner.rmsd if needed for inspection

    mobila =  u_in_mem.select_atoms(arg.sela)
    resnumsa = mobila.resnums

    mobilb =  u_in_mem.select_atoms(arg.selb)
    resnumsb = mobilb.resnums

    sela_res = mobila.residues
    sela_seg = sela_res.segids
    sela_resids = sela_res.resids
    sela_nums = sela_res.resnums

    selb_res = mobilb.residues
    selb_seg = selb_res.segids
    selb_resids = selb_res.resids
    selb_nums = selb_res.resnums

    dotlista = []
    dotlistb = []
    dotlistab = []

    resa = []
    resb = []
    fi=0
    for f in u_in_mem.trajectory:
        backbone_a = u_in_mem.select_atoms(arg.sela)
        cooraf = backbone_a.atoms.positions

        ref_a = reference.select_atoms(arg.sela)
        refaf = ref_a.atoms.positions

        posit_a = []
        for i in range(len(cooraf)):
            posit_a.append(cooraf[i]-refaf[i])

        backbone_b = u_in_mem.select_atoms(arg.selb)
        coorbf = backbone_b.atoms.positions

        ref_b = reference.select_atoms(arg.selb)
        refbf = ref_b.atoms.positions

        posit_b = []
        for i in range(len(coorbf)):
            posit_b.append(coorbf[i]-refbf[i])

        dota = []
        dotb = []
        dotab = []
        for i in range(len(posit_a)):
            for j in range(len(posit_b)):
              dota.append(np.dot(posit_a[i],posit_a[i]))
              dotb.append(np.dot(posit_b[j],posit_b[j]))
              dotab.append(np.dot(posit_a[i],posit_b[j]))
              if fi==0:
                  resa.append(sela_resids[i])
                  resb.append(selb_resids[j])

        dotlista.append(dota)
        dotlistb.append(dotb)
        dotlistab.append(dotab)
        fi=fi+1      

    dotlista_array = np.asarray(dotlista)
    dotlistb_array = np.asarray(dotlistb)
    dotlistab_array = np.asarray(dotlistab)

    dotlista_array_t = dotlista_array.T
    dotlistb_array_t = dotlistb_array.T
    dotlistab_array_t = dotlistab_array.T
        
    return dotlista_array_t, dotlistb_array_t, dotlistab_array_t, resa, resb

U = mda.Universe(arg.pdb_path, arg.traj_path, in_memory=True)	

outputfile = open(arg.out_path,"w")
dotlista_array_t, dotlistb_array_t, dotlistab_array_t, resa, resb = cross_corr(U, arg.pdb_path)

for i in range(len(dotlista_array_t)):
    xa_avg = sum(dotlista_array_t[i])/len(dotlista_array_t[i])
    xb_avg = sum(dotlistb_array_t[i])/len(dotlistb_array_t[i])
    xab_avg = sum(dotlistab_array_t[i])/len(dotlistab_array_t[i])
    cross_corr = xab_avg/((xa_avg*xb_avg)**0.5)
    print(resa[i],resb[i],cross_corr,file=outputfile)

outputfile.close()
