#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# B Dutagaci, 2022 - Modified from the Original script from MDAnalysis.analysis.nuclinfo module of MDAnalysis


import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PDB_small
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import nuclinfo
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
    '--angle',type=str, default='alpha',
    help='nucleic acid dihedral angle - select one of the followings: alpha, beta, gamma, delta, epsilon, zeta, chi')
parser.add_argument(
    '--out_path',type=str, default='dihedral.dat',
    help='Output path.')
arg = parser.parse_args()


def select_dihedral_a(u,seg,i):
    a = u.select_atoms(" atom {0!s} {1!s} O3\' ".format(seg, i - 1),
                              " atom {0!s} {1!s} P  ".format(seg, i),
                              " atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i))
    return a
def select_dihedral_b(u,seg,i):
    b = u.select_atoms(" atom {0!s} {1!s} P    ".format(seg, i),
                              " atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i))
    return b
def select_dihedral_g(u,seg,i):
    g = u.select_atoms(" atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i))
    return g
def select_dihedral_d(u,seg,i):
    d = u.select_atoms(" atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i))
    return d
def select_dihedral_e(u,seg,i):
    e = u.select_atoms(" atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i),
                              " atom {0!s} {1!s} P    ".format(seg, i + 1))
    return e
def select_dihedral_z(u,seg,i):
    z = u.select_atoms(" atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i),
                              " atom {0!s} {1!s} P    ".format(seg, i + 1),
                              " atom {0!s} {1!s} O5\' ".format(seg, i + 1))
    return z
def select_dihedral_c(u,seg,i):
    c = u.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                              " atom {0!s} {1!s} C1\' ".format(seg, i),
                              " atom {0!s} {1!s} N9 ".format(seg, i),
                              " atom {0!s} {1!s} C4  ".format(seg, i))
    if len(c) < 4:
        c = u.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                              " atom {0!s} {1!s} C1\' ".format(seg, i),
                              " atom {0!s} {1!s} N1 ".format(seg, i),
                              " atom {0!s} {1!s} C2  ".format(seg, i))
    return c

dcds = arg.traj_path 

U = mda.Universe(arg.psf_path, dcds, in_memory=True)

outputfile = open(arg.out_path,"w")

for frame in U.trajectory:
 if arg.angle=="alpha":
  alpha = ""
  seg = "N01N"
  ires=1
  fres=12
  for i in range(ires+1,fres):
    a = select_dihedral_a(U,seg,i)
    alpha = alpha+str(a.dihedral.value())+" "
  seg = "N01T"
  ires=32
  fres=43
  for i in range(ires+1,fres):
    a = select_dihedral_a(U,seg,i)
    alpha = alpha+str(a.dihedral.value())+" "
  print(alpha,file=outputfile)

 elif arg.angle=="beta":
  beta = ""
  seg = "N01N"
  ires=1
  fres=12
  for i in range(ires+1,fres):
    b = select_dihedral_b(U,seg,i)
    beta = beta+str(b.dihedral.value())+" "
  seg = "N01T"
  ires=32
  fres=43
  for i in range(ires+1,fres):
    b = select_dihedral_b(U,seg,i)
    beta = beta+str(b.dihedral.value())+" "
  print(beta,file=outputfile)

 elif arg.angle=="gamma":
  gamma = ""
  seg = "N01N"
  ires=1
  fres=12
  for i in range(ires+1,fres):
    g = select_dihedral_g(U,seg,i)
    gamma = gamma+str(g.dihedral.value())+" "
  seg = "N01T"
  ires=32
  fres=43
  for i in range(ires+1,fres):
    g = select_dihedral_g(U,seg,i)
    gamma = gamma+str(g.dihedral.value())+" "
  print(gamma,file=outputfile)

 elif arg.angle=="delta":
  delta = ""
  seg = "N01N"
  ires=1
  fres=12
  for i in range(ires+1,fres):
    d = select_dihedral_d(U,seg,i)
    delta = delta+str(d.dihedral.value())+" "
  seg = "N01T"
  ires=32
  fres=43
  for i in range(ires+1,fres):
    d = select_dihedral_d(U,seg,i)
    delta = delta+str(d.dihedral.value())+" "
  print(delta,file=outputfile)

 elif arg.angle=="epsilon":
  epsilon = ""
  seg = "N01N"
  ires=1
  fres=12
  for i in range(ires+1,fres):
    e = select_dihedral_e(U,seg,i)
    epsilon = epsilon+str(e.dihedral.value())+" "
  seg = "N01T"
  ires=32
  fres=43
  for i in range(ires+1,fres):
    e = select_dihedral_e(U,seg,i)
    epsilon = epsilon+str(e.dihedral.value())+" "
  print(epsilon,file=outputfile)

 elif arg.angle=="zeta":
  zeta = ""
  seg = "N01N"
  ires=1
  fres=12
  for i in range(ires+1,fres):
    z = select_dihedral_z(U,seg,i)
    zeta = zeta+str(z.dihedral.value())+" "
  seg = "N01T"
  ires=32
  fres=43
  for i in range(ires+1,fres):
    z = select_dihedral_z(U,seg,i)
    zeta = zeta+str(z.dihedral.value())+" "
  print(zeta,file=outputfile)

 elif arg.angle=="chi":
  chi = ""
  seg = "N01N"
  ires=1
  fres=12
  for i in range(ires+1,fres):
    c = select_dihedral_c(U,seg,i)
    chi = chi+str(c.dihedral.value())+" "
  seg = "N01T"
  ires=32
  fres=43
  for i in range(ires+1,fres):
    c = select_dihedral_c(U,seg,i)
    chi = chi+str(c.dihedral.value())+" "
  print(chi,file=outputfile)


outputfile.close()

