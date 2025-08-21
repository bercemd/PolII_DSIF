#!/usr/bin/env python

# B Dutagaci, 2025

import numpy as np
from sklearn.decomposition import PCA
import math
import os,sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '--input_path',type=str, default='coor.dat',
    help='Coordinates for each frame.')
parser.add_argument(
    '--output_path',type=str, default='angle.dat',
    help='Output file for BH bending.')
arg = parser.parse_args()


"""

Bending Angle was calculated between the principal axes of helices with residue selections of 835-852 and 853-870 

"""

def angle_between_vectors(v1, v2):
    uv1 = v1 / np.linalg.norm(v1)  # Normalize v1
    uv2 = v2 / np.linalg.norm(v2)  # Normalize v2
    angle_rad = np.arccos(np.clip(np.dot(uv1, uv2), -1.0, 1.0))  # Calculate angle in radians
    angle_deg = np.degrees(angle_rad)  # Convert to degrees
    return angle_deg,angle_rad

def determine_principal_axis(points):
    centroid = np.mean(points, axis=0)
    centered_points = points - centroid
    covariance_matrix = np.cov(centered_points, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
    max_eigenvalue_index = np.argmax(eigenvalues)
    principal_axis = eigenvectors[:, max_eigenvalue_index]
    return principal_axis

inputfile = open(arg.input_path,"r")
lines = inputfile.readlines()
inputfile.close()

outputfile = open(arg.output_path,"w")
for l in range(len(lines)):
        points1 = []
        for j in range(0,54,3):
              points1.append([float(lines[l].split()[j]),float(lines[l].split()[j+1]),float(lines[l].split()[j+2])])
        pa1 = determine_principal_axis(points1)

        points2 = []

        for j in range(54,len(lines[l].split()),3):
              points2.append([float(lines[l].split()[j]),float(lines[l].split()[j+1]),float(lines[l].split()[j+2])])
        pa2 = determine_principal_axis(points2)

        angle_deg,angle_rad = angle_between_vectors(pa1, pa2)
        angle_deg_cor = angle_deg
        if angle_deg<90:
            angle_deg_cor = 180 - angle_deg

        print(angle_deg, angle_deg_cor, angle_rad,file=outputfile)

outputfile.close()


