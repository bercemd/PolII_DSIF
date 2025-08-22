# PolII_DSIF
## Investigation of Pol II - DSIF Elongation Complex Dynamics using MD Simulations
This repository contains the data for the simulations of the RNA polymerase II (Pol II) elongation complexes with and without DSIF elongation factor. 

The ```input``` folder provides the scripts to run equilibration (```openmm.equil.full.slurm```) and production (```openmm.prod.full.slurm```) simulations using 
OpenMM simulation package and initial solvated structures of the two systems (```structures.tar.gz```). 

The ```input/system_preparation``` folder contains the bash script (```system_prep.sh```) that have the stepwise instructions to prepare the solvated Pol II systems from the structure with PDB ID: 5OIK, and all the additional scripts that are called by this script. 

The folder ```input/openmm``` contains the Python scripts, which are originated by the CHARMM-GUI server and slightly modified, for running MD simulations. It also contains input parameters for equilibration (```step4_equilibration.#.inp```) and production ```step5_production.inp```) runs. 

The folder ```input/toppar``` contains the force field parameters, which are CHARMM-c36(m) parameters and we modified them to repartition Hydrogen atom masses as described in the preprint provided in the citation section below.

The ```analysis``` folder contains the scripts for the analysis of the simulations under the ```analysis/scripts``` folder and the output data obtained from the analysis under the ```analysis/data``` folder. 

*** Requirements:
```
Python versions 3+

MDAnalysis

Linux environment

CHARMM

MMTSB

OpenMM

MODELLER

SLURM job scheduler
```
*** Usage:

To prepare the initial systems:

```
./system_prep.sh

```
To run equilibration and production simulations:

```
sbatch openmm.equil.full.slurm
sbatch openmm.prod.full.slurm
```
To run analysis scripts:

```
python [options] filename.py
```
An example for running an analysis script:

To calculate RMSD of biomolecules (proteins or nucleic acids):
```
python rmsd_mda.py --traj_path=md.dcd --psf_path=rnap.psf --pdb_path=rnap.pdb --fit_selection='(segid P01A and name CA)' --rmsd_selection='(segid P01A and name CA)' --out_path=rmsd.dat
```
*** Citation:

```
Amith, W. D., Bogart B., Dutagaci B., Molecular Basis for Impacts of DSIF on the Dynamics of RNA Polymerase II Elongation Complex, bioRxiv, doi: https://doi.org/10.1101/2025.08.09.669504, 2025
```
