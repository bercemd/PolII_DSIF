#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --partition #specify partition
#SBATCH --mem=16G
#SBATCH --nodelist #specify node 

#Python Set up
export PATH="$HOME/miniconda3/bin:$PATH"

#OpenMM Set up
export OPENMM_PLUGIN_DIR=$HOME/miniconda3/pkgs/openmm-8.0.0-py310h466fd29_0/lib/plugins
export LD_LIBRARY_PATH=$HOME/miniconda3/pkgs/openmm-8.0.0-py310h466fd29_0/lib:$OPENMM_PLUGIN_DIR:$LD_LIBRARY_PATH

# Specify the replicate number
n=1 

# Folders Set up
odir=`pwd`
proddir=$odir/md/full/prod.$n
wdir=$odir/openmm
cd $wdir
equi_prefix=step4_equilibration
equi_out=$proddir/equil

# Equilibration
python -u openmm_run.py -i ${equi_prefix}.1.inp -t toppar.str -p ${odir}/system.full.solv.complete.psf -c ${odir}/system.full.solv.complete.crd -b sysinfo.full.dat -r prot_pos.full.txt -orst ${equi_out}.1.rst -odcd ${equi_out}.1.dcd > ${equi_out}.1.out

for eq_step in {2..8}; do
	prev=`expr ${eq_step} - 1`

	python -u openmm_run.py -i ${equi_prefix}.${eq_step}.inp -t toppar.str -p ${odir}/system.full.solv.complete.psf -c ${odir}/system.full.solv.complete.crd -b sysinfo.full.dat -r prot_pos.full.txt -irst ${equi_out}.${prev}.rst -orst ${equi_out}.${eq_step}.rst -odcd ${equi_out}.${eq_step}.dcd > ${equi_out}.${eq_step}.out
done
exit 0
