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
init=$odir/system.full.solv.complete
#equi_prefix=step4_equilibration
equi_out=$proddir/equil.8
prod_prefix=step5_production

# Production

cnt=1
cntmax=2
while [ ${cnt} -le ${cntmax} ]
do
    pcnt=`echo ${cnt} | awk '{print $1-1}'`
    istep=$proddir/prod_${cnt}
    pstep=$proddir/prod_${pcnt}
    if [ ${cnt} == 1 ];then
      pstep=${equi_out}
    fi
    if [ ! -r ${istep}.dcd ]; then
      python -u openmm_run.py -i ${prod_prefix}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -irst ${pstep}.rst -orst ${istep}.rst -odcd ${istep}.dcd > ${istep}.out
    fi
    cnt=`echo ${cnt} | awk '{print $1+1}'`
    #@ cnt += 1

done

exit 0

