#!/bin/bash

#Charmm Executable Set up
export CHARMMEXEC="$HOME/charmm/bin/charmm -chsize 800000"

#Modeller Set up
export MODINSTALL10v1="$HOME/modeller-10.1"
export KEY_MODELLER10v1="MODELIRANJE"
export LIBS_LIB10v1="$MODINSTALL10v1/modlib/libs.lib"
alias  mod="mod10.1"
export MODELLEREXEC="$HOME/modeller-10.1/bin/mod10.1_x86_64-intel8"
export LD_LIBRARY_PATH="$HOME/modeller-10.1/lib/x86_64-intel8":$LD_LIBRARY_PATH

#Python Set up
export PATH="$HOME/miniconda3/bin:$PATH"
export PYTHONPATH="$HOME/miniconda3/lib/python3.7"
export PYTHONHOME="$HOME/miniconda3/lib"

#MMTSB Set up
export MMTSBDIR="$HOME/mmtsb"
export CHARMMDATA=$MMTSBDIR/data/charmm
export PATH=$PATH:$MMTSBDIR/perl:$MMTSBDIR/bin


# Organize PDB file: add segment names, convert ZN to ZN2 for forcefield 
convpdb.pl -segnames 5oik.pdb > 5oik.seg.pdb
sed -i 's/ZN  ZN /ZN  ZN2/g' 5oik.seg.pdb
sed -i 's/HETATM/ATOM  /g' 5oik.seg.pdb

# Make the RNA chain R for MMTSB to recognize
sed 's/N01P/N01R/g' 5oik.seg.pdb | convpdb.pl -readseg -chainfromseg > 5oik.seg.chain.pdb

# Make W segment for ions manually for CHARMM to handle

# Generate missing residues of Trigger Loop ##########
loopModel.pl -models 5 -loop A1103:TLNTFHYAGVSA 5oik.seg.chain.pdb
# It will generate 5 models, follow the next step using the model with the minimum energy.

# Segnames for generated model
convpdb.pl -segnames model.5.pdb > model5.seg.pdb

############# FOR -DSIF SYSTEM ###############

# Remove subunits SPT4 and SPT5 (chains Y and Z)
convpdb.pl -nsel !Y+Z: -readseg model.5.pdb > spt4_5.ex.seg.pdb
# Continue next step with spt4_5.ex.seg.pdb for -DSIF system

##############################################

# Let CHARMM build H atoms
enerCHARMM.pl -par xtop=toppar/top_all36_na.rtf,pstream=toppar/toppar_all36_na_nad_ppi.str,parflex,terlist=N01R:5TER:3TRM,hsp=A1108,noperiodic -custom custompdb.inp -log ener.log -cmd ener.cmd model5.seg.pdb
convpdb.pl -chainfromseg system.pdb > system.chain.pdb

# Rotate to optimum orientation found by orient.py
convpdb.pl system.chain.pdb -readseg -rotatex 0 | convpdb.pl -readseg -rotatey 30 | convpdb.pl -readseg -rotatez 30 > system.rotate.pdb

# Solvate the box
convpdb.pl -solvate -cubic -cutoff 9 system.rotate.pdb > system.solv.pdb

# Calculate the charge and add ions to neutralize
enerCHARMM.pl -par param=x,xtop=toppar/top_all36_prot.rtf:toppar/top_all36_na.rtf,xpar=toppar/par_all36m_prot.prm:toppar/par_all36_na.prm,pstream=toppar/toppar_all36_na_nad_ppi.str,parflex,terlist=N01R:5TER:3TRM,hsp=A1108,noperiodic -charge -log charge.log system.pdb

convpdb.pl -ions SOD:110 system.solv.pdb > system.solv.ions.pdb

# Organize water for CHARMM to generate psf file later

convpdb.pl -nsel water -segnames system.solv.ions.pdb > system.solv.water.pdb
convpdb.pl -nsel ^water -segnames system.solv.ions.pdb > system.solv.nowater.pdb
convpdb.pl -readseg -renumber 1 system.solv.water.pdb > system.solv.water.renumber.pdb
./water.sh system.solv.water.renumber.pdb
mv temp.pdb system.solv.water.mod.pdb

cat system.solv.nowater.pdb system.solv.water.mod.pdb | grep -v END > system.full.solv.complete.pdb
echo END >> system.full.solv.complete.pdb

# Recenter the molecule for OpenMM
enerCHARMM.pl -par param=x,xtop=toppar/c36m.hv/c36m_hv.rtf:toppar/c36m.hv/c36m_hv_na.rtf,xpar=toppar/c36m.hv/c36m_hv.prm:toppar/c36m.hv/c36m_hv_na.prm,pstream=toppar/c36m.hv/c36m_hv_na_nad_ppi.str,parflex,terlist=N01R:5TER:3TRM,hsp=A1108 -custom center.full.inp -log center.full.log system.full.solv.complete.pdb

# Generate psf
genPSF.pl -par param=x,xtop=toppar/c36m.hv/c36m_hv.rtf:toppar/c36m.hv/c36m_hv_na.rtf,xpar=toppar/c36m.hv/c36m_hv.prm:toppar/c36m.hv/c36m_hv_na.prm,pstream=toppar/c36m.hv/c36m_hv_na_nad_ppi.str,parflex,terlist=N01R:5TER:3TRM,hsp=A1108 -log genpsf.log -crdout system.full.solv.complete.crd system.full.solv.complete.center.pdb > system.full.solv.complete.psf


