#!/bin/sh
# script that can be used to run a box simulation on a cluster with a batch system
date
export MPDROOT_DIR="" # path to your mpdroot installation
source ${MPDROOT_DIR}/mpdroot_config.sh
export OUTDIR="" # your output path here: path to directory where results will be stored
# make a separate directory for each worker
export S=`printf "%02d" ${SLURM_PROCID}`
export D=`printf "%04d" $I`
export SEED=$S$D # set fixed seeds for each worker for reproducibility
mkdir -p $OUTDIR/$D/$S
cd $OUTDIR/$D/$S
# symlink ROOT macro
ln -s $OUTDIR/sim.C
ln -s $OUTDIR/rec.C
ln -s $OUTDIR/geometry.C
# run Geant simulation and track reconstruction
root -b -q sim.C >> sim.log 2>&1
root -b -q rec.C > rec.log 2>&1
# remove any residual files
rm *.txt
rm *.dat
rm sim.C
rm rec.C
rm geometry.C
date
