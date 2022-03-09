#!/bin/sh
# script that can be used to run a conversion of reconstruction results into
# ROOT trees in a cluster with a batch system
date
export MPDROOT_DIR="" # path to your mpdroot installation
source ${MPDROOT_DIR}/mpdroot_config.sh
export OUTDIR="" # your output path here: path to directory where results will be stored
export INDIR="" # path to reconstruction results
# make a separate directory for each worker
export S=`printf "%02d" ${SLURM_PROCID}`
export D=`printf "%04d" $I`
mkdir -p $OUTDIR/$D/$S
cd $OUTDIR/$D/$S
# symlink ROOT macro
ln -s $OUTDIR/tree.C
ln -s $OUTDIR/conv.C
ln -s $INDIR/$D/$S/rec.root
ln -s $INDIR/$D/$S/sim.root
root -b -q conv.C >> conv.log 2>&1
rm tree.C
rm conv.C
rm rec.root
rm sim.root
rm *.pcm
rm *.d
rm *.so
rm *.txt
rm *.dat
date
