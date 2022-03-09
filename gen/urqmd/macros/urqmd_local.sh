#!/bin/sh
# a script that can be used to run a simulation locally
date
export URQMD_DIR="" # URQMD installation path
export OUTDIR=""    # your output path here: path to directory where results will be stored
export inputfile="" # your inputfile for URQMD -- see examples in gen/urqmd/inputfiles
# make a separate directory for results
mkdir $OUTDIR/results
cd $OUTDIR/results
# symlink urqmd executable and essentials
ln -s ${URQMD_DIR}/urqmd # make sure that this corresponds to your executable, e.g. urqmd.x86_64
ln -s ${URQMD_DIR}/eosfiles
ln -s ${URQMD_DIR}/runqmd.bash
ln -s ${OUTDIR}/read_urqmd.C
cp ${inputfile} .
sed -e "s|randomnumber|$RANDOM|" -i inputfile # set random seed for a run
# run URQMD
./runqmd.bash
# remove any residual files
rm urqmd.f13
rm urqmd.f15
rm urqmd.f16
rm urqmd.f19
rm urqmd.f20
rm tables.dat
rm eosfiles
rm runqmd.bash
rm urqmd
date

