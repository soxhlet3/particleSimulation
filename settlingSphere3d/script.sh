#!/bin/sh
cd ${0%/*} || exit 1                        # Run from this directory

make clean

make

mpirun --oversubscribe -np 12  settlingSphere3d | tee -a terminal-output.txt

paste pos_output vel_output force_output > _output #comment out when resolved

g++ resultAnalysis-subgrid.cpp -o resultAnalysis-subgrid
#g++ resultAnalysis-resolved.cpp -o resultAnalysis-resolved

./resultAnalysis-subgrid | tee -a terminal-output.txt
#./resultAnalysis-resolved | tee -a terminal-output.txt

gnuplot subgrid.p
#gnuplot resolved.p

eog v_a_t.png & 



#------------------------------------------------------------------------------
