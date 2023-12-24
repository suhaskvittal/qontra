#!/bin/sh
#
# INTENDED FOR USE WITH GEORGIA TECH'S PACE COMPUTING CLUSTER.

nodes=$1
proc=$2
walltime=$3

salloc -A gts-mqureshi4-rg -qinferno -N $nodes --ntasks-per-node $proc -t$walltime
