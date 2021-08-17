#!/usr/bin/env bash
# Author Jonathan Lombardino, Burton Lab, University of Wisconsin Madison
# make windows from a genome file  
# Usage: bash makewindows.bash mygenome.genome windowsize windowstep

genomefile=$1
name=${genomefile%.*}
windowsfile=""$name"_windows_size"$2"_step"$3".bed"
size=$2
step=$3
bedtools makewindows -g $genomefile -w $2 -s $3 > $windowsfile

