#!/bin/bash
### job Name
#PBS -N yfff
### Output Files
#PBS -o 121.out
#PBS -e 121.err
### Queue Name
#PBS -q low
### Number of nodes
#PBS -l nodes=1:ppn=24
root_dir=/home/dongxq/Desktop/mol2-file
cd $root_dir/align_code
sastbx.python nprocess_align_onebyone.py -output $root_dir/align_code/np_onbyone -mol2dir $root_dir/align_code/data/split-dude-actives/ -processnum 24

