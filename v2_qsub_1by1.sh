#!/bin/bash
### job Name
#PBS -N yfff
### Output Files
#PBS -o 0226.out
#PBS -e 0226.err
### Queue Name
#PBS -q low
### Number of nodes
#PBS -l nodes=1:ppn=24
root_dir=/home/dongxq/Desktop/mol2-file
cd $root_dir/align_code
sastbx.python v2_nprocess_nlm_1by1.py -output $root_dir/align_code/cc_res_10.npy -input $root_dir/align_code/off10 -processnum 24

