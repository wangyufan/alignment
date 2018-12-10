cat > $1_pdb_nl_job.pbs << END_TEXT
#!/bin/bash
### job Name
#PBS -N yf_nn
### Output Files
#PBS -o nl.out
#PBS -e nl.err
### Queue Name
#PBS -q low
### Number of nodes
#PBS -l nodes=1:ppn=24 
cd $2/align_code
sastbx.python nl.py -processnum 24 -dbpath $2/align_code/myDB/myDB -target_shape $2/pdb_test/$1.pdb -target_type pdb -output $2/pdb_out_nl/output_nn_$1 --nmax 20
