cat > $1_nn_job.pbs << END_TEXT
#!/bin/bash
### job Name
#PBS -N yf_nn
### Output Files
#PBS -o nn.out
#PBS -e nn.err
### Queue Name
#PBS -q low
### Number of nodes
#PBS -l nodes=1:ppn=24 
cd $2/align_code
sastbx.python nn.py -processnum 24 -dbpath $2/align_code/myDB/$1 -target_shape $2/mol2_520/cavity-all.mol2 -target_type mol2 -output $2/output_nn_$1 --nmax 20

