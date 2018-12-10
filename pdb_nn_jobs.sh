root_dir=/home/dongxq/Desktop/mol2-file
filename=_pdb_nn_job.pbs
files=$(ls $1)
for file_pdb in $files;
do
    chmod 777 generate_pdb_nn_job.sh
	pdb_basename=${file_pdb%.*}
    ./generate_pdb_nn_job.sh $pdb_basename $root_dir
	qsub $pdb_basename$filename
done
