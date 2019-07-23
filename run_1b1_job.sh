root_dir=/home/dongxq/Desktop/mol2-file
for(( i=1; i<=20; i++ ))
do
	chmod 777 generate_1b1_job.sh
	./generate_1b1_job.sh $i $root_dir
done
filename=_nlm_1by1_job.pbs
for(( i=1; i<=20; i++ ))
do
	echo $i
	qsub $i$filename
done