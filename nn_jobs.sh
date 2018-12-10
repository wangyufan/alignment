#dblist=(dude1-db dude2-db dude-actives)
root_dir=/home/dongxq/Desktop/mol2-file
for i in $@;
do
    chmod 777 generate_nn_job.sh
    ./generate_nn_job.sh $i $root_dir
done
filename=_nn_job.pbs
for i in $@;
do
    echo $i
    qsub $i$filename
done

