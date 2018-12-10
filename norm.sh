for i in {1..2}  
do  
outputfile = "c$i.data'
awk '{print $4,$5}' $i |awk -F')' '{print $1}' |awk -F',' '{print $1,$2}' > outputfile
sed -i '$d' $outputfile
sed -i '$d' $outputfile
sed -i '$d' $outputfile
sed -i '1,2d' $outputfile
done  
