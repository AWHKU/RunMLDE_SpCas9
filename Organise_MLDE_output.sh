cd /path/to/MLDE/results/Belper/ExecuteMlde_para1

samplesname=($(ls -d *))
for i in ${samplesname[@]}
do
echo "${i}"
files=($(ls ${i}/))
mv ${i}/${files}/* ${i}/
rm -r ${i}/${files} 
done

cd /path/to/MLDE/results/Belper/ExecuteMlde_para2

samplesname=($(ls -d *))
for i in ${samplesname[@]}
do
echo "${i}"
files=($(ls ${i}/))
mv ${i}/${files}/* ${i}/
rm -r ${i}/${files} 
done

cd /path/to/MLDE/results/Georgiev/ExecuteMlde_para1
samplesname=($(ls -d *))
echo "${samplesname[@]}"
for i in ${samplesname[@]}
do
echo "${i}"
files=($(ls ${i}/))
mv ${i}/${files}/* ${i}/
rm -r ${i}/${files} 
done

cd /path/to/MLDE/results/Georgiev/ExecuteMlde_para2
samplesname=($(ls -d *))
echo "${samplesname[@]}"
for i in ${samplesname[@]}
do
echo "${i}"
files=($(ls ${i}/))
mv ${i}/${files}/* ${i}/
rm -r ${i}/${files} 
done
