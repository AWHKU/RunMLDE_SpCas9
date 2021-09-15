python3 2020_01_22_EncodingGenerator_learned.py \
'bepler' SaCas9 \
--fasta /ssd/users/athchu/SaCas9_MLDE_input.fasta \
--out /ssd/users/athchu/SaCas9_full_length_bepler

python3 2020_05_14_EncodingGenerator_georgiev.py \
'georgiev' SaCas9 \
--fasta /ssd/users/athchu/SaCas9_MLDE_input.fasta \
--out /ssd/users/athchu/SaCas9_full_length_georgiev


wkdir=/path/to/Belper_encoding/
cd /path/to/SpCas9_diverse_variants
fitness_files=($(ls *max*))
echo "${fitness_files[@]}"
cd ~/MLDE
for i in ${fitness_files[@]}
do
mkdir ${wkdir}/ExecuteMlde_para1/${i/.csv/}
echo "${i/.csv/}"
python ExecuteMlde.py \
 /path/to/SpCas9_diverse_variants/${i} \
 ${wkdir}/Encodings/Cas9_bepler_Normalized.npy  \
 ${wkdir}/Encodings/Cas9_bepler_ComboToIndex.pkl \
 --model_params ${wkdir}/ExecuteMlde_para1/testMldeParam.csv \
 --output ${wkdir}/ExecuteMlde_para1/${i/.csv/}/ \
 --hyperopt 
done

for i in ${fitness_files[@]}
do
mkdir ${wkdir}/ExecuteMlde_para2/${i/.csv/}
echo "${i/.csv/}"
python ExecuteMlde.py \
 /path/to/SpCas9_diverse_variants/${i} \
 ${wkdir}/Encodings/Cas9_bepler_Normalized.npy  \
 ${wkdir}/Encodings/Cas9_bepler_ComboToIndex.pkl \
 --model_params ${wkdir}/ExecuteMlde_para2/testMldeParam.csv \
 --output ${wkdir}/ExecuteMlde_para2/${i/.csv/}/ \
 --hyperopt 
done


wkdir=/path/to/Georgiev_encoding
cd /path/to/SpCas9_diverse_variants
fitness_files=($(ls *fitness*))
echo "${fitness_files[@]}"
cd ~/MLDE

for i in ${fitness_files[@]}
do
mkdir ${wkdir}/ExecuteMlde_para1/${i/.csv/}
echo "${i/.csv/}"
python ExecuteMlde.py \
 /path/to/SpCas9_diverse_variants/${i} \
 ${wkdir}/Encodings/Cas9_georgiev_Normalized.npy  \
 ${wkdir}/Encodings/Cas9_georgiev_ComboToIndex.pkl \
 --model_params ${wkdir}/ExecuteMlde_para1/testMldeParam.csv \
 --output ${wkdir}/ExecuteMlde_para1/${i/.csv/}/ \
 --hyperopt 
done

for i in ${fitness_files[@]}
do
mkdir ${wkdir}/ExecuteMlde_para2/${i/.csv/}
echo "${i/.csv/}"
python ExecuteMlde.py \
 /path/to/SpCas9_diverse_variants/${i} \
 ${wkdir}/Encodings/Cas9_georgiev_Normalized.npy  \
 ${wkdir}/Encodings/Cas9_georgiev_ComboToIndex.pkl \
 --model_params ${wkdir}/ExecuteMlde_para2/testMldeParam.csv \
 --output ${wkdir}/ExecuteMlde_para2/${i/.csv/}/ \
 --hyperopt 
done
