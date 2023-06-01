#!/bin/bash

#Go to folder
cd /home/kris/bp_ceroxyloideae/trimal/

#Activate enviroment
source /home/kris/miniconda3/etc/profile.d/conda.sh
conda activate amas


#Calculate AMAS statistics for each gt value. However, if you send the command to run in the whole folder will give you an out-of-range error. Then create a summary for subsets and concatenate them later
cd /home/kris/bp_ceroxyloideae/trimal/

for i in 0.*/*aligned.fasta; do AMAS.py summary -c 16 -f fasta -d dna -i ${i} -o ${i}_summary.txt; done

#Delete the first line in the summary files
for i in 0.*/*summary.txt; do sed '1d' ${i} > ${i}_summary_for_real.txt; done

#Save all summary for real in a statistics for each gt
for i in 0.*; do cat ${i}/*summary_for_real.txt > ${i}/summary_${i}.txt; done

#Move all summary gt to the main trimming folder
mv */*summary_0* .

#Put all information in a new line
for i in summary_0*; do sed -e 's/HEY/\nHEY/g' -e 's/EGU/\nEGU/g' ${i} > ${i}_2.txt; done

#Remove all older summary gt information
rm summary_0.15.txt
rm summary_0.1.txt
rm summary_0.25.txt
rm summary_0.2.txt
rm summary_0.35.txt
rm summary_0.3.txt
rm summary_0.45.txt
rm summary_0.4.txt
rm summary_0.55.txt
rm summary_0.5.txt
rm summary_0.65.txt
rm summary_0.6.txt
rm summary_0.75.txt
rm summary_0.7.txt
rm summary_0.85.txt
rm summary_0.8.txt
rm summary_0.95.txt
rm summary_0.9.txt

mv summary_0.15.txt_2.txt summary_0.15.txt
mv summary_0.1.txt_2.txt summary_0.1.txt
mv summary_0.25.txt_2.txt summary_0.25.txt
mv summary_0.2.txt_2.txt summary_0.2.txt
mv summary_0.35.txt_2.txt summary_0.35.txt
mv summary_0.45.txt_2.txt summary_0.45.txt
mv summary_0.3.txt_2.txt summary_0.3.txt
mv summary_0.4.txt_2.txt summary_0.4.txt
mv summary_0.5.txt_2.txt summary_0.5.txt
mv summary_0.55.txt_2.txt summary_0.55.txt
mv summary_0.65.txt_2.txt summary_0.65.txt
mv summary_0.6.txt_2.txt summary_0.6.txt
mv summary_0.75.txt_2.txt summary_0.75.txt
mv summary_0.7.txt_2.txt summary_0.7.txt
mv summary_0.85.txt_2.txt summary_0.85.txt
mv summary_0.8.txt_2.txt summary_0.8.txt
mv summary_0.95.txt_2.txt summary_0.95.txt
mv summary_0.9.txt_2.txt summary_0.9.txt

#Remove all previous summary
rm */*_summary_for_real.txt
