#copy raw data into processing folder so we don't alter raw data
cp -r Raw_Data/. Processed_Data

#insert country labels to be used for identification later
sed -i 's/>/>IND_/g' Processed_Data/IND_flu_reads.fa
sed -i 's/>/>AUS_/g' Processed_Data/AUS_flu_reads.fa
sed -i 's/>/>CAN_/g' Processed_Data/CAN_flu_reads.fa

#rename files so as to not confuse with raw data
for f in Processed_Data/*.fa; do
    mv -- "$f" "${f%.fa}_labeled.fa"
done

#create reference file needed for alignment
touch Processed_Data/flu_ref.fa
awk "/^>/ {n++} n>1 {exit} {print}" Processed_Data/IND_flu_reads_labeled.fa >> Processed_Data/flu_ref.fa

#create 1 file with all reads needed for alignment
touch Processed_Data/all_flu_reads.fa 
cat Processed_Data/CAN_flu_reads_labeled.fa >> Processed_Data/all_flu_reads.fa
cat Processed_Data/AUS_flu_reads_labeled.fa >> Processed_Data/all_flu_reads.fa
cat Processed_Data/IND_flu_reads_labeled.fa >> Processed_Data/all_flu_reads.fa