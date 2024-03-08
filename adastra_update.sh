#!/bin/bash

file_list='/home/subpolare/adastra/bill_cipher_TFs.txt'
genome='/home/subpolare/genome/GRCh38.primary_assembly.genome.fa'
output='/home/subpolare/adastra/adastra_pwm_0.01.tsv'
p_value=0.001
threads=10
discr=1000

echo -e "Factor\tAdastra_concordance\tConcordance\tAdastra_discordance\tDiscordance\tAdastra_hit\tHit\tAdastra_ho_hit\tNo_hit\t0\t1\t2\t3" > adastra_pwm_0.tsv
echo -e "Factor\tAdastra_concordance\tConcordance\tAdastra_discordance\tDiscordance\tAdastra_hit\tHit\tAdastra_ho_hit\tNo_hit\t0\t1\t2\t3" > adastra_pwm_1.tsv
echo -e "Factor\tAdastra_concordance\tConcordance\tAdastra_discordance\tDiscordance\tAdastra_hit\tHit\tAdastra_ho_hit\tNo_hit\t0\t1\t2\t3" > adastra_pwm_2.tsv
echo -e "Factor\tAdastra_concordance\tConcordance\tAdastra_discordance\tDiscordance\tAdastra_hit\tHit\tAdastra_ho_hit\tNo_hit\t0\t1\t2\t3" > adastra_pwm_3.tsv

mkdir pwm_results_0
mkdir pwm_results_1
mkdir pwm_results_2
mkdir pwm_results_3
mkdir filtered

process_file() {
    tf=$1
    index=$2
    discr=$3
    p_value=$4
    genome=$5
    factor=$(basename ${tf} | cut -f1 -d '.' | cut -f1 -d '_')

    scripts/asb.py ${tf} 0.05 ${factor}
    scripts/make_snps_list.sh filtered/${factor}.filtered ${genome} > ${factor}.snps

    java -cp scripts/ape-3.0.6.jar ru.autosome.perfectosape.SNPScan hocomoco/v12/pwm/${factor}.H12RSNP.${index}.*.pwm ${factor}.snps --single-motif -F 1 -P 1 -d ${discr} > pwm_results_${index}/${factor}.perfectos

    conc_adastra=$(grep -E 'Concordant' filtered/${factor}.filtered | wc -l)
    disc_adastra=$(grep -E 'Discordant' filtered/${factor}.filtered | wc -l)
    nohit_adastra=$(grep -E 'No Hit' filtered/${factor}.filtered | wc -l)
    scripts/concordance.py pwm_results_${index}/${factor}.perfectos ${tf} ${conc_adastra} ${disc_adastra} ${nohit_adastra} ${factor} ${p_value} >> adastra_pwm_${index}.tsv
}

export -f process_file

parallel -j $threads process_file ::: $(cat ${file_list}) ::: 0 ::: ${discr} ::: ${p_value} ::: ${genome}
parallel -j $threads process_file ::: $(cat ${file_list}) ::: 1 ::: ${discr} ::: ${p_value} ::: ${genome}
parallel -j $threads process_file ::: $(cat ${file_list}) ::: 2 ::: ${discr} ::: ${p_value} ::: ${genome}
parallel -j $threads process_file ::: $(cat ${file_list}) ::: 3 ::: ${discr} ::: ${p_value} ::: ${genome}

scripts/adastra_pwm.py ${output}

rm *.snps
find pwm_results_* -maxdepth 1 -type f -size 0 -delete
