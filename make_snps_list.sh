#!/bin/bash
  
input=$1
genome=$2

while IFS=$'\t' read -r chr pos ID ref alt repeat_type n_peak_calls n_peak_callers mean_BAD mean_SNP_per_segment n_aggregated total_cover es_mean_ref es_mean_alt fdrp_bh_ref fdrp_bh_alt motif_log_pref motif_log_palt motif_fc motif_pos motif_orient motif_conc; do
        left_sequence=$(bedtools getfasta -fi "${genome}" -bed <(echo -e "${chr}\t$((pos-30))\t$((pos-1))") -name | grep -v ">")
        right_sequence=$(bedtools getfasta -fi "${genome}" -bed <(echo -e "${chr}\t$((pos))\t$((pos+30))") -name | grep -v ">")

	echo -e "${ID}\t${left_sequence}[${ref}/${alt}]${right_sequence}"
done < $input | tail -n +2
