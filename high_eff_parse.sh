


grep "#" master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.ann.gene_list.recode.vcf >> master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.ann.gene_list.HIGH_EFF.recode.vcf
grep "|HIGH|" master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.ann.gene_list.recode.vcf >> master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.ann.gene_list.HIGH_EFF.recode.vcf
grep -v "##" master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.ann.gene_list.HIGH_EFF.recode.vcf | grep "#" | awk '{$3=$6=$7=$8=$9=""; print $0}' > high_eff.tab; for i in $(seq 1 1 $(grep -v "#" master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.ann.gene_list.HIGH_EFF.recode.vcf | wc -l)); do line=$(grep -v "#" master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.ann.gene_list.HIGH_EFF.recode.vcf | sed -n ${i}p); echo "$line" | awk '{$3=$6=$7=$8=$9=""; print $0}'; echo "$line" | awk '{print $8}' | sed 's/.*ANN=//g' | sed "s/,/\n/g" | awk 'BEGIN {FS="|"} {print $2,$7}'; done >> high_eff.tab
