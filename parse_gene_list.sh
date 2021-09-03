for i in `grep "Ec-" blast_results/candidate_genes_vs_SlaSLCT1FG3_1_GeneCatalog_proteins_20210608.prot.blast.tab | awk '{print $2}' | sort -u`; do printf "%s\t%s\n" "$i" "Ecto" >> gene.list; done
for i in `grep "PGA" blast_results/candidate_genes_vs_SlaSLCT1FG3_1_GeneCatalog_proteins_20210608.prot.blast.tab | awk '{print $2}' | sort -u`; do printf "%s\t%s\n" "$i" "Macro" >> gene.list; done

sed -i "s/-/|/g" gene.list
awk '{print $1}' gene.list | awk 'BEGIN {FS="|"} ; {print $3}' | paste - gene.list > gene2.list

for i in $(grep "DE" gene2.list | awk '{print $1}'); do line=$(grep $'\tgene\t' *gff3 | grep "transcriptId=${i};"); coords=$(echo "$line" | awk '{print $1,$4,$5}'); printf "%s\t%s\n" "$i" "$coords"; done > DE.coords
for i in $(grep "Spo11" gene2.list | awk '{print $1}'); do line=$(grep $'\tgene\t' *gff3 | grep "transcriptId=${i}"); coords=$(echo "$line" | awk '{print $1,$4,$5}'); printf "%s\t%s\n" "$i" "$coords"; done > Spo11.coords
for i in $(grep "Macro\|Ecto" gene2.list | awk '{print $1}'); do line=$(grep $'\tgene\t' *gff3 | grep "proteinId=${i};"); coords=$(echo "$line" | awk '{print $1,$4,$5}'); printf "%s\t%s\n" "$i" "$coords"; done > macro_ecto.coords
cat DE.coords Spo11.coords macro_ecto.coords > all.coords
sed -i "s/ /\t/g" coords/all.coords


echo "#ID	full_ID	Source	Chr	Start	End" > final_gene.list
awk '{print $2,$3,$4}' coords/all.coords | paste gene_lists/gene2.list - > final_gene.list
sed -i "s/ /\t/g" final_gene.list
