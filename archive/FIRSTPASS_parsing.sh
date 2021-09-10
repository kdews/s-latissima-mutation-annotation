#for i in $(grep "Ec-" $in_blast | awk '{print $1,$2}' | sort -u)
#do
#	printf "%s\t%s\n" "$i" "Ecto" >> gene.list
#done
#for i in $(grep "PGA" $in_blast | awk '{print $1,$2}' | sort -u)
#do
#	printf "%s\t%s\n" "$i" "Macro" >> gene.list
#done

#for i in $(grep "DE" gene2.list | awk '{print $1}')
#do
#	line=$(grep $'\tgene\t' *gff3 | grep "transcriptId=${i};")
#	coords=$(echo "$line" | awk '{print $1,$4,$5}')
#	printf "%s\t%s\n" "$i" "$coords" >> DE.coords
#done
#for i in $(grep "Spo11" gene2.list | awk '{print $1}')
#do
#	line=$(grep $'\tgene\t' *gff3 | grep "transcriptId=${i}")
#	coords=$(echo "$line" | awk '{print $1,$4,$5}')
#	printf "%s\t%s\n" "$i" "$coords" >> Spo11.coords
#done
#for i in $(grep "Macro\|Ecto" gene2.list | awk '{print $1}')
#do
#	line=$(grep $'\tgene\t' *gff3 | grep "proteinId=${i};")
#	coords=$(echo "$line" | awk '{print $1,$4,$5}')
#	printf "%s\t%s\n" "$i" "$coords" >> > macro_ecto.coords
#done 

# Extract ID nums
#awk '{print $1}' gene.list | awk 'BEGIN {FS="|"} {print $3}' | paste - gene.list > gene2.list

#awk '{print $2,$3,$4}' coords/all.coords | paste gene_lists/gene2.list - > final_gene.list
#sed -i "s/ /\t/g" final_gene.list
