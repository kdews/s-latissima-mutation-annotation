coding_region_positions_file = "assembly_ST_collapse_with_short_genes.complete_ORFs.cds_coords"
coding_region_positions_file_no_ext = coding_region_positions_file.rsplit('.', 1)[0]

f=open(coding_region_positions_file,'r')
lines=f.readlines()
f.close()

contig_id_dict = {}


for line in lines:
    contig_line = line.strip().split(" ")[0]
    contig_id = contig_line.split()[0]
    start_position = contig_line.split()[1]
    end_position = contig_line.split()[2]
    strand = contig_line.split()[3]
    contig_id_list = []
    contig_id_list.append(start_position)
    contig_id_list.append(end_position)
    contig_id_list.append(strand)
    print(contig_id_list)
    contig_id_dict[contig_id] = contig_id_list    
print(contig_id_dict)

