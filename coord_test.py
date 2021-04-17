coding_region_positions_file = "assembly_ST_collapse_with_short_genes.complete_ORFs.cds_coords"
coding_region_positions_file_no_ext = coding_region_positions_file.rsplit('.', 1)[0]

f=open(coding_region_positions_file,'r')
coords=f.readlines()
f.close()

coords=coords[1:5]

cds_dict = {}
for line in coords:
    cds_line = line.strip().split("\t")
    print(cds_line)
    contig_id = cds_line[0]
    print(contig_id)
    start_position = cds_line[1]
    print(start_position)
    end_position = cds_line[2]
    print(end_position)
    strand = cds_line[3]
    print(strand)
    protein = cds_line[4]
    print(protein)
    cds_list = [start_position, end_position, strand, contig_id]
    cds_dict[protein] = cds_list

