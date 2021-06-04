import sys
import math
import gzip
from Bio import SeqIO

# SANDBOX PARAMETERS
# "data/assembly_ST_collapse_with_short_genes.complete_ORFs.cds_coords"
# "data/assembly_ST_collapse_with_short_genes.fa"
# "data/s_latissima_wgs_50_subset_FILT_qual20_SNPs_biallelic.calls.recode.vcf"

# USAGE #
if len(sys.argv) < 2 or sys.argv[1] == "-h":
    print("\n-------------------------FIND DELETERIOUS MUTATIONS IN VARIATION DATA-------------------------\n\n"
          "Takes a CDS coordinate file (format: contig_id, first coord., last coord., strand, protein_id), \n"
          "a transcriptome or genome FASTA file, and a VCF file, and outputs effect \n"
          "of variants in an annotated VCF file (original filename + annot.vcf) with \"EFF=...\" appended to \n"
          "the INFO column, as well as a summary statistics file (substitution_statistics.txt).\n\n" +
          "Usage: python " + sys.argv[0] + " input.cds_coords input.fasta input.vcf" + "\n\n" +
          "Requires: BioPython")
    sys.exit(0)

# BASIC OBJECTS AND FUNCTIONS #
# define codon dictionary
codon_dict = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W"}

# define amino acid dictionaries of volume and polarity
all_amino_acids = ["G", "P", "A", "V", "L", "I", "M", "C", "F", "Y",
                   "W", "H", "K", "R", "Q", "N", "E", "D", "S", "T"]
high_vol = ["L", "I", "F", "M", "Y", "W", "H", "K", "R", "E", "Q"]
low_vol = [AA for AA in all_amino_acids if AA not in high_vol]
polar = ["Y", "W", "H", "K", "R", "E", "Q", "T", "D", "N", "S", "C"]
nonpolar = [AA for AA in all_amino_acids if AA not in polar]
volume_dict = dict.fromkeys(high_vol, "high volume")
volume_dict.update(dict.fromkeys(low_vol, "low volume"))
polarity_dict = dict.fromkeys(polar, "polar")
polarity_dict.update(dict.fromkeys(nonpolar, "nonpolar"))


def compbase(base):
    """Returns complement of base"""
    if base == "A":
        compbase = "T"
    elif base == "C":
        compbase = "G"
    elif base == "G":
        compbase = "C"
    elif base == "T":
        compbase = "A"
    else:
        compbase = base  # if there's an N in the sequence for example, leave it alone
    return compbase


def comp(seq):
    """Returns the complement of seq"""
    # generate empty list
    comp = []
    for base in seq:
        nuc = compbase(base)
        comp.append(nuc)
    cp = "".join(comp)
    # return reverse complement
    return cp


def readfasta(filename):
    """Reads in a FASTA file (zipped or unzipped) and converts to a dictionary
    where the read ID is the key and the sequence is the value"""
    fasta_dict = {}
    try:
        fasta = open(filename, "r")
        # trying to read zipped file will produce an error
        fasta.read()
        fasta.close()
    except:
        with gzip.open(filename, "rt") as reads:
            for read in SeqIO.parse(reads, "fasta"):
                fasta_dict[read.id] = read.seq
    else:
        with open(filename, "r") as reads:
            for read in SeqIO.parse(reads, "fasta"):
                fasta_dict[read.id] = read.seq
    return fasta_dict


def non_coding_line(vcf_line):
    """Takes a line of input VCF file (vcf_line) and returns annotated VCF line
    (See: # DEFINE OUTPUT FILES #)"""
    # appends effect information to the VCF "INFO" field
    annot_vcf_line = vcf_line
    annot_vcf_line[7] = annot_vcf_line[7] + ";EFF=NON_CODING"
    return annot_vcf_line


def synon_line(codon, alt_codon, codon_dict,
               snp_pos_in_cds, protein,
               vcf_contig_id, vcf_line):
    """Takes values for codon and alt_codon (in nucleic acid),
    a line of input VCF file (vcf_line) and returns annotated VCF line
    (See: # DEFINE OUTPUT FILES #)"""
    # appends effect information to the VCF "INFO" field
    annot_vcf_line = vcf_line
    annot_vcf_line[7] = annot_vcf_line[7] + ";EFF=SYNONYMOUS(LOW|SILENT|" + \
                  codon + "/" + alt_codon + "|" + \
                  codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                  codon_dict[alt_codon] + "|" + \
                  str(math.ceil(snp_pos_in_cds / 3)) + "|" + \
                  vcf_contig_id + "|" + "|" + "CODING" + "|" + \
                  protein + "|)"
    return annot_vcf_line


def nonsynon_line(codon, alt_codon, codon_dict,
                  cds_dict, snp_pos_in_cds, protein,
                  vcf_contig_id, vcf_line):
    """Takes values for codon and alt_codon (in nucleic acid),
        a codon to amino acid dictionary (codon_dict), the position of the SNP
        in the CDS (snp_pos_in_cds), a CDS dictionary (cds_dict), a line of input
        VCF file (vcf_line) and returns annotated VCF line
        (See: # DEFINE OUTPUT FILES #)"""
    # annotate radical and conservative amino acid changes
    if codon_dict[codon] == "*" or codon_dict[alt_codon] == "*":
        sub_class = "CONSERVATIVE"
    elif volume_dict[codon_dict[codon]] == volume_dict[codon_dict[alt_codon]] and \
            polarity_dict[codon_dict[codon]] == polarity_dict[codon_dict[alt_codon]]:
        sub_class = "CONSERVATIVE"
    else:
        sub_class = "RADICAL"
    # appends effect information to the VCF "INFO" field
    annot_vcf_line = vcf_line
    annot_vcf_line[7] = annot_vcf_line[7] + ";EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|" + \
                  codon + "/" + alt_codon + "|" + \
                  codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                  codon_dict[alt_codon] + "|" + \
                  str(math.ceil((int(cds_dict[protein][1]) -
                                 int(cds_dict[protein][0])) / 3)) + "|" + \
                  vcf_contig_id + "|" + "|" + "CODING" + "|" + \
                  protein + "|" + sub_class + "|)"
    nonsynon_info = [annot_vcf_line, sub_class]
    return nonsynon_info


def nonsense_line(codon, alt_codon, codon_dict,
                  snp_pos_in_cds, protein,
                  vcf_contig_id, vcf_line):
    """Takes values for codon and alt_codon (in nucleic acid),
    a codon to amino acid dictionary (codon_dict),
    a CDS dictionary (cds_dict), the position of the
    SNP in the CDS (snp_pos_in_cds), a line of input
    VCF file (vcf_line) and returns annotated VCF line
    (See: # DEFINE OUTPUT FILES #)"""
    # appends effect information to the VCF "INFO" field
    annot_vcf_line = vcf_line
    annot_vcf_line[7] = annot_vcf_line[7] + ";EFF=STOP_GAINED(HIGH|NONSENSE|" + \
                  codon + "/" + alt_codon + "|" + \
                  codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                  codon_dict[alt_codon] + "|" + \
                  str(math.ceil(snp_pos_in_cds / 3)) + "|" + \
                  vcf_contig_id + "|" + "|" + "CODING" + "|" + \
                  protein + "|)"
    return annot_vcf_line


def deletion_line(vcf_line):
    """Takes a line of input VCF file (vcf_line) and returns annotated VCF line
     (See: # DEFINE OUTPUT FILES #)"""
    # appends effect information to the VCF "INFO" field
    annot_vcf_line = vcf_line
    annot_vcf_line[7] = annot_vcf_line[7] + ";EFF=DELETION"
    return annot_vcf_line


def insertion_line(vcf_line):
    """Takes a line of input VCF file (vcf_line) and returns annotated VCF line
     (See: # DEFINE OUTPUT FILES #)"""
    # appends effect information to the VCF "INFO" field
    annot_vcf_line = vcf_line
    annot_vcf_line[7] = annot_vcf_line[7] + ";EFF=INSERTION"
    return annot_vcf_line


def coding_line(codon_dict, cds_dict, seq_dict,
                protein, vcf_contig_id, vcf_snp_position,
                vcf_alt_base, vcf_line):
    """Takes values for a codon to amino acid dictionary (codon_dict),
    a CDS dictionary (cds_dict), a line of input VCF file (vcf_line)
    as well as its associated protein and contig IDs (protein, vcf_contig_id),
    SNP position (vcf_snp_position) and alternative base call (vcf_alt_base),
    and returns annotated VCF line
    (See: # DEFINE OUTPUT FILES #)"""
    global nonsynon_info
    snp_info = []
    # handles CDSs translated in forward direction
    if cds_dict[protein][2] == "+":
        # add 1 because subtracting sequences
        snp_pos_in_cds = vcf_snp_position - int(cds_dict[protein][0]) + 1
        if snp_pos_in_cds % 3 == 1:
            snp_info.append("first")
            # all positions must be corrected for Python numbering (starts at 0, so subtract 1)
            # reference codon sequence
            codon = seq_dict[vcf_contig_id][vcf_snp_position - 1] + \
                    seq_dict[vcf_contig_id][vcf_snp_position] + \
                    seq_dict[vcf_contig_id][vcf_snp_position + 1]
            # alternate codon sequence
            alt_codon = vcf_alt_base + \
                        seq_dict[vcf_contig_id][vcf_snp_position] + \
                        seq_dict[vcf_contig_id][vcf_snp_position + 1]
        elif snp_pos_in_cds % 3 == 2:
            snp_info.append("second")
            codon = seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                    seq_dict[vcf_contig_id][vcf_snp_position - 1] + \
                    seq_dict[vcf_contig_id][vcf_snp_position]
            alt_codon = seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                        vcf_alt_base + \
                        seq_dict[vcf_contig_id][vcf_snp_position]
        elif snp_pos_in_cds % 3 == 0:
            snp_info.append("third")
            codon = seq_dict[vcf_contig_id][vcf_snp_position - 3] + \
                    seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                    seq_dict[vcf_contig_id][vcf_snp_position - 1]
            alt_codon = seq_dict[vcf_contig_id][vcf_snp_position - 3] + \
                        seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                        vcf_alt_base
        else:
            print("Error calculating SNP position.")
    elif cds_dict[protein][2] == "-":
        # subtract SNP position from start position (higher #) to find SNP position in codon
        # add 1 because subtracting sequences
        snp_pos_in_cds = int(cds_dict[protein][1]) - vcf_snp_position + 1
        # handles SNPs in first position in codon
        if snp_pos_in_cds % 3 == 1:
            snp_info.append("first")
            codon = seq_dict[vcf_contig_id][vcf_snp_position - 1] + \
                    seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                    seq_dict[vcf_contig_id][vcf_snp_position - 3]
            codon = comp(codon)
            alt_codon = vcf_alt_base + \
                        seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                        seq_dict[vcf_contig_id][vcf_snp_position - 3]
            alt_codon = comp(alt_codon)
        # handles SNPs in second position in codon
        elif snp_pos_in_cds % 3 == 2:
            snp_info.append("second")
            codon = seq_dict[vcf_contig_id][vcf_snp_position] + \
                    seq_dict[vcf_contig_id][vcf_snp_position - 1] + \
                    seq_dict[vcf_contig_id][vcf_snp_position - 2]
            codon = comp(codon)
            alt_codon = seq_dict[vcf_contig_id][vcf_snp_position] + \
                        vcf_alt_base + \
                        seq_dict[vcf_contig_id][vcf_snp_position - 2]
            alt_codon = comp(alt_codon)
        elif snp_pos_in_cds % 3 == 0:
            snp_info.append("third")
            # reference codon sequence
            codon = seq_dict[vcf_contig_id][vcf_snp_position + 1] + \
                    seq_dict[vcf_contig_id][vcf_snp_position] + \
                    seq_dict[vcf_contig_id][vcf_snp_position - 1]
            codon = comp(codon)
            # alternate codon sequence
            alt_codon = seq_dict[vcf_contig_id][vcf_snp_position + 1] + \
                        seq_dict[vcf_contig_id][vcf_snp_position] + \
                        vcf_alt_base
            alt_codon = comp(alt_codon)
        else:
            print("Error calculating SNP position.")
    else:
        print("Error detecting strand direction (+/- expected).")
    if codon_dict[alt_codon] == codon_dict[codon]:
        snp_info.append("synonymous")
        annot_vcf_line = synon_line(codon, alt_codon, codon_dict,
                                    snp_pos_in_cds, protein,
                                    vcf_contig_id, vcf_line)
    elif codon_dict[alt_codon] == "*":
        snp_info.append("nonsense")
        annot_vcf_line = nonsense_line(codon, alt_codon, codon_dict,
                                             snp_pos_in_cds, protein,
                                             vcf_contig_id, vcf_line)
    else:
        snp_info.append("nonsynonymous")
        subclass = nonsynon_info[1]
        snp_info.append(subclass)
        nonsynon_info = nonsynon_line(codon, alt_codon, codon_dict,
                                             cds_dict, snp_pos_in_cds, protein,
                                             vcf_contig_id, vcf_line)
        annot_vcf_line = nonsynon_info[0]
    snp_info = [annot_vcf_line, snp_info]
    return snp_info


def get_stats(coding, first, second, third, synon, nonsynon, nonsense, deletion, insertion, radical, conservative):
    """Writes variant statistics to stats output file (out_stats)
    (See: # DEFINE OUTPUT FILES #)"""
    stats = str("Variant Statistics" + "\n\n" +
                "Total SNPs: " + str(snp) + "\n" +
                "Non-coding: " + str(non_cod) + "\n" +
                "Coding: " + str(coding) + "\n" +
                "% Coding: " + str(round(coding / (non_cod + coding) * 100, 2)) + "%\n\n" +
                "Synonymous substitutions: " + str(synon) + "\n" +
                "Non-synonymous substitutions: " + str(nonsynon) + "\n" +
                "Nonsense (early stop) substitutions: " + str(nonsense) + "\n" +
                "Insertions: " + str(insertion) + "\n" +
                "Deletions: " + str(deletion) + "\n"
                "Radical: " + str(radical) + "\n"
                "Conservative: " + str(conservative) + "\n\n" +
                "Position in codon:\n" +
                "1st - " + str(first) + "\n" +
                "2nd - " + str(second) + "\n" +
                "3rd - " + str(third) + "\n")
    return stats


# INPUT FILES #
# CDS coordinate file
coding_region_positions_file = sys.argv[1]
# FASTA file
fasta_file = sys.argv[2]
# VCF file
vcf_file = sys.argv[3]
# save VCF name without extension for output filename
vcf_file_no_ext = vcf_file.rsplit(".", 1)[0]

# PARSE INPUT #
# parse CDS coordinate file
f = open(coding_region_positions_file, "r")
coords = f.readlines()
f.close()
# define dictionary of protein IDs keys with values of
# CDS coordinates, strandedness, and contig IDs
# **keep in mind, ONLY contigs with a CDS are in this dictionary**
cds_dict = {}
for line in coords:
    cds_line = line.strip().split()
    contig_id = cds_line[0]
    start_position = cds_line[1]
    end_position = cds_line[2]
    strand = cds_line[3]
    protein = cds_line[4]
    cds_list = [start_position, end_position, strand, contig_id]
    cds_dict[protein] = cds_list
# parse FASTA file
seq_dict = readfasta(fasta_file)
# parse VCF file
g = open(vcf_file, "r")
vcf = g.readlines()
g.close()

# DEFINE OUTPUT FILES #
# open output annotated VCF file
out_vcf = open(vcf_file_no_ext + ".annot.vcf", "w")
out_stats = "VCF_statistics.txt"

# VCF ANNOTATION PIPELINE #
# create counters to save coding SNP statistics
snp = 0
coding = 0
non_cod = 0
first = 0
second = 0
third = 0
synon = 0
nonsynon = 0
nonsense = 0
deletion = 0
insertion = 0
radical = 0
conservative = 0
# begin annotating by iterating through input VCF and saving each line to output
for line in vcf:
    # save header lines to output VCF file without editing
    if line.lstrip().startswith("#"):
        out_vcf.write(line)
        continue
    # split non-header lines of VCF file into a list + named components
    else:
        snp += 1
        vcf_line = line.strip().split("\t")
        vcf_contig_id = vcf_line[0]
        vcf_snp_position = int(vcf_line[1])
        vcf_ref_base = vcf_line[3]
        vcf_alt_base = vcf_line[4]
        annot_vcf_line = []
        # write deletion line to annotated VCF and go back to the start of the loop
        if len(vcf_ref_base) > 1:
            deletion += 1
            annot_vcf_line = deletion_line(vcf_line)
            continue
        # create a list of alt_bases to iterate through and check for protein effect in each of those mutations
        alt_base_list = vcf_line[4].strip().split(",")
        # verify that contig has CDS
        if any(vcf_contig_id in s for s in cds_dict.keys()):
            protein_list = [s for s in cds_dict.keys() if vcf_contig_id in s]
            correct_proteins = []
            for protein in protein_list:
                # verify that SNP position is within a CDS
                if int(cds_dict[protein][0]) <= int(vcf_snp_position) <= int(cds_dict[protein][1]):
                    # make list to hold protein IDs of potentially overlapping protein CDSs
                    correct_proteins.append(protein)
            # define annot_vcf_line and update SNP position/effect information for SNP_statistics output
            if correct_proteins:
                coding += 1
                # if the variant is in more than one protein CDS, it will iterate through each protein ID
                # and add as many annotations to the SAME annotated VCF line as there are protein IDs
                # in the format: "EFF=protein1;EFF=protein2;..."
                for protein in correct_proteins:
                    for vcf_alt_base in alt_base_list:
                        # if alt base is bigger than one, means it's an insertion ex: ref == A alt == ATC
                        if len(vcf_alt_base) > 1:
                            insertion += 1
                            annot_vcf_line = insertion_line(vcf_line)
                        elif vcf_alt_base == "*":
                            deletion += 1
                            annot_vcf_line = deletion_line(vcf_line)
                        else:
                            snp_info = coding_line(codon_dict, cds_dict, seq_dict,
                                                 protein, vcf_contig_id,
                                                 vcf_snp_position, vcf_alt_base,
                                                 vcf_line)
                            annot_vcf_line = snp_info[0]
                            if snp_info[1][0] == "first":
                                first += 1
                            elif snp_info[1][0] == "second":
                                second += 1
                            elif snp_info[1][0] == "third":
                                third += 1
                            if snp_info[1][1] == "synonymous":
                                synon += 1
                            elif snp_info[1][1] == "nonsynonymous":
                                nonsynon += 1
                                if snp_info[2] == "RADICAL":
                                    radical += 1
                                elif snp_info[2] == "CONSERVATIVE":
                                    conservative += 1
                            elif snp_info[1][1] == "nonsense":
                                nonsense += 1
                    # set vcf_line to annot_vcf_line
                    # (only used if iterating through again, e.g., variant is in multiple protein CDSs)
                    vcf_line = annot_vcf_line
                    if len(correct_proteins) > 1:
                        print(annot_vcf_line)
            else:
                non_cod += 1
                annot_vcf_line = non_coding_line(vcf_line)
        # else annotate as non-coding mutation
        else:
            non_cod += 1
            annot_vcf_line = non_coding_line(vcf_line)
        if annot_vcf_line:
            annot_vcf_line = "\t".join(annot_vcf_line)
            out_vcf.write(annot_vcf_line + "\n")
        else:
            print("Error creating annotated VCF line.")
        # write updated stats to out_stats summary file after every 100th line
        if snp % 100 == 0:
            h = open(out_stats, "w")
            h.write(get_stats(coding, first, second, third, synon, nonsynon, nonsense, deletion, insertion))
            h.close()
            # print(get_stats(coding, first, second, third, synon, nonsynon, nonsense, deletion, insertion))


# close output VCF file
out_vcf.close()

# Filtering out duplicates that come up with multiallelic calls
fileIn = open(vcf_file_no_ext + ".annot.vcf", "r")
out_vcf = fileIn.readlines()
fileIn.close()

fileOut = open(vcf_file_no_ext + "_annot_filtered.vcf", "w")

counter = 0
for line in out_vcf:
    if line.lstrip().startswith("#"):
        fileOut.write(line)
        continue
    elif counter == 0:
        previous_vcf_line_full = line
        vcf_line = line.strip().split("\t")
        previous_vcf_contig_id = vcf_line[0]
        previous_vcf_snp_position = int(vcf_line[1])
        counter += 1
    else:
        vcf_line = line.strip().split("\t")
        vcf_contig_id = vcf_line[0]
        vcf_snp_position = int(vcf_line[1])
        if vcf_contig_id == previous_vcf_contig_id:
            if int(vcf_snp_position) == int(previous_vcf_snp_position):
                previous_vcf_line_full = line
            else:
                fileOut.write(previous_vcf_line_full)
                previous_vcf_line_full = line
                previous_vcf_contig_id = vcf_contig_id
                previous_vcf_snp_position = vcf_snp_position
fileOut.write(line)
fileOut.close()
### INTERMEDIATE FILE THAT IS NOT CORRECTED FOR COPIES CAN BE DELETED HERE!! ###
