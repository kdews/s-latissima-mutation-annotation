import sys
import math
import gzip
from Bio import SeqIO

# USAGE #
if len(sys.argv) < 2 or sys.argv[1] == "-h":
    print("\n-------------------------FIND DELETERIOUS MUTATIONS IN VARIATION DATA-------------------------\n\n"
          "Takes a CDS coordinate file (format: contig_id, first coord., last coord., strand, protein_id), \n"
          "a transcriptome or genome FASTA file, and a VCF file (filtered for SNPs), and outputs effect \n"
          "of SNPs in an annotated VCF file (original filename + annot.vcf) with \"EFF=...\" appended to \n"
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


def write_non_coding_line(vcf_line, out_vcf):
    """Takes a line of input VCF file (vcf_line) and an opened
     output VCF file (out_vcf), annotates the VCF line, and
     writes the annotated line to the output VCF
     (See: # DEFINE OUTPUT FILES #)"""
    # appends effect information to the VCF "INFO" field
    vcf_line[7] = vcf_line[7] + ";EFF=NON_CODING"
    vcf_line = "\t".join(vcf_line)
    out_vcf.write(vcf_line + "\n")


def write_synon_line(codon, alt_codon,
                     vcf_line, out_vcf):
    """Takes values for codon and alt_codon (in nucleic acid),
    a line of input VCF file (vcf_line), and an opened output VCF file
    (out_vcf), annotates the VCF line, and
    writes the annotated line to the output VCF
    (See: # DEFINE OUTPUT FILES #)"""
    # appends effect information to the VCF "INFO" field
    vcf_line[7] = vcf_line[7] + ";EFF=SYNONYMOUS(SILENT|" + \
                  codon + "/" + alt_codon + "|)"
    vcf_line = "\t".join(vcf_line)
    out_vcf.write(vcf_line + "\n")


def write_nonsynon_line(codon, alt_codon, codon_dict,
                        cds_dict, snp_pos_in_cds, protein,
                        vcf_contig_id, vcf_line, out_vcf):
    """Takes values for codon and alt_codon (in nucleic acid),
        a codon to amino acid dictionary (codon_dict), the position of the SNP
        in the CDS (snp_pos_in_cds), a CDS dictionary (cds_dict), a line of input
        VCF file (vcf_line), and an opened output VCF file (out_vcf), annotates
        the VCF line, and writes the annotated line to the output VCF
        (See: # DEFINE OUTPUT FILES #)"""
    # appends effect information to the VCF "INFO" field
    vcf_line[7] = vcf_line[7] + ";EFF=NON_SYNONYMOUS_CODING(MISSENSE|" + \
                  codon + "/" + alt_codon + "|" + \
                  codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                  codon_dict[alt_codon] + "|" + \
                  str(math.ceil((int(cds_dict[protein][1]) -
                                 int(cds_dict[protein][0])) / 3)) + "|" + \
                  vcf_contig_id + "|" + "CODING" + "|" + \
                  cds_dict[protein][3] + "|)"
    vcf_line = "\t".join(vcf_line)
    out_vcf.write(vcf_line + "\n")


def write_nonsense_line(codon, alt_codon, codon_dict,
                        cds_dict, snp_pos_in_cds, protein,
                        vcf_contig_id, vcf_line, out_vcf):
    """Takes values for codon and alt_codon (in nucleic acid),
    a codon to amino acid dictionary (codon_dict),
    a CDS dictionary (cds_dict), the position of the
    SNP in the CDS (snp_pos_in_cds), a line of input
    VCF file (vcf_line), and an opened output VCF file
    (out_vcf), annotates the VCF line, and writes
    the annotated line to the output VCF
    (See: # DEFINE OUTPUT FILES #)"""
    # appends effect information to the VCF "INFO" field
    vcf_line[7] = vcf_line[7] + ";EFF=STOP_GAINED(NONSENSE|" + \
                  codon + "/" + alt_codon + "|" + \
                  codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                  codon_dict[alt_codon] + "|" + \
                  str(math.ceil(snp_pos_in_cds / 3)) + "|" + \
                  vcf_contig_id + "|" + "CODING" + "|" + \
                  cds_dict[protein][3] + "|)"
    vcf_line = "\t".join(vcf_line)
    out_vcf.write(vcf_line + "\n")


def write_coding_line(codon_dict, cds_dict, seq_dict,
                      protein, vcf_contig_id, vcf_snp_position,
                      vcf_alt_base, vcf_line, out_vcf):
    """Takes values for a codon to amino acid dictionary (codon_dict),
    a CDS dictionary (cds_dict), a line of input VCF file (vcf_line)
    as well as its associated protein and contig IDs (protein, vcf_contig_id),
    SNP position (vcf_snp_position) and alternative base call (vcf_alt_base),
    and an opened output VCF file (out_vcf), annotates the VCF line, and writes
    the annotated line to the output VCF
    (See: # DEFINE OUTPUT FILES #)"""
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
        write_synon_line(codon, alt_codon, vcf_line, out_vcf)
    elif codon_dict[alt_codon] == "*":
        snp_info.append("nonsense")
        write_nonsense_line(codon, alt_codon, codon_dict,
                            cds_dict, snp_pos_in_cds, protein,
                            vcf_contig_id, vcf_line, out_vcf)
    else:
        snp_info.append("nonsynonymous")
        write_nonsynon_line(codon, alt_codon, codon_dict,
                            cds_dict, snp_pos_in_cds, protein,
                            vcf_contig_id, vcf_line, out_vcf)
    return snp_info


def get_stats(coding, first, second, third, synon, nonsynon, nonsense):
    """Writes SNP statistics to stats output file (out_stats)
    (See: # DEFINE OUTPUT FILES #)"""
    stats = str("SNP Statistics" + "\n\n" +
                 "Total SNPs: " + str(snp) + "\n" +
                 "Non-coding: " + str(non_cod) + "\n" +
                 "Coding: " + str(coding) + "\n" +
                 "% Coding: " + str(round(coding / (non_cod + coding) * 100, 2)) + "%\n\n" +
                 "Synonymous substitutions: " + str(synon) + "\n" +
                 "Non-synonymous substitutions: " + str(nonsynon) + "\n" +
                 "Nonsense (early stop) substitutions: " + str(nonsense) + "\n\n" +
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
    cds_line = line.strip().split("\t")
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
out_stats = "substitution_statistics.txt"


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
        # verify that contig has CDS
        if any(vcf_contig_id in s for s in cds_dict.keys()):
            protein_list = [s for s in cds_dict.keys() if vcf_contig_id in s]
            for protein in protein_list:
                # verify that SNP position is within a CDS
                if int(cds_dict[protein][0]) <= int(vcf_snp_position) <= int(cds_dict[protein][1]):
                    correct_protein = protein
                else:
                    correct_protein = "none"
            if correct_protein != "none":
                coding += 1
                # write to VCF and return a list (snp_info) containing SNP position and effect information
                snp_info = write_coding_line(codon_dict, cds_dict, seq_dict,
                                             correct_protein, vcf_contig_id,
                                             vcf_snp_position, vcf_alt_base,
                                             vcf_line, out_vcf)
                if snp_info[0] == "first":
                    first += 1
                elif snp_info[0] == "second":
                    second += 1
                elif snp_info[0] == "third":
                    third += 1
                if snp_info[1] == "synonymous":
                    synon += 1
                elif snp_info[1] == "nonsynonymous":
                    nonsynon += 1
                elif snp_info[1] == "nonsense":
                    nonsense += 1
            else:
                non_cod += 1
                write_non_coding_line(vcf_line, out_vcf)
        # else annotate as non-coding mutation
        else:
            non_cod += 1
            write_non_coding_line(vcf_line, out_vcf)
        # write updated stats to out_stats summary file with each iteration
        h = open(out_stats, "w")
        h.write(get_stats(coding, first, second, third, synon, nonsynon, nonsense))
        h.close()

# close output VCF file
out_vcf.close()
