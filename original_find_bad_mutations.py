# Takes VCF file (filtered for SNPs), CDS coordinate file (format: contig_id, first coord., last coord., strand),
# and transcriptome or genome FASTA file, and outputs effect of SNPs in annotated VCF file

import sys
import math
from Bio import SeqIO
import gzip


# USAGE #
if len(sys.argv) < 2:
    print("Usage: python " + sys.argv[0] + " input.cds_coords input.fasta input.vcf" + "\n" +
          "Requires: BioPython")
    sys.exit(0)


# BASIC OBJECTS AND FUNCTIONS #
# define codon dictionary
codon_dict = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}


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
    comp = []  # generate empty list
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
        fasta = open(filename, 'r')
        fasta.read()  # trying to read zipped file will produce an error
        fasta.close()
    except:
        with gzip.open(filename, 'rt') as reads:
            for read in SeqIO.parse(reads, 'fasta'):
                fasta_dict[read.id] = read.seq
    else:
        with open(filename, 'r') as reads:
            for read in SeqIO.parse(reads, 'fasta'):
                fasta_dict[read.id] = read.seq
    return fasta_dict


# INPUT FILES #
# CDS coordinate file
coding_region_positions_file = sys.argv[1]
# FASTA file
fasta_file = sys.argv[2]
# VCF file
vcf_file = sys.argv[3]
# save VCF name without extension for output filename
vcf_file_no_ext = vcf_file.rsplit('.', 1)[0]


# PARSE INPUT #
# parse CDS coordinate file
f = open(coding_region_positions_file, 'r')
lines = f.readlines()
f.close()
# define dictionary of contig IDs and CDS coordinates + strandedness
contig_id_dict = {}
for line in lines:
    contig_line = line.strip().split(" ")[0]
    contig_id = contig_line.split()[0]
    start_position = contig_line.split()[1]
    end_position = contig_line.split()[2]
    strand = contig_line.split()[3]
    protein = contig_line.split()[4]
    contig_id_list = []
    contig_id_list.append(start_position)
    contig_id_list.append(end_position)
    contig_id_list.append(strand)
    contig_id_list.append(protein)
    contig_id_dict[contig_id] = contig_id_list


# parse FASTA file
seq_dict = readfasta(fasta_file)


# parse VCF file
g = open(vcf_file, 'r')
lines = g.readlines()
g.close()

# OUTPUT FILE #
# open output annotated VCF file
out_vcf = open(vcf_file_no_ext + ".annot.vcf", 'w')


# save coding SNP statistics
coding = 0
non_cod = 0
pos = 0
neg = 0
first = 0
second = 0
third = 0
# begin annotating by iterating through input VCF and saving each line to output
for line in lines:
    # save header lines to output VCF file
    if line.lstrip().startswith('#'):
        out_vcf.write(line)
        continue
    # split non-header lines of VCF file into a list + named components
    if not line.lstrip().startswith('#'):
        vcf_line = line.strip().split(" ")[0]
        vcf_line = vcf_line.split()
        vcf_contig_id = vcf_line[0]
        vcf_snp_position = int(vcf_line[1])
        vcf_ref_base = vcf_line[3]
        vcf_alt_base = vcf_line[4]
        # verify that contig has CDS
        if vcf_contig_id in contig_id_dict.keys():
            # verify that SNP position is within CDS
            if int(vcf_snp_position) >= int(contig_id_dict[vcf_contig_id][0]):
                if int(vcf_snp_position) <= int(contig_id_dict[vcf_contig_id][1]):
                    coding += 1
                    # handles CDSs translated in forward direction
                    if contig_id_dict[vcf_contig_id][2] == "+":
                        pos += 1
                        # subtract SNP position from start position (lower #) to find SNP position in codon
                        # add one because subtracting sequences
                        snp_pos_in_cds = vcf_snp_position - int(contig_id_dict[vcf_contig_id][0]) + 1
                        # handles SNPs in first position in codon
                        if snp_pos_in_cds % 3 == 1:
                            first += 1
                            # all positions must be corrected for Python numbering (starts at 0, so subtract 1)
                            # reference codon sequence
                            codon = seq_dict[vcf_contig_id][vcf_snp_position - 1] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position + 1]
                            # alternate codon sequence
                            alt_codon = vcf_alt_base + \
                                seq_dict[vcf_contig_id][vcf_snp_position] + \
                                seq_dict[vcf_contig_id][vcf_snp_position + 1]
                            if codon_dict[alt_codon] == codon_dict[codon]:
                                vcf_line[7] = vcf_line[7] + ";EFF=SYNONYMOUS(SILENT|" + \
                                              codon + "/" + alt_codon + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            elif codon_dict[alt_codon] == "*":
                                vcf_line[7] = vcf_line[7] + ";EFF=STOP_GAINED(NONSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil(snp_pos_in_cds / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            else:
                                vcf_line[7] = vcf_line[7] + ";EFF=NON_SYNONYMOUS_CODING(MISSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil((int(contig_id_dict[vcf_contig_id][1]) -
                                                  int(contig_id_dict[vcf_contig_id][0])) / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                        # handles SNPs in second position in codon
                        elif snp_pos_in_cds % 3 == 2:
                            second += 1
                            # reference codon sequence
                            codon = seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position - 1] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position]
                            # alternate codon sequence
                            alt_codon = seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                                vcf_alt_base + \
                                seq_dict[vcf_contig_id][vcf_snp_position]
                            if codon_dict[alt_codon] == codon_dict[codon]:
                                vcf_line[7] = vcf_line[7] + ";EFF=SYNONYMOUS(SILENT|" + \
                                              codon + "/" + alt_codon + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            elif codon_dict[alt_codon] == "*":
                                vcf_line[7] = vcf_line[7] + ";EFF=STOP_GAINED(NONSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil(snp_pos_in_cds / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            else:
                                vcf_line[7] = vcf_line[7] + ";EFF=NON_SYNONYMOUS_CODING(MISSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil((int(contig_id_dict[vcf_contig_id][1]) -
                                                  int(contig_id_dict[vcf_contig_id][0])) / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                        # handles SNPs in third position in codon
                        elif snp_pos_in_cds % 3 == 0:
                            third += 1
                            # reference codon sequence
                            codon = seq_dict[vcf_contig_id][vcf_snp_position - 3] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position - 1]
                            # alternate codon sequence
                            alt_codon = seq_dict[vcf_contig_id][vcf_snp_position - 3] + \
                                seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                                vcf_alt_base
                            if codon_dict[alt_codon] == codon_dict[codon]:
                                vcf_line[7] = vcf_line[7] + ";EFF=SYNONYMOUS(SILENT|" + \
                                              codon + "/" + alt_codon + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            elif codon_dict[alt_codon] == "*":
                                vcf_line[7] = vcf_line[7] + ";EFF=STOP_GAINED(NONSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil(snp_pos_in_cds / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            else:
                                vcf_line[7] = vcf_line[7] + ";EFF=NON_SYNONYMOUS_CODING(MISSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil((int(contig_id_dict[vcf_contig_id][1]) -
                                                  int(contig_id_dict[vcf_contig_id][0])) / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                        # alerts user to error
                        else:
                            print("broken")
                    # handles CDSs translated in reverse direction
                    elif contig_id_dict[vcf_contig_id][2] == "-":
                        neg += 1
                        # subtract SNP position from start position (higher #) to find SNP position in codon
                        # add one because subtracting sequences
                        snp_pos_in_cds = int(contig_id_dict[vcf_contig_id][1]) - vcf_snp_position + 1
                        # handles SNPs in first position in codon
                        if snp_pos_in_cds % 3 == 1:
                            first += 1
                            # reference codon sequence
                            codon = seq_dict[vcf_contig_id][vcf_snp_position - 1] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position - 3]
                            codon = comp(codon)
                            # alternate codon sequence
                            alt_codon = vcf_alt_base + \
                                seq_dict[vcf_contig_id][vcf_snp_position - 2] + \
                                seq_dict[vcf_contig_id][vcf_snp_position - 3]
                            alt_codon = comp(alt_codon)
                            if codon_dict[alt_codon] == codon_dict[codon]:
                                vcf_line[7] = vcf_line[7] + ";EFF=SYNONYMOUS(SILENT|" + \
                                              codon + "/" + alt_codon + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            elif codon_dict[alt_codon] == "*":
                                vcf_line[7] = vcf_line[7] + ";EFF=STOP_GAINED(NONSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil(snp_pos_in_cds / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            else:
                                vcf_line[7] = vcf_line[7] + ";EFF=NON_SYNONYMOUS_CODING(MISSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil((int(contig_id_dict[vcf_contig_id][1]) -
                                                  int(contig_id_dict[vcf_contig_id][0])) / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                        # handles SNPs in second position in codon
                        elif snp_pos_in_cds % 3 == 2:
                            second += 1
                            # reference codon sequence
                            codon = seq_dict[vcf_contig_id][vcf_snp_position] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position - 1] + \
                                    seq_dict[vcf_contig_id][vcf_snp_position - 2]
                            codon = comp(codon)
                            # alternate codon sequence
                            alt_codon = seq_dict[vcf_contig_id][vcf_snp_position] + \
                                vcf_alt_base + \
                                seq_dict[vcf_contig_id][vcf_snp_position - 2]
                            alt_codon = comp(alt_codon)
                            if codon_dict[alt_codon] == codon_dict[codon]:
                                vcf_line[7] = vcf_line[7] + ";EFF=SYNONYMOUS(SILENT|" + \
                                              codon + "/" + alt_codon + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            elif codon_dict[alt_codon] == "*":
                                vcf_line[7] = vcf_line[7] + ";EFF=STOP_GAINED(NONSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil(snp_pos_in_cds / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            else:
                                vcf_line[7] = vcf_line[7] + ";EFF=NON_SYNONYMOUS_CODING(MISSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil((int(contig_id_dict[vcf_contig_id][1]) -
                                                  int(contig_id_dict[vcf_contig_id][0])) / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                        # handles SNPs in third position in codon
                        elif snp_pos_in_cds % 3 == 0:
                            third += 1
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
                            if codon_dict[alt_codon] == codon_dict[codon]:
                                vcf_line[7] = vcf_line[7] + ";EFF=SYNONYMOUS(SILENT|" + \
                                              codon + "/" + alt_codon + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            elif codon_dict[alt_codon] == "*":
                                vcf_line[7] = vcf_line[7] + ";EFF=STOP_GAINED(NONSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil(snp_pos_in_cds / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                            else:
                                vcf_line[7] = vcf_line[7] + ";EFF=NON_SYNONYMOUS_CODING(MISSENSE|" + \
                                              codon + "/" + alt_codon + "|" + \
                                              codon_dict[codon] + str(math.ceil(snp_pos_in_cds / 3)) + \
                                              codon_dict[alt_codon] + "|" + \
                                              str(math.ceil((int(contig_id_dict[vcf_contig_id][1]) -
                                                  int(contig_id_dict[vcf_contig_id][0])) / 3)) + "|" + \
                                              vcf_contig_id + "|" + "CODING" + "|" + \
                                              contig_id_dict[vcf_contig_id][3] + "|)"
                                vcf_line = "\t".join(vcf_line)
                                out_vcf.write(vcf_line + "\n")
                        # alerts user to error
                        else:
                            print("broken")
                    else:
                        print("broken")
                else:
                    non_cod += 1
                    vcf_line[7] = vcf_line[7] + ";EFF=NON_CODING"
                    vcf_line = "\t".join(vcf_line)
                    out_vcf.write(vcf_line + "\n")
            else:
                non_cod += 1
                vcf_line[7] = vcf_line[7] + ";EFF=NON_CODING"
                vcf_line = "\t".join(vcf_line)
                out_vcf.write(vcf_line + "\n")
        else:
            non_cod += 1
            vcf_line[7] = vcf_line[7] + ";EFF=NON_CODING"
            vcf_line = "\t".join(vcf_line)
            out_vcf.write(vcf_line + "\n")
    else:
        print("broken")


# close output VCF file
out_vcf.close()
# print SNP statistics to stdout
print("SNP Statsitics" + "\n" +
      "Non-coding: " + str(non_cod) + "\n" +
      "Coding: " + str(coding) + "\n" +
      "% Coding: " + str(round(coding/(non_cod+coding)*100, 2)) + "%\n\n" +
      "Positive strand: " + str(pos) + "\n" +
      "Negative strand: " + str(neg) + "\n\n" +
      "Position in codon:\n" +
      "1st - " + str(first) + "\n" +
      "2nd - " + str(second) + "\n" +
      "3rd - " + str(third))

