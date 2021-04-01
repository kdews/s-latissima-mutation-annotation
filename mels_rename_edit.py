import os,re,sys
import subprocess

# This script renames all files in a folder to preferred format.

## the files we want to rename are in the format:
## alpha_NanAmp_Sac1_species_sampleid_plate_bases_Sac2_sequencer_lane_read.fastq.gz
## examples:
### IZXU_NanoAmplified_Saccharina_latissima_SL-JS-19-FG-2_1_TCCAACGC_Saccharina_I997_L1_R1.fastq.gz
### IZXU_NanoAmplified_Saccharina_latissima_SL-JS-19-FG-2_1_TCCAACGC_Saccharina_I997_L1_R2.fastq.gz
### IZYN_NanoAmplified_Saccharina_angustissima_SA-CB-5-MG-3_1_GCAATGCA_Saccharina_I997_L1_R1.fastq.gz
### IZYN_NanoAmplified_Saccharina_angustissima_SA-CB-5-MG-3_1_GCAATGCA_Saccharina_I997_L1_R2.fastq.gz

## Some files (sequenced by machine I1018 or I1019) have a _num_ before sampleid

### Renamed format:
## sampleid_sequencer_plate_lane_read.fastq.gz

print("Creating list of all FASTQs in folder...")
cmd = "ls *fastq.gz > original_filenames.txt"
os.system(cmd)

print("Renaming...")

with open('original_filenames.txt','r') as files:
    s = list()
    for file in files:
        file_stripped = file.strip()
        alts = ['I1018','I1019']
        if any(item in file_stripped for item in alts):
            sampleid = file_stripped.rsplit('_')[5]
            sequencer = file_stripped.rsplit('_')[9]
#            plate = file_stripped.rsplit('_')[6]
            lane = file_stripped.rsplit('_')[10]
            read_fq = file_stripped.rsplit('_')[11]
        else:
            sampleid = file_stripped.rsplit('_')[4]
            sequencer = file_stripped.rsplit('_')[8]
#            plate = file_stripped.rsplit('_')[5]
            lane = file_stripped.rsplit('_')[9]
            read_fq = file_stripped.rsplit('_')[10]
        newname = sampleid + '_' + sequencer + '_' + lane + '_' + read_fq
        s.append(newname)
        cmd = 'mv ' + file_stripped + ' ' + newname
        os.system(cmd)
    files.close()

## Define function to convert a list to string using join() function
def list2string(mylist):
    # initialize an empty string
    str1 = '\n'
    # return the new joined string
    return str1.join(mylist)

## Print new filenames to text file in same order as original filenames
with open('new_filenames.txt', 'w') as newfilenames:
    newfilenames.write(list2string(s))  # convert list to string before writing to file


