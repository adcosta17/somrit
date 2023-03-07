import argparse
import sys
import re
from collections import defaultdict
import pysam
import mappy as mp
from intervaltree import Interval, IntervalTree


def get_annotation(read_sequence, control_aligner, to_print):
    # Map the read sequence to each of the controls
    default_annotation = ["No_Mapping", "+", "0","0.0","0.0","0","0","0","0"]
    for hit in control_aligner.map(read_sequence):
        if hit.is_primary:
            strand = '+'
            if hit.strand < 0:
                strand = '-'
            default_annotation = [hit.ctg, strand, str(hit.mapq), str(len(read_sequence)/hit.mlen), str(hit.ctg_len/hit.mlen), str(hit.q_st), str(hit.q_en), str(hit.r_st), str(hit.r_en)]
            if to_print:
                print(hit.ctg)
            return default_annotation
    return default_annotation

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def get_consensus(result_dict):
    # Determine the conensus
    max_count = 0
    max_annotation = "No_Mapping"
    total_count = result_dict["No_Mapping"] + result_dict["LINE"] + result_dict["LINE"] + result_dict["SVA"] + result_dict["ERV"]
    if result_dict["LINE"] > max_count:
        max_annotation = "LINE"
        max_count = result_dict["LINE"]
    if result_dict["SINE"] > max_count:
        max_annotation = "SINE"
        max_count = result_dict["SINE"]
    if result_dict["SVA"] > max_count:
        max_annotation = "SVA"
        max_count = result_dict["SVA"]
    if result_dict["ERV"] > max_count:
        max_annotation = "ERV"
        max_count = result_dict["ERV"]
    if total_count > 0:
        return (max_annotation, max_count/total_count)
    else:
        return (max_annotation, 0)

def classify_insertions(tsv_list, control_sequences_file, output_tsv, fastq_list, realign_tsv, samples):
    seen = defaultdict(IntervalTree)
    realigned_reads = defaultdict(int)
    if realign_tsv is not None:
        with open(realign_tsv, 'r') as in_tsv:
            count = 0
            for line in in_tsv:
                if count == 0:
                    count = 1
                    continue
                row = line.strip().split('\t')
                if len(row) < 6:
                    # Concensus not supported by any reads. Likely computed but not more favorable than the original reference
                    continue
                chrom = row[0]
                start = int(row[1])
                end =  int(row[2])
                for read_insert in row[5].split(','):
                    read = read_insert.split(':')[1]
                    seen[chrom][start:end] = read
                    realigned_reads[read] = 1
    for i in range(len(tsv_list.split(','))):
        file = tsv_list.split(',')[i]
        with open(file, 'r') as in_tsv:
            count = 0
            for line in in_tsv:
                if count == 0:
                    count = 1
                    continue
                row = line.strip().split('\t')
                chrom = row[0]
                start = int(row[1]) - 1000
                end =  int(row[2]) + 1000
                if row[3] in realigned_reads:
                    # Have re-aligned this read
                    continue
                else:
                    #print(line.strip())
                    seen[chrom][int(row[1]):int(row[2])] = row[3]
    # Load in control_sequences. Assume this is a fasta file
    control_aligner = mp.Aligner(control_sequences_file)
    fastqs = []
    for file in fastq_list.strip().split(','):
        fastqs.append(pysam.FastaFile(filename=file))
    # For each insertion in tsv, read in insert sequence, map it using minimap2 to control sequences
    # If we reach some min percent identity then consider it mapped
    with open(output_tsv, 'w') as out_tsv:
        # Classify Re-aligned inserts first and then do individual ones that did not get re-aligned
        out_tsv.write("Chromosome\tStart\tEnd\tName\tInsertSequence\tSupportingReads\tAnnotation\tRTFamily\tRTOrientation\tScore\tPercentIdentityRead\tPercentIdentityConsensus\tInsertAnnotationStart\tInsertAnnotationEnd\tAnnotationStart\tAnnotationEnd\n")
        if realign_tsv is not None:
            with open(realign_tsv, 'r') as in_tsv:
                count = 0
                for line in in_tsv:
                    if count == 0:
                        count = 1
                        continue
                    row = line.strip().split('\t')
                    if len(row) < 6:
                        # Concensus not supported by any reads. Likely computed but not more favorable than the original reference
                        continue
                    result = get_annotation(row[4], control_aligner, False)
                    annotation = "PASS"
                    if "No_Mapping" in result:
                        annotation = "no_rt_mapping"
                    chrom = row[0]
                    start = int(row[1])
                    end =  int(row[2])
                    out_tsv.write("\t".join(row)+"\t"+annotation+"\t"+"\t".join(result)+"\n")
        for i in range(len(tsv_list.split(','))):
            file = tsv_list.split(',')[i]
            sample = samples[i]
            with open(file, 'r') as in_tsv:
                count = 0
                for line in in_tsv:
                    if count == 0:
                        count = 1
                        continue
                    row = line.strip().split('\t')
                    chrom = row[0]
                    start = int(row[1]) - 1000
                    end =  int(row[2]) + 1000
                    if row[3] in realigned_reads:
                        # Have re-aligned this read already
                        continue
                    result = get_annotation(row[7], control_aligner, False)
                    annotation = row[8]
                    if "No_Mapping" in result:
                        if annotation == "PASS":
                            annotation = "no_rt_mapping"
                        else: 
                            annotation += ",no_rt_mapping"
                    out_tsv.write(row[0]+"\t"+row[1]+"\t"+row[2]+"\t"+row[3]+":"+row[4]+"-"+row[5]+"\t"+row[7]+"\t"+sample+":"+row[3]+":"+row[6]+":"+row[4]+"-"+row[5]+"\t"+annotation+"\t"+"\t".join(result)+"\n")


