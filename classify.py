import argparse
import sys
import re
from collections import defaultdict
import pysam
import mappy as mp


def get_annotation(read_sequence, control_aligner, to_print):
    # Map the read sequence to each of the controls
    default_annotation = ["No_Mapping","0","0.0","0.0","0","0","0","0"]
    for hit in control_aligner.map(read_sequence):
        if hit.is_primary:
            default_annotation = [hit.ctg, str(hit.mapq), str(len(read_sequence)/hit.mlen), str(hit.ctg_len/hit.mlen), str(hit.q_st), str(hit.q_en), str(hit.r_st), str(hit.r_en)]
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

def classify_insertions(input_tsv, control_sequences_file, output_tsv, fastq_list):
    # First load in control_sequences. Assume this is a fasta file
    control_aligner = mp.Aligner(control_sequences_file)
    fastqs = []
    for file in fastq_list.strip().split(','):
        fastqs.append(pysam.FastaFile(filename=file))
    # For each insertion in tsv, read in insert sequence, map it using parasail to control sequences
    # If we reach some min percent identity then consider it mapped
    with open(output_tsv, 'w') as out_tsv:
        with open(input_tsv, 'r') as in_tsv:
            count = 0
            for line in in_tsv:
                if count == 0:
                    count = 1
                    out_tsv.write(line.strip()+"\tAnnotation\tConsensustRTFamily\tFraction\tRTFamily\tScore\tPercentIdentityRead\tPercentIdentityControl\tInsertAnnotationStart\tInsertAnnotationEnd\tAnnotationStart\tAnnotationEnd\n")
                    continue
                row = line.strip().split('\t')
                results = defaultdict(int)
                read_positions = defaultdict(list)
                if len(row) > 5:
                    for read_insert in row[5].split(','):
                        sample = read_insert.split(':')[0]
                        read = read_insert.split(':')[1]
                        orientation = read_insert.split(':')[2]
                        insert_start = int(read_insert.split(':')[3].split('-')[0])
                        insert_end = int(read_insert.split(':')[3].split('-')[1])
                        read_positions[read] = [insert_start, insert_end, orientation]
                else:
                    row.append("NA")
                # Get sequence from reads
                for fastq in fastqs:
                    for read in read_positions:
                        try:
                            seq = fastq.fetch(read)
                            if read_positions[read][2] == '-':
                                seq = reverse_complement(seq)
                            insert_seq = seq[read_positions[read][0]:read_positions[read][1]+1]
                            # Align insert seq
                            result = get_annotation(insert_seq, control_aligner, False)
                            if float(result[3]) < 0.5:
                                results["No_Mapping"] += 1
                            elif "LINE" in result[0]:
                                results["LINE"] += 1
                            elif "SINE" in result[0]:
                                results["SINE"] += 1
                            elif "ERV" in result[0]:
                                results["ERV"] += 1
                            elif "SVA" in result[0]:
                                results["SVA"] += 1
                            else:
                                results["No_Mapping"] += 1
                        except KeyError:
                            # Read not found in this fastq
                            pass
                # Get annotation for representative haplotype sequence
                result = get_annotation(row[4], control_aligner, False)
                # Check result percent identity against minimum
                con_ann, con_frac = get_consensus(results)
                annotation = "PASS"
                if "No_Mapping" in result:
                    annotation = "no_rt_mapping"
                out_tsv.write("\t".join(row)+"\t"+annotation+"\t"+con_ann+"\t"+str(con_frac)+"\t"+"\t".join(result)+"\n")


