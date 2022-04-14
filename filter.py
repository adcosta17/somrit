import argparse
import sys
import re
from collections import defaultdict
import parasail
import mappy as mp
from intervaltree import Interval, IntervalTree
import pysam

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
user_matrix = parasail.matrix_create("ACGT", 2, -1)

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def filter_control(sample_dict, control_sample):
    if control_sample is not None:
        if control_sample in sample_dict:
            return "In_Control_Sample"
    return None

def filter_centromere_telomere(row, centromere_positions, telomere_positions):
    if centromere_positions[row[0]][int(row[1]):int(row[2])]:
        return "In_Centromere"
    if telomere_positions[row[0]][int(row[1]):int(row[2])]:
        return "In_Teloomere"
    return None

def filter_insert_map(chrom, start, end, read_positions, sequences, bam_read_positions, ref_aligner, min_mapq):
    # Then for each read that has a secondary mapping other than the insert pos, extract the insert sequence
    if len(bam_read_positions) == 0:
        return None
    read_failed = {}
    for read in read_positions:
        bam_positions = bam_read_positions[read]
        pos = read_positions[read]
        read_start = pos[0]
        read_end = pos[1]
        insert_seq = sequences[read][read_start:read_end+1]
        # Map the insert sequence to the ref and check position. If position overlaps a secondary mapping, flag
        for hit in ref_aligner.map(insert_seq):
            for bam_pos  in bam_positions:
                if bam_pos[0] == chrom and bam_pos[1] <= start and bam_pos[2] >= end:
                    continue
                if hit.ctg == bam_pos[0] and hit.r_st <= bam_pos[1] and hit.r_en >= bam_pos[2]:
                    read_failed[read] = 1
    if len(read_failed) > 0:
        return "In_Secondary_Mapping_"+str(len(read_failed))
    return None

def filter_mapping_qual(chrom, start, end, bam_ref_positions, min_mapq, chrom_list):
    # Look at all the reads that map in the region at the insertion. If more than half are of mapq < 20. Flag
    count = 0
    low_mapq = 0
    if chrom not in bam_ref_positions or chrom not in chrom_list:
        return None
    for pos in bam_ref_positions[chrom]:
        pos_start = pos.split('-')[0]
        pos_end = pos.split('-')[1]
        if len(pos_start) == 0:
            continue
        if len(pos_end) == 0:
            continue
        pos_start = int(pos_start)
        pos_end = int(pos_end)
        if start > pos_start and end < pos_end:
            for item in bam_ref_positions[chrom][pos]:
                if item[2] < start and item[3] > end:
                    count += 1
                    if item[1] < min_mapq:
                        low_mapq += 1
    if count == 0:
        return None
    if low_mapq/count > 0.5:
        return "Low_Mapping_Quality_Region"
    return None

def filter_poly_AT(read_positions, sequences):
    has_poly_AT = 0
    total = 0
    for read in read_positions:
        pos = read_positions[read]
        read_start = pos[0]
        read_end = pos[1]
        insert_seq = sequences[read][read_start:read_end+1]
        # Get first and last 50 bases
        start = insert_seq[:50]
        end = insert_seq[len(insert_seq)-50:]
        pA = "A"*10
        pT = "T"*10
        total += 1
        if pA in start or pT in start or pA in end or pT in end:
            has_poly_AT += 1
    if total == 0:
        return None
    if has_poly_AT/total < 0.5:
        return "Low_Poly_AT_Percentage"
    return None

def check_match(cigar):
    cgs = re.findall('[0-9]*[A-Z=]', cigar)
    for cg in cgs:
        if cg.endswith("M"):
            if int(cg[:cg.find("M")]) >= 5:
                return True
        if cg.endswith("="):
            if int(cg[:cg.find("=")]) >= 5:
                return True
    return False

def filter_TSD(read_positions, sequences):
    tsd = []
    for read in read_positions:
        pos = read_positions[read]
        read_start = pos[0]
        read_end = pos[1]
        if read_start+25 > read_end-25:
            # Have too small an insertion
            continue
        if (read_start-25) < 0 or (read_end -25) < 0 :
            continue
        if (read_start+25) >= len(sequences[read]) or (read_end+25) >= len(sequences[read]):
            continue
        insert_start_seq = sequences[read][read_start-25:read_start+25]
        insert_end_seq = sequences[read][read_end-25:read_end+25]
        # Need to look for duplication between the two
        # Align the two sequences together and check for at least 5 matches in a row as a possible TSD
        result = parasail.sw_trace(insert_start_seq, insert_end_seq, 11, 1, user_matrix)
        cigar = result.cigar.decode.decode('utf-8')
        if check_match(cigar):
            tsd.append(read)
    if len(tsd) == 0:
        return "No_TSD_Found"
    return None

def update_annotation(annotation, filters):
    ret = annotation
    for f in filters:
        if f is not None:
            if ret == "PASS":
                ret = f
            else:
                ret = ","+f
    return ret

def filter_insertions(input_tsv, output_tsv, bam_list, fastq_list, contromeres, telomeres, control_sample, min_mapq, ref, cluster_window, chrs_to_use):
    # Open and read in centromere and telomere positions
    telomere_positions= defaultdict(IntervalTree)
    centromere_positions = defaultdict(IntervalTree)
    chrom_list = chrs_to_use.split(',')
    with open(telomeres) as in_tf:
        count = 0
        for line in in_tf:
            if count == 0:
                count = 1
                continue
            line_args = line.strip().split('\t')
            chrom = line_args[1]
            start = int(line_args[2])
            end = int(line_args[3])
            key = chrom + ":" + str(start) + "-" + str(end)
            telomere_positions[chrom][start:end] = key
    with open(contromeres) as in_cf:
        count = 0
        for line in in_cf:
            if count == 0:
                count = 1
                continue
            line_args = line.strip().split('\t')
            chrom = line_args[1]
            start = int(line_args[2])
            end = int(line_args[3])
            key = chrom + ":" + str(start) + "-" + str(end)
            centromere_positions[chrom][start:end] = key
    fastqs = []
    for file in fastq_list.strip().split(','):
        fastqs.append(pysam.FastaFile(filename=file))
    print("Fastqs")
    all_positions_list = defaultdict(list)
    with open(input_tsv, 'r') as in_tsv:
        count = 0
        for line in in_tsv:
            if count == 0:
                count = 1
                continue
            row = line.strip().split('\t')
            all_positions_list[row[0]].append((int(row[1])-cluster_window, int(row[2])+cluster_window))
    print("All Positions")
    for chrom in all_positions_list:
        sorted_by_lower_bound = sorted(all_positions_list[chrom], key=lambda tup: tup[0])
        merged = []
        for higher in sorted_by_lower_bound:
            if not merged:
                merged.append(higher)
            else:
                lower = merged[-1]
                # test for intersection between lower and higher:
                # we know via sorting that lower[0] <= higher[0]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    merged[-1] = (lower[0], upper_bound)  # replace by merged interval
                else:
                    merged.append(higher)
        all_positions_list[chrom] = merged
    print("Merged")
    all_positions = defaultdict(IntervalTree)
    for chrom in all_positions_list:
        for pos in all_positions_list[chrom]:
            start = pos[0]
            end = pos[1]
            all_positions[chrom][start:end] = 1
    print("IntervalTree")
    # Parse through the list of bams, get the positions of all alignments
    bam_read_positions = defaultdict(list)
    bam_ref_positions = {}
    for bam in bam_list.split(','):
        reader = pysam.AlignmentFile(bam)
        for record in reader.fetch():
            if record.mapping_quality >= min_mapq:
                bam_read_positions[record.query_name].append([record.reference_name, record.reference_start, record.reference_end])
            nearby = all_positions[record.reference_name][record.reference_start:record.reference_end]
            if len(nearby) > 0:
                # Have a record in a region we care about
                if record.reference_name not in bam_ref_positions:
                    bam_ref_positions[record.reference_name] = defaultdict(list)
                for item in nearby:
                    #print(record.query_name+" "+record.reference_name+":"+str(item.begin)+"-"+str(item.end))
                    bam_ref_positions[record.reference_name][str(item.begin)+"-"+str(item.end)].append([record.query_name, record.mapping_quality, record.reference_start, record.reference_end])
    print("Setup Dicts")
    ref_aligner = mp.Aligner(ref)
    # Read in tsv and then run each filter
    with open(input_tsv, 'r') as in_tsv:
        with open(output_tsv, 'w') as out_tsv:
            count = 0
            for line in in_tsv:
                if count == 0:
                    out_tsv.write(line)
                    count = 1
                    continue
                row = line.strip().split('\t')
                sequences = defaultdict(str)
                read_positions = defaultdict(list)
                samples = defaultdict(int)
                print(row[0:3])
                for read_insert in row[5].split(','):
                    if read_insert == "NA":
                        continue
                    sample = read_insert.split(':')[0]
                    samples[sample] += 1
                    read = read_insert.split(':')[1]
                    orientation = read_insert.split(':')[2]
                    insert_start = int(read_insert.split(':')[3].split('-')[0])
                    insert_end = int(read_insert.split(':')[3].split('-')[1])
                    read_positions[read] = [insert_start, insert_end, orientation]
                # Get sequence from reads
                for fastq in fastqs:
                    for read in read_positions:
                        try:
                            seq = fastq.fetch(read)
                            if read_positions[read][2] == '-':
                                seq = reverse_complement(seq)
                            sequences[read] = seq
                        except KeyError:
                            # Read not found in this fastq
                            pass
                ret = []
                # Control sample if requested
                ret.append(filter_control(samples, control_sample))
                print("control")
                # Centromere/telomere
                ret.append(filter_centromere_telomere(row, centromere_positions, telomere_positions))
                print("filter_centromere_telomere")
                # Insert maps to same location as other mapping for supporting reads
                ret.append(filter_insert_map(row[0], int(row[1]), int(row[2]), read_positions, sequences, bam_read_positions, ref_aligner, min_mapq))
                print("filter_insert_map")
                # Low mapping quality at area of insertion
                ret.append(filter_mapping_qual(row[0], int(row[1]), int(row[2]), bam_ref_positions, min_mapq, chrom_list))
                print("filter_mapping_qual")
                # Poly A/T Tail
                ret.append(filter_poly_AT(read_positions, sequences))
                print("filter_poly_AT")
                # Target Site Duplications
                ret.append(filter_TSD(read_positions, sequences))
                print("filter_TSD")
                # Update the row's annotation based on which filters it failed
                row[6] = update_annotation(row[6], ret)
                out_tsv.write("\t".join(row)+'\n')

