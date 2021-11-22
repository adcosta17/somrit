import pysam
import argparse
import sys
import re
from collections import defaultdict
from .merge import merge_records

def calculate_sdust_score(seq):
    if seq == "":
        return 0
    triplet_counts = defaultdict(int)
    for i in range(0, len(seq) - 2):
        triplet_counts[seq[i:i+3]] += 1
    sum_score = 0
    for triplet in triplet_counts:
        count = triplet_counts[triplet]
        s = float(count * (count - 1)) / 2.0
        sum_score += s
    if len(seq) - 1 == 0:
        return 0
    sum_score /= (len(seq) - 1)
    return sum_score


def get_deletion_pos(cigarstring):
    # Count up the position on the read until we get to a deletion
    deletion_positions = []
    count = 0
    for cg in re.findall('[0-9]*[A-Z]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            count += int(cg[:cg.find("H")])
        elif cg.endswith('D'):
            if int(cg[:cg.find("D")]) >= 100:
                deletion_positions.append((count,int(cg[:cg.find("D")])))
    return deletion_positions

def update_annotation(annotation, update):
    if annotation == "PASS":
        annotation = update
    else:
        annotation = annotation+","+update
    return annotation


def get_tsv_record(record, max_ref_gap_at_candidate, min_detected_inclusion_length, min_mapq, min_insertion_length, min_flank_size, merged):
    records_to_output =[]
    read_annotation = "PASS"
    if record.mapq < min_mapq:
        read_annotation = "mapq<20" 
    # check for any long insertions
    aligned_pairs = record.get_aligned_pairs(matches_only=True)
    for idx in range(0, len(aligned_pairs) - 1):
        read_gap = aligned_pairs[idx + 1][0] - aligned_pairs[idx][0]
        ref_gap = aligned_pairs[idx + 1][1] - aligned_pairs[idx][1]
        if read_gap >= min_detected_inclusion_length and ref_gap <= max_ref_gap_at_candidate:
            ref_start = aligned_pairs[idx][1]
            ref_end = aligned_pairs[idx+1][1]
            read_start = aligned_pairs[idx][0]
            read_end = aligned_pairs[idx+1][0]
            insertion_sequence = ""
            if record.query_sequence is not None:
                insertion_sequence = record.query_sequence[read_start:read_end]
            is_merged = "0"
            if merged:
                is_merged = "1"
            annotation = read_annotation
            if not ((abs(ref_start - record.reference_start) > int(min_flank_size)) and 
                (abs(record.reference_end - ref_end) > int(min_flank_size))):
                annotation = update_annotation(annotation, "flank_size")
            if not read_gap >= min_insertion_length:
                annotation = update_annotation(annotation, "min_insertion_length")
            # convert read start and end positions if the alignment is rc
            if record.is_reverse and insertion_sequence != "":
                tmp = len(record.query_sequence) - read_start
                read_start = len(record.query_sequence) - read_end
                read_end = tmp
            deletion_positions = get_deletion_pos(record.cigarstring)
            if len(deletion_positions) > 0:
                for del_pos in deletion_positions:
                    if record.query_sequence is not None:
                        read_start_tmp = read_start
                        read_end_tmp = read_end
                        if record.is_reverse:
                            tmp = len(record.query_sequence) - read_start
                            read_start_tmp = len(record.query_sequence) - read_end
                            read_end_tmp = tmp
                        if ((abs(del_pos[0] - read_start_tmp) < 2*del_pos[1] or abs(del_pos[0] - read_end_tmp) < 2*del_pos[1]) and 
                            (read_end_tmp - read_start_tmp)*0.75 < del_pos[1]):
                            # update the annotation to indicate a possible nearby deletion => mapping artifact
                            annotation = update_annotation(annotation, "deletion_possible_mapping_artifact")
            if insertion_sequence == "":
                n = read_end - read_start
                insertion_sequence = ''.join([c*n for c in "N"])
            records_to_output.append(["%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s" % (record.reference_name, ref_start, ref_end, record.query_name, read_start, read_end, is_merged, insertion_sequence), annotation])
    return records_to_output


def extract_candidate_insertions(bam, output_tsv, min_insertion_length, min_detected_inclusion_length, min_flank_size, min_mapq, reference_gap_minimum, fastq_file):
    with open(output_tsv, 'w') as out_tsv:
        out_tsv.write("\t".join(["chromosome", "reference_insertion_start", "reference_insertion_end", "read_name", "read_insertion_start", "read_insertion_end", "merged", "insertion_sequence","pass_fail"])+"\n")
        header = pysam.AlignmentFile(bam).header
        for sq in header['SQ']:
            #print(sq['SN'])
            sam_reader = pysam.AlignmentFile(bam)
            tmp_sam_reader = sam_reader.fetch(contig=sq['SN'])
            read_to_reference_alignments = defaultdict(list)
            for record in tmp_sam_reader:
                read_to_reference_alignments[record.query_name].append(record)
            for read_id in read_to_reference_alignments:
                records = sorted(read_to_reference_alignments[read_id], key = lambda x: (x.reference_id, x.reference_start))
                # Go read by read based on sorted records. Merge first and then extract
                new_records = merge_records(records, header, reference_gap_minimum, min_mapq, min_detected_inclusion_length, fastq_file)
                for record in new_records:
                    for item in get_tsv_record(record[0], reference_gap_minimum, min_detected_inclusion_length, min_mapq, min_insertion_length, min_flank_size, record[1]):
                        out_tsv.write(item[0]+"\t"+item[1]+"\n")

