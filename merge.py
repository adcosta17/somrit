#---------------------------------------------------------
# Copyright 2023 Ontario Institute for Cancer Research
# Written by Alister D'Costa (alister.d'costa@oicr.on.ca)
#---------------------------------------------------------
# merge.py


import argparse
import pysam
import sys
import re
import gzip
from collections import defaultdict
from os import listdir
from os.path import isfile, join

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def is_number(n):
    try:
        float(n)
    except ValueError:
        return False
    return True

def to_sam(query_name, strand, target_name, position, qual, cigar):
    if(strand == '+'):
        bits = 0
    else:
        bits = 16
    print(query_name+"\t"+str(bits)+"\t"+target_name+"\t"+str(position)+"\t"+str(qual)+"\t"+cigar+"\t*\t0\t0\t*\t*")

# Removes soft and hard clips off the end of a cigar string
def clean_cigar(cigar):
    seperator = ''
    cgs = re.findall('[0-9]*[A-Z]', cigar)
    start = False
    end = False
    if cgs[0].endswith('H') or cgs[0].endswith('S'):
        start = True
    if cgs[len(cgs)-1].endswith('H') or cgs[len(cgs)-1].endswith('S'):
        end = True
    if start and end:
        return seperator.join(cgs[1:len(cgs)-1])
    elif start:
        return seperator.join(cgs[1:])
    elif end:
        return seperator.join(cgs[:len(cgs)-1])
    else:
        return seperator.join(cgs)

# Checks for an insertion larger than the given distance in a bam record
def check_insertion_bam(record, distance):
    #Look for an insertion in the read
    #Parse cigar string looking for distance sized or more insertion
    if record.mapping_quality < args.minimum_mapping_qual:
        return False
    for cg in re.findall('[0-9]*[A-Z]', str(record.cigarstring)):
        if cg.endswith('I'):
            if int(cg[:cg.find("I")]) >= int(distance):
                return True
                break;
    return False

# Checks for an insertion or hard clip larger than the given distance in a bam record
# Needed as hard clips don't appear in the seq in a record. Need to get them from the fastq
def check_insertion_hard_clip_bam(record, distance):
    #Look for an insertion in the read
    #Parse cigar string looking for distance sized or more insertion
    #if record.mapping_quality < args.minimum_mapping_qual:
    #    return False
    for cg in re.findall('[0-9]*[A-Z]', str(record.cigarstring)):
        if cg.endswith('H'):
            if int(cg[:cg.find("H")]) > 0:
                return True
                break;
    return False

# Checks for an insertion, hardclip or softclip larger than the given distance in a bam record
def check_insertion_or_soft_clip_bam(record, distance):
    #Look for an insertion in the read
    #Parse cigar string looking for distance sized or more insertion
    #if record.mapping_quality < args.minimum_mapping_qual:
    #    return False
    for cg in re.findall('[0-9]*[A-Z]', str(record.cigarstring)):
        if cg.endswith('I'):
            if int(cg[:cg.find("I")]) >= int(distance):
                return True
        if cg.endswith('S'):
            if int(cg[:cg.find("S")]) >= int(distance):
                return True
        if cg.endswith('H'):
            if int(cg[:cg.find("H")]) >= int(distance):
                return True
    return False

# Gets the total read length based on the cigar string
def get_read_length(cigarstring):
    # Gets the read length based on the Cigar String
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
    return count

def get_id(header, name):
    i = 0
    while i < len(header['SQ']):
        if header['SQ'][i]['SN'] == name:
            return i
        else:
            i += 1
    return 0

def get_start_end_string(cigar1, cigar2, read_id, ref1_start, ref1_end, ref2_start, ref2_end):
    # Assume the cigars are oriented such the alignment is cigar1 - indel - cigar2
    # Remove any hard or soft clipped bases from the start of the read. This is our starting position
    # Next iterate over a cleaned cigar string to get how many bases in first half on read, add the indel size and then iterate over second cleaned cigar string
    # now we know length of read to traverse. Gen ends pos by adding length to start pos
    cg1_start = 0
    cgs = re.findall('[0-9]*[A-Z]', cigar1)
    if cgs[0].endswith('S'):
        cg1_start += int(cgs[0][:cgs[0].find("S")])
    elif cgs[0].endswith('H'):
        cg1_start += int(cgs[0][:cgs[0].find("H")])
    cg1_cleaned = clean_cigar(cigar1)
    cg1_len = get_read_length(cg1_cleaned)
    cg1_end = cg1_start + cg1_len
    cg2_start = 0
    cgs = re.findall('[0-9]*[A-Z]', cigar2)
    if cgs[0].endswith('S'):
        cg2_start += int(cgs[0][:cgs[0].find("S")])
    elif cgs[0].endswith('H'):
        cg2_start += int(cgs[0][:cgs[0].find("H")])
    cg2_cleaned = clean_cigar(cigar2)
    cg2_len = get_read_length(cg2_cleaned)
    cg2_end = cg2_start + cg2_len
    if (cg1_start >= cg2_end) or ((cg2_start - cg1_end) < 0):
        return "0_0_0"
    if ref2_start - ref1_end == 0:
        # no issues
        return str(cg1_start) + "_" + str(cg2_end) + "_" + cg1_cleaned + str(cg2_start - cg1_end) + "I" + cg2_cleaned
    elif ref2_start - ref1_end < 0:
        n_count1 = abs(ref2_start - ref1_end)/2
        n_count2 = n_count1
        if abs(ref2_start - ref1_end) % 2 == 1:
            n_count2 += 1
        return str(cg1_start) + "_" + str(cg2_end) + "_" + cg1_cleaned + str(n_count2) + "N" + str((cg2_start - cg1_end)) + "I" + str(n_count1) + "N" + cg2_cleaned
    else:
        # Need to mofidy end and insert
        n_count1 = abs(ref2_start - ref1_end)/2
        n_count2 = n_count1
        if abs(ref2_start - ref1_end) % 2 == 1:
            n_count2 += 1
        return str(cg1_start) + "_" + str(cg2_end) + "_" + cg1_cleaned + str(n_count2) + "N" + str((cg2_start - cg1_end)) + "I" + str(n_count1) + "N" + cg2_cleaned


# Gets the start and end postions of any any hard and soft clips on a read
# if both at one end combine them together into one segment
def get_read_pos(record, min_insert, front):
    cigars = re.findall('[0-9]*[A-Z]', record.cigarstring)
    if not front:
        cigars = cigars[::-1]
    # Search until we find a softclip. Will happen at start or end. May have a hard clip before
    start = 0
    end = 0
    hard = False
    for cg in cigars:
        if cg.endswith('H'):
            step = int(cg[:cg.find("H")])
            end += step
            hard = True
        elif cg.endswith('S'):
            step = int(cg[:cg.find("S")])
            end += step
        else:
            break
    if end != 0 and hard:
        # Things are in the correct orientation
        end += start
        if abs(end - start) >= int(min_insert):
            return(start, end)
        else:
            return None
    else:
        return None

def merge_records(records, header, reference_gap_minimum, minimum_mapping_qual, min_detected_inclusion_length, fastq_file):
    # Assume that records has been sorted by chromsosme and position
    records_to_return = []
    alignments_to_output = []
    alignment_start_end = []
    need_hard_clip = []
    if len(records) > 1:
        # Check for split reads
        record_1 = records[0]
        i = 1
        record_1_merged = False
        while i < len(records):
            record_2 = records[i]
            if (record_1.query_name != record_2.query_name or 
                record_1.reference_name != record_2.reference_name or
                record_1.is_reverse != record_2.is_reverse or
                (record_1.reference_length + record_2.reference_length) > 1.2* get_read_length(record_1.cigarstring) or
                (abs(record_1.reference_end - record_2.reference_start) > int(reference_gap_minimum) and 
                abs(record_2.reference_end - record_1.reference_start) > int(reference_gap_minimum)) or
                record_1.mapping_quality < minimum_mapping_qual or record_2.mapping_quality < minimum_mapping_qual):
                if check_insertion_or_soft_clip_bam(record_1, min_detected_inclusion_length):
                    if check_insertion_hard_clip_bam(record_1, min_detected_inclusion_length):
                        need_hard_clip.append(record_1)
                    else:
                        records_to_return.append([record_1, False])
                record_1 = record_2
                i += 1
                continue
            # Need to now combine the two records together
            a = pysam.AlignedSegment(header=header)
            a.query_name = record_1.query_name
            a.flag = min(record_1.flag, record_2.flag)
            a.reference_id = get_id(header, record_1.reference_name)
            a.reference_start = min(record_1.reference_start, record_2.reference_start)
            a.mapping_quality = record_1.mapping_quality
            insert_size = 0
            tmp_str = ""
            if record_1.reference_start > record_2.reference_start:
                tmp_str = get_start_end_string(record_2.cigarstring, record_1.cigarstring, 
                    record_1.query_name, record_2.reference_start, record_2.reference_end, 
                    record_1.reference_start, record_1.reference_end)
                a.cigarstring = tmp_str.split("_")[2]
            else:
                tmp_str = get_start_end_string(record_1.cigarstring, record_2.cigarstring, 
                    record_1.query_name, record_1.reference_start, record_1.reference_end, 
                    record_2.reference_start, record_2.reference_end)
                a.cigarstring = tmp_str.split("_")[2]
            if tmp_str != "0_0_0":
                alignment_start_end.append(tmp_str)
                alignments_to_output.append(a)
                tmp = [a]
                tmp.extend(records[i+1:])
                records = tmp
                print(a)
                i = 1
            else:
                # Need to check insertions in each of the records together
                # Don't want to exclude them just because they aren't properly split
                if check_insertion_or_soft_clip_bam(record_1, min_detected_inclusion_length):
                    if check_insertion_hard_clip_bam(record_1, min_detected_inclusion_length):
                        need_hard_clip.append(record_1)
                    else:
                        records_to_return.append([record_1, False])
                if check_insertion_or_soft_clip_bam(record_2, min_detected_inclusion_length):
                    if check_insertion_hard_clip_bam(record_2, min_detected_inclusion_length):
                        need_hard_clip.append(record_2)
                    else:
                        records_to_return.append([record_2, False])
                i += 1
                if i < len(records):
                    record_1_merged = False
                    record_1 = records[i]
                    i += 1
        #Check insertion for record_1 now
        if check_insertion_or_soft_clip_bam(record_1, min_detected_inclusion_length) and not record_1_merged:
            if check_insertion_hard_clip_bam(record_1, min_detected_inclusion_length):
                need_hard_clip.append(record_1)
            else:
                records_to_return.append([record_1, False])
    else:
        # Single mapping for read, Check for insert
        if check_insertion_or_soft_clip_bam(records[0], min_detected_inclusion_length):
            if check_insertion_hard_clip_bam(records[0], min_detected_inclusion_length):
                need_hard_clip.append(records[0])
            else:
                records_to_return.append([records[0], False])
    with pysam.FastaFile(filename = fastq_file) as fq:
        i = 0
        seen = {}
        for i in range(len(alignments_to_output)):
            if alignments_to_output[i].query_name not in seen:
                seen[alignments_to_output[i].query_name] = []
            try:
                start = int(alignment_start_end[i].split("_")[0])
                end = int(alignment_start_end[i].split("_")[1])
                seq = fq.fetch(alignments_to_output[i].query_name)
                # Check if the SEQ needs to be reverse complimented based on the bitflag
                if end - start != get_read_length(alignments_to_output[i].cigarstring):
                   #print(alignments_to_output[i].query_name + " " + alignments_to_output[i].reference_name + " " + 
                   # str(len(seq)) + " " + str(end - start) + " " + str(get_read_length(alignments_to_output[i].cigarstring)))
                   continue 
                if(alignments_to_output[i].is_reverse):
                    alignments_to_output[i].query_sequence = reverse_complement(seq)[start:end]
                else:
                    alignments_to_output[i].query_sequence = seq[start:end]
                seen[alignments_to_output[i].query_name].append(alignments_to_output[i])
                out_count += 1
            except KeyError:
                pass
        for read in seen:
            i = 0
            while i < len(seen[read]):
                if i + 1 < len(seen[read]):
                    if seen[read][i].reference_start == seen[read][i+1].reference_start or seen[read][i].reference_end == seen[read][i+1].reference_end:
                        pass
                    else:
                        records_to_return.append([alignments_to_output[i], True])
                else:
                    records_to_return.append([alignments_to_output[i], True])
                i += 1
        seen = []
        for i in range(len(need_hard_clip)):
            if i in seen:
                continue
            try:
                record = need_hard_clip[i]
                read_seq = fq.fetch(record.query_name)
                # Look at the front of the sequence for a soft/hard clip
                read_pos = get_read_pos(record, min_detected_inclusion_length, True)
                if read_pos != None:
                    #get sequence of clip and add it accordingly
                    #Have to remove any soft clip off the front and replace it with the soft and hardclip sequence. 
                    #Adjust the cigar string to now represent a softclip of the correct size
                    cigars = re.findall('[0-9]*[A-Z]', record.cigarstring)
                    adjust = False
                    soft = 0
                    if cigars[0].endswith('H') and cigars[1].endswith('S'):
                        # have to adjust cigar string
                        soft = int(cigars[1][:cigars[1].find("S")])
                        cigars = cigars[1:]
                        cigars[0] = str(read_pos[1])+"S"
                        adjust = True
                    elif cigars[0].endswith('H'):
                        cigars[0] = str(read_pos[1])+"S"
                        adjust = True
                    if adjust:
                        record.cigarstring = "".join(cigars)
                        # remove soft bases from seq
                        record.query_sequence = record.query_sequence[soft:]
                        # now add region from read pos in
                        record.query_sequence = read_seq[read_pos[0]:read_pos[1]] + record.query_sequence
                # Look at the end of the sequence for a soft/hard clip
                read_pos = get_read_pos(record, min_detected_inclusion_length, False)
                if read_pos != None:
                    #get sequence of clip and add it accordingly
                    #Have to remove any soft clip off the end and replace it with the soft and hardclip sequence. 
                    #Adjust the cigar string to now represent a softclip of the correct size
                    cigars = re.findall('[0-9]*[A-Z]', record.cigarstring)[::-1]
                    adjust = False
                    soft = 0
                    if cigars[0].endswith('H') and cigars[1].endswith('S'):
                        # have to adjust cigar string
                        soft = int(cigars[1][:cigars[1].find("S")])
                        cigars = cigars[1:]
                        cigars[0] = str(read_pos[1]-read_pos[0])+"S"
                        adjust = True
                    elif cigars[0].endswith('H'):
                        cigars[0] = str(read_pos[1]-read_pos[0])+"S"
                        adjust = True
                    if adjust:
                        record.cigarstring = "".join(cigars[::-1])
                        # remove soft bases from seq
                        record.query_sequence = record.query_sequence[:record.query_length - soft]
                        # now add region from read pos in
                        record.query_sequence = record.query_sequence + read_seq[read_pos[0]:read_pos[1]]
                    # If its just a soft clip we don't have to do anything.
                seen.append(i)
                records_to_return.append([record, False])
                out_count += 1
            except KeyError:
                pass
    return records_to_return
