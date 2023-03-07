import pysam
import argparse
import sys
import re
from collections import defaultdict, namedtuple
import gzip
from os import listdir
from os.path import isfile, join
from multiprocessing import Lock, Queue
import threading
import time
import queue

Record = namedtuple('Record', ['query_name', 'orientation', 'ref_name', 'ref_start', 'ref_end', 'mapq', 'cigarstring', 'read_start', 'read_end', 'merged'])

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

def strip_soft_hard(cigarstring):
    read_start = 0
    read_end= 0
    new_cig = ""
    count = 0
    seen_match = False
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
            seen_match = True
            new_cig += cg
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
            seen_match = True
            new_cig += cg
        elif cg.endswith('='):
            count += int(cg[:cg.find("=")])
            seen_match = True
            new_cig += cg
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
            new_cig += cg
        elif cg.endswith('D'):
            new_cig += cg
        elif cg.endswith('P'):
            new_cig += cg
        elif cg.endswith('N'):
            new_cig += cg
        elif cg.endswith('S'):
            if not seen_match:
                read_start += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            if not seen_match:
                read_start += int(cg[:cg.find("H")])
    read_end = read_start+count
    return read_start, read_end, new_cig

def get_record_info(record):
    query_name = record.query_name
    orientation = '+'
    if record.is_reverse:
        orientation = '-'
    ref_name = record.reference_name
    ref_start = record.reference_start
    ref_end = record.reference_end
    mapq = record.mapping_quality
    read_start, read_end, cigarstring = strip_soft_hard(record.cigarstring)
    return Record(query_name, orientation, ref_name, ref_start, ref_end, mapq, cigarstring, read_start, read_end, False)

def merge_records(records, header, reference_gap_minimum, minimum_mapping_qual, min_detected_inclusion_length, f_lock):
    # Assume that records has been sorted by chromsosme and position
    records_to_return = {}
    alignments_to_output = []
    alignment_start_end = []
    need_hard_clip = []
    # Get the read sequence
    count = 0
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
                records_to_return[count] = get_record_info(record_1)
                count += 1
                record_1 = record_2
                i += 1
                continue
            # Need to now combine the two records together
            query_name = record_1.query_name
            orientation = '+'
            if record_1.is_reverse:
                orientation = '-'
            ref_name = record_1.reference_name
            ref_start = min(record_1.reference_start, record_2.reference_start)
            ref_end = max(record_1.reference_end, record_2.reference_end)
            mapq = record_1.mapping_quality
            tmp_str = ""
            cigarstring = ""
            if record_1.reference_start > record_2.reference_start:
                tmp_str = get_start_end_string(record_2.cigarstring, record_1.cigarstring, 
                    record_1.query_name, record_2.reference_start, record_2.reference_end, 
                    record_1.reference_start, record_1.reference_end)
                cigarstring = tmp_str.split("_")[2]
            else:
                tmp_str = get_start_end_string(record_1.cigarstring, record_2.cigarstring, 
                    record_1.query_name, record_1.reference_start, record_1.reference_end, 
                    record_2.reference_start, record_2.reference_end)
                cigarstring = tmp_str.split("_")[2]
            if tmp_str != "0_0_0":
                read_start = int(tmp_str.split("_")[0])
                read_end = int(tmp_str.split("_")[1])
                # Manually add details to records_to_return
                records_to_return[count] = Record(query_name, orientation, ref_name, ref_start, ref_end, mapq, cigarstring, read_start, read_end, True)
                count += 1
                records = records[i+1:]
                i = 1
            else:
                # Need to check insertions in each of the records together
                # Don't want to exclude them just because they aren't properly split
                records_to_return[count] = get_record_info(record_1)
                count += 1
                records_to_return[count] = get_record_info(record_2)
                count += 1
                i += 1
                if i < len(records):
                    record_1_merged = False
                    record_1 = records[i]
                    i += 1
        #Check insertion for record_1 now
        records_to_return[count] = get_record_info(record_1)
        count += 1
    else:
        # Single mapping for read, Check for insert
        records_to_return[count] = get_record_info(records[0])
        count += 1
    return records_to_return


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

def get_insertion_pos(cigarstring, min_detected_inclusion_length):
    # Count up the position on the read until we get to a deletion
    insert_positions = []
    read_count = 0
    ref_count = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
            ref_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
            ref_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            if int(cg[:cg.find("I")]) > min_detected_inclusion_length:
                insert_positions.append([read_count, int(cg[:cg.find("I")]), ref_count])
            read_count += int(cg[:cg.find("I")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('N'):
            ref_count += int(cg[:cg.find("N")])
    return insert_positions

def get_tsv_record(record, max_ref_gap_at_candidate, min_detected_inclusion_length, min_mapq, min_insertion_length, min_flank_size, min_read_len, in_fq, f_lock):
    records_to_output = {}
    read_annotation = "PASS"
    if record.mapq < min_mapq:
        read_annotation = "mapq<20" 
    # check for any long insertions
    insert_positions = get_insertion_pos(record.cigarstring, min_detected_inclusion_length)
    deletion_positions = get_deletion_pos(record.cigarstring)
    count = 0
    merged = False
    read_seq = None
    if len(insert_positions) > 0 :
        f_lock.acquire()
        try:
            read_seq = in_fq.fetch(record.query_name)
            if len(read_seq) < min_read_len:
                read_annotation = update_annotation(read_annotation, "min_read_length")
        except:
            pass
        f_lock.release()
    for insert in insert_positions:
        ref_start = insert[2]+record.ref_start
        ref_end = insert[2]+1+record.ref_start
        read_start = insert[0]
        read_end = insert[0]+insert[1]
        if record.orientation == '-' and insertion_sequence != "":
            tmp = len(read_seq) - read_start
            read_start = len(read_seq) - read_end
            read_end = tmp
        insertion_sequence = ""
        if read_seq is not None:
            insertion_sequence = read_seq[read_start:read_end]
        if record.merged:
            merged = True
        annotation = read_annotation
        if not ((abs(ref_start - record.ref_start) > int(min_flank_size)) and 
            (abs(record.ref_end - ref_end) > int(min_flank_size))):
            annotation = update_annotation(annotation, "flank_size")
        if not read_end-read_start >= min_insertion_length:
            annotation = update_annotation(annotation, "min_insertion_length")
        if len(deletion_positions) > 0:
            for del_pos in deletion_positions:
                if read_seq is not None:
                    read_start_tmp = read_start
                    read_end_tmp = read_end
                    if record.orientation == '-':
                        tmp = len(read_seq) - read_start
                        read_start_tmp = len(read_seq) - read_end
                        read_end_tmp = tmp
                    if ((abs(del_pos[0] - read_start_tmp) < 2*del_pos[1] or abs(del_pos[0] - read_end_tmp) < 2*del_pos[1]) and 
                        (read_end_tmp - read_start_tmp)*0.75 < del_pos[1]):
                        # update the annotation to indicate a possible nearby deletion => mapping artifact
                        annotation = update_annotation(annotation, "deletion_possible_mapping_artifact")
        if insertion_sequence == "":
            n = read_end - read_start
            insertion_sequence = "N"*n
        output = "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s" % (record.ref_name, ref_start, ref_end, record.query_name, read_start, read_end, record.orientation, insertion_sequence,annotation)
        records_to_output[count] = output
        count += 1
    return records_to_output, merged


class myThread (threading.Thread):
   def __init__(self, i, q,header, reference_gap_minimum, min_mapq, min_detected_inclusion_length, in_fq, min_insertion_length, min_flank_size, min_read_len, t_lock, f_lock, out_tsv, out_merged):
      threading.Thread.__init__(self)
      self.i = i
      self.q = q
      self.header = header
      self.reference_gap_minimum = reference_gap_minimum
      self.min_mapq = min_mapq
      self.min_detected_inclusion_length = min_detected_inclusion_length
      self.in_fq = in_fq
      self.min_insertion_length = min_insertion_length
      self.min_flank_size = min_flank_size
      self.min_read_len = min_read_len
      self.t_lock = t_lock
      self.f_lock = f_lock
      self.out_tsv = out_tsv
      self.out_merged = out_merged
   def run(self):
        count = 0
        merge_time = 0
        record_time = 0
        out_time = 0
        while(True):
            if self.q.empty():
                return
            data = self.q.get()
            read_id = data[0]
            record_list = data[1]
            #if count % 1000 == 0:
            #    print(str(self.i)+"\t"+str(count))
            #count += 1
            if len(record_list) > 1:
                records = sorted(record_list, key = lambda x: (x.reference_id, x.reference_start))
                # Go read by read based on sorted records. Merge first and then extract
            else:
                records = record_list
            new_records = merge_records(records, self.header, self.reference_gap_minimum, self.min_mapq, self.min_detected_inclusion_length, self.f_lock)
            for c in new_records:
                inserts, merged = get_tsv_record(new_records[c], self.reference_gap_minimum, self.min_detected_inclusion_length, self.min_mapq, self.min_insertion_length, self.min_flank_size, self.min_read_len, self.in_fq, self.f_lock)
                self.t_lock.acquire()
                for i in inserts:
                    self.out_tsv.write(inserts[i]+"\n")
                if merged:
                    self.out_merged.write(read_id+"\n")
                self.t_lock.release()


def extract_candidate_insertions(bam, output_tsv, output_merged, min_insertion_length, min_detected_inclusion_length, min_flank_size, min_read_len, min_mapq, reference_gap_minimum, fastq_file, threads):
    t_lock = threading.Lock()
    f_lock = threading.Lock()
    with open(output_tsv, 'w') as out_tsv, open(output_merged, 'w') as out_merged, pysam.FastaFile(filename = fastq_file) as in_fq:
        out_tsv.write("\t".join(["chromosome", "reference_insertion_start", "reference_insertion_end", "read_name", "read_insertion_start", "read_insertion_end", "orientation", "insertion_sequence","pass_fail"])+"\n")
        header = pysam.AlignmentFile(bam).header
        sam_reader = pysam.AlignmentFile(bam)
        for sq in header['SQ']:
            #print(sq['SN'])
            read_to_reference_alignments = defaultdict(list)
            for record in sam_reader.fetch(contig=sq['SN']):
                if record.mapping_quality < min_mapq:
                    continue
                read_to_reference_alignments[record.query_name].append(record)
            count = 0
            q_list = [None]*threads
            for i in range(threads):
                q_list[i] = queue.Queue()
            for read_id in read_to_reference_alignments:
                q_list[count%threads].put([read_id, read_to_reference_alignments[read_id]])
                count += 1
            thread_list = [None] *threads
            for i in range(threads):
                thread_list[i] = myThread(i, q_list[i], header, reference_gap_minimum, min_mapq, min_detected_inclusion_length, in_fq, min_insertion_length, min_flank_size, min_read_len, t_lock,f_lock,out_tsv,out_merged)
                thread_list[i].start()
            for i in range(threads):
                thread_list[i].join()

                
                


