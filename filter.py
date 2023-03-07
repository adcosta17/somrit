import argparse
import sys
import re
from collections import defaultdict
import parasail
import mappy as mp
from intervaltree import Interval, IntervalTree
import pysam
from multiprocessing import Lock, Queue
import threading
import time
import queue

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
user_matrix = parasail.matrix_create("ACGT", 2, -1)

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def get_nearby(chrom, start, end, seen):
    i = start
    count = 0
    nearby = {}
    if chrom not in seen:
        return nearby
    while i < end:
        if i in seen[chrom]:
            nearby[count] = seen[chrom][i]
            count += 1
        i += 1
    return nearby

def get_nearby_with_position(chrom, start, end, seen):
    i = start
    nearby = {}
    count = 0
    if chrom not in seen:
        return nearby
    while i < end:
        if i in seen[chrom]:
            nearby[count] = [seen[chrom][i], i]
            count += 1
        i += 1
    return nearby

def filter_control(chrom, start, end, seen, control_sample):
    if control_sample is not None:
        nearby = get_nearby(chrom, start-500, end+500, seen)
        for item in nearby:
            for j in nearby[item]:
                sample = j.split(':')[0]
            if sample == control_sample:
                return "In_Control_Sample"
    return None

def get_polymorphic_in_window(chrom, start, end, seen, softclips, control_sample, bam_ref_positions, window, cutoff, total_read_count, nearby):
    insert_reads = {}
    for i in nearby:
        for j in nearby[i][0]:
            if abs(nearby[i][1]-start) < window or abs(nearby[i][1]-end) < window:
                sample = j.split(':')[0]
                if sample == control_sample:
                    continue
                read = j.split(':')[1]
                insert_reads[read] = 1
    insert_count = len(insert_reads)
    softclip_count = 0
    if chrom in softclips:
        softclip_count = len(softclips[chrom][str(start)+":"+str(end)])
    if total_read_count == 0:
        return "NA"
    if (insert_count+softclip_count)/total_read_count > cutoff:
        return "Polymorphic_"+str(window)+"_"+str((insert_count+softclip_count)/total_read_count)
    return "NA"

def check_polymorphic(chrom, start, end, seen, softclips, control_sample, bam_ref_positions):
    total_reads = {}
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
                    total_reads[item[0]] = 1
    total_read_count = len(total_reads)
    nearby = get_nearby_with_position(chrom, start-500, end+500, seen)
    # 500, 0.8
    ret = get_polymorphic_in_window(chrom, start, end, seen, softclips, control_sample, bam_ref_positions, 500, 0.8, total_read_count, nearby)
    if ret != "NA":
        return ret
    # 200, 0.5
    ret = get_polymorphic_in_window(chrom, start, end, seen, softclips, control_sample, bam_ref_positions, 200, 0.5, total_read_count, nearby)
    if ret != "NA":
        return ret
    # 100, 0.3
    ret = get_polymorphic_in_window(chrom, start, end, seen, softclips, control_sample, bam_ref_positions, 100, 0.3, total_read_count, nearby)
    if ret != "NA":
        return ret
    return "NA"

def filter_centromere_telomere(row, centromere_positions, telomere_positions):
    if centromere_positions[row[0]][int(row[1]):int(row[2])]:
        return "In_Centromere"
    if telomere_positions[row[0]][int(row[1]):int(row[2])]:
        return "In_Teloomere"
    return None

def filter_insert_map(chrom, start, end, read_positions, sequence, bam_read_positions, ref_aligner, min_mapq):
    # Then for each read that has a secondary mapping other than the insert pos, extract the insert sequence
    if len(bam_read_positions) == 0:
        return None
    read_failed = {}
    for hit in ref_aligner.map(sequence):
        for read in read_positions:
            bam_positions = bam_read_positions[read]
            # Map the insert sequence to the ref and check position. If position overlaps a secondary mapping, flag
            for bam_pos in bam_positions:
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

def filter_poly_AT(sequence):
    if len(sequence) < 50:
        return "No_Poly_A_Tail:Insert_to_small"
    start = sequence[:50]
    end = sequence[len(sequence)-50:]
    pA = "A"*10
    pT = "T"*10
    if pA in start or pT in start or pA in end or pT in end:
        return str("Has_Poly_A_Tail")
    return "No_Poly_A_Tail"

def check_match(cigar, sequence):
    cgs = re.findall('[0-9]*[A-Z=]', cigar)
    count = 0
    for cg in cgs:
        if cg.endswith('M'):
            if int(cg[:cg.find("M")]) >= 5:
                return sequence[count:count+int(cg[:cg.find("M")])]
            count += int(cg[:cg.find("M")])
        elif cg.endswith('='):
            if int(cg[:cg.find("=")]) >= 5:
                return sequence[count:count+int(cg[:cg.find("=")])]
            count += int(cg[:cg.find("=")])
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
    return None

def filter_TSD(chrom, start, end, sequence, ref_contig):
    upstream_ref = ref_contig[start-25:start]
    downstream_ref = ref_contig[end:end+25]
    # Align the first 25 bp of the insert sequence to the upstream and downstream ref
    insert_start_seq = sequence[:25]
    insert_end_seq = sequence[len(sequence)-25:]
    # Need to look for duplication between the ref_seqs and insert seq
    # Align the two sequences together and check for at least 5 matches in a row as a possible TSD
    result_start_up = parasail.sw_trace(insert_start_seq, upstream_ref, 11, 1, user_matrix)
    cigar = result_start_up.cigar.decode.decode('utf-8')
    tsd = check_match(cigar, insert_start_seq)
    if tsd is not None:
        return tsd
    result_start_down = parasail.sw_trace(insert_start_seq, downstream_ref, 11, 1, user_matrix)
    cigar = result_start_down.cigar.decode.decode('utf-8')
    tsd = check_match(cigar, insert_start_seq)
    if tsd is not None:
        return tsd
    result_end_up = parasail.sw_trace(insert_end_seq, upstream_ref, 11, 1, user_matrix)
    cigar = result_end_up.cigar.decode.decode('utf-8')
    tsd = check_match(cigar, insert_end_seq)
    if tsd is not None:
        return tsd
    result_end_down = parasail.sw_trace(insert_end_seq, downstream_ref, 11, 1, user_matrix)
    cigar = result_end_down.cigar.decode.decode('utf-8')
    tsd = check_match(cigar, insert_end_seq)
    if tsd is not None:
        return tsd
    return "No_TSD_Found"

def filter_min_reads(chrom, start, end, seen, min_reads):
    reads = {}
    nearby = get_nearby(chrom, start-50, end+50, seen)
    for item in nearby:
        for j in nearby[item]:
            read = j.split(':')[1]
            reads[read] = 1
    if len(reads) < min_reads:
        return "min_reads"
    return None


def get_motif(chrom, start, end, sequence, ref_contig, orientation, tsd):
    if tsd == "No_TSD_Found":
        return "NA"
    # Have a TSD, find it in the reference
    seq_to_search = ref_contig[start-25:start+25]
    i = seq_to_search.find(tsd)
    if orientation == '+':
        return reverse_complement(ref_contig[start-25+i-2:start-25+i+4])[::-1]
    else:
        i = i+len(tsd)
        return ref_contig[start-25+i-4:start-25+i+2]


def update_annotation(annotation, filters):
    ret = annotation
    new = []
    for f in filters:
        if f is not None:
            new.append(f)
    if len(new) > 0:
        if ret == "PASS":
            ret = ','.join(new)
        else:
            ret = annotation+','+','.join(new)
    return ret


class myThread2(threading.Thread):
   def __init__(self, i, q, seen, control_sample, centromere_positions, telomere_positions, bam_read_positions, ref_aligner, t_lock, min_mapq, bam_ref_positions, chrom_list, min_reads, softclips, ref_contigs, out_tsv):
      threading.Thread.__init__(self)
      self.i = i
      self.q = q
      self.seen = seen
      self.control_sample = control_sample
      self.centromere_positions = centromere_positions
      self.telomere_positions = telomere_positions
      self.bam_read_positions = bam_read_positions
      self.ref_aligner = ref_aligner
      self.min_mapq = min_mapq
      self.t_lock = t_lock
      self.bam_ref_positions = bam_ref_positions
      self.out_tsv = out_tsv
      self.chrom_list = chrom_list
      self.min_reads = min_reads
      self.softclips = softclips
      self.ref_contigs = ref_contigs
      self.string = ""
   def run(self):
        times = defaultdict(float)
        count = 0 
        while(True):
            if self.q.empty():
                self.t_lock.acquire()
                self.out_tsv.write(self.string)
                self.t_lock.release()
                return
            data = self.q.get()
            row = data[0]
            if row[6] != "PASS":
                # Only filter things that are possible insertions to save time
                row.extend(["NA","NA","NA","NA"])
                self.string += "\t".join(row)+'\n'
                count += 1
                continue
            sequence = row[4]
            read_positions = data[1]
            ret = []
            # Control sample if requested
            #to_print = False
            #if self.i == 0 and count % 1000 == 0:
            #    to_print = True
            #if to_print:
            #    print(row[0]+"\t"+row[1]+"\t"+row[2]+"\t"+row[3])
            #start = time.time()
            ret.append(filter_control(row[0], int(row[1]), int(row[2]), self.seen, self.control_sample))
            #end = time.time()
            #if to_print:
            #    times["filter_control"] += (end-start)
            #    print("filter_control " + str(times["filter_control"]))
            # Centromere/telomere
            #start = time.time()
            ret.append(filter_centromere_telomere(row, self.centromere_positions, self.telomere_positions))
            #end = time.time()
            #if to_print:
            #    times["Centromere/Telomere"] += (end-start)
            #    print("Centromere/Telomere " + str(times["Centromere/Telomere"]))
            # Insert maps to same location as other mapping for supporting reads
            #start = time.time()
            #ret.append(filter_insert_map(row[0], int(row[1]), int(row[2]), read_positions, sequence, self.bam_read_positions, self.ref_aligner, self.min_mapq))
            #end = time.time()
            #if to_print:
            #    times["InsertMap"] += (end-start)
            #    print("InsertMap "  + str(times["InsertMap"]))
            # Low mapping quality at area of insertion
            #start = time.time()
            ret.append(filter_mapping_qual(row[0], int(row[1]), int(row[2]), self.bam_ref_positions, self.min_mapq, self.chrom_list))
            #end = time.time()
            #if to_print:
            #    times["Mapping Qual"] += (end-start)
            #    print("Mapping Qual " + str(times["Mapping Qual"]))
            # Check for min_read support
            #start = time.time()
            ret.append(filter_min_reads(row[0], int(row[1]), int(row[2]), self.seen, self.min_reads))
            #end = time.time()
            #if to_print:
            #    times["MinReads"] += (end-start)
            #    print("MinReads " + str(times["MinReads"]))
            row[6] = update_annotation(row[6], ret)
            # Check to see if insertion is polymorphic or not
            #start = time.time()
            row.append(check_polymorphic(row[0], int(row[1]), int(row[2]), self.seen, self.softclips, self.control_sample, self.bam_ref_positions))
            #end = time.time()
            #if to_print:
            #    times["Polymorphic"] += (end-start)
            #    print("Polymorphic " + str(times["Polymorphic"]))
            # Poly A/T Tail
            #start = time.time()
            row.append(filter_poly_AT(sequence))
            #end = time.time()
            #if to_print:
            #    times["PolyA"] += (end-start)
            #    print("PolyA " + str(times["PolyA"]))
            # Target Site Duplications
            #start = time.time()
            tsd = filter_TSD(row[0], int(row[1]), int(row[2]), sequence, self.ref_contigs[row[0]])
            #end = time.time()
            #if to_print:
            #    times["TSD"] += (end-start)
            #    print("TSD " + str(times["TSD"]))
            row.append(tsd)
            # Update the row's annotation based on which filters it failed
            #start = time.time()
            row.append(get_motif(row[0], int(row[1]), int(row[2]), sequence, self.ref_contigs[row[0]], row[9], tsd))
            #end = time.time()
            #if to_print:
            #    times["Motif"] += (end-start)
            #    print("Motif" + str(times["Motif"]))
            self.string += "\t".join(row)+'\n'
            count += 1
            if count % 10000 == 0:
                self.t_lock.acquire()
                self.out_tsv.write(self.string)
                self.t_lock.release()
                self.string = ""


def filter_insertions(input_tsv, output_tsv, bam_list, fastq_list, contromeres, telomeres, control_sample, min_mapq, ref, cluster_window, chrs_to_use, min_reads, threads):
    # Open and read in centromere and telomere positions
    telomere_positions= defaultdict(IntervalTree)
    centromere_positions = defaultdict(IntervalTree)
    seen = {}
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
    reads_in_tsv = defaultdict(int)
    all_positions_list = defaultdict(list)
    with open(input_tsv, 'r') as in_tsv:
        count = 0
        for line in in_tsv:
            if count == 0:
                count = 1
                continue
            row = line.strip().split('\t')
            all_positions_list[row[0]].append((int(row[1])-cluster_window, int(row[2])+cluster_window))
            chrom = row[0]
            start = int(row[1])
            end =  int(row[2])
            for read_insert in row[5].split(','):
                if read_insert == "NA":
                    continue
                read = read_insert.split(':')[1]
                if chrom not in seen:
                    seen[chrom] = defaultdict(list)
                updated_insert = read_insert+":"+chrom+":"+row[1]+"-"+row[2]
                seen[chrom][start].append(updated_insert)
                reads_in_tsv[read] = 1
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
    softclips = {}
    count = 0
    for bam in bam_list.split(','):
        reader = pysam.AlignmentFile(bam)
        for record in reader.fetch():
            count += 1
            #if count % 10000 == 0:
            #    print(count)
            cg_tuples = record.cigartuples
            if cg_tuples[0][0] == 4 or cg_tuples[0][0] == 5:
                # Record starts with a hard or softclip
                if cg_tuples[0][1] >= 100:
                    # Have a clip that is larger than 100bp, check to see if it intersects a position we care about
                    nearby = get_nearby(record.reference_name, record.reference_start-50, record.reference_start+50, seen)
                    if len(nearby) > 0:
                        # Have a nearby softclip, add it to the list for the position
                        if record.reference_name not in softclips:
                            softclips[record.reference_name] = defaultdict(set)
                        for item in nearby:
                            for j in nearby[item]:
                                pos = j.split(':')[5]
                                softclips[record.reference_name][pos].add(record.query_name)
            elif cg_tuples[len(cg_tuples)-1][0] == 4 or cg_tuples[len(cg_tuples)-1][0] == 5:
                # Record starts with a hard or softclip
                if cg_tuples[len(cg_tuples)-1][1] >= 100:
                    # Have a clip that is larger than 100bp, check to see if it intersects a position we care about
                    nearby = get_nearby(record.reference_name, record.reference_end-50, record.reference_end+50, seen)
                    if len(nearby) > 0:
                        # Have a nearby softclip, add it to the list for the position
                        if record.reference_name not in softclips:
                            softclips[record.reference_name] = defaultdict(set)
                        for item in nearby:
                            for j in nearby[item]:
                                pos = j.split(':')[5]
                                softclips[record.reference_name][pos].add(record.query_name)
            nearby = all_positions[record.reference_name][record.reference_start:record.reference_end]
            if len(nearby) > 0:
                if record.mapping_quality >= min_mapq:
                    bam_read_positions[record.query_name].append([record.reference_name, record.reference_start, record.reference_end])
                # Have a record in a region we care about
                if record.reference_name not in bam_ref_positions:
                    bam_ref_positions[record.reference_name] = defaultdict(list)
                for item in nearby:
                    #print(record.query_name+" "+record.reference_name+":"+str(item.begin)+"-"+str(item.end))
                    bam_ref_positions[record.reference_name][str(item.begin)+"-"+str(item.end)].append([record.query_name, record.mapping_quality, record.reference_start, record.reference_end])
    print("Setup Dicts")
    ref_aligner = mp.Aligner(ref)
    ref_file = pysam.FastaFile(filename=ref)
    ref_contigs = {}
    # Read in tsv and then run each filter
    t_lock = threading.Lock()
    q_list = [None]*threads
    for i in range(threads):
        q_list[i] = queue.Queue()
    with open(input_tsv, 'r') as in_tsv:
        with open(output_tsv, 'w') as out_tsv:
            count = 0
            for line in in_tsv:
                if count == 0:
                    out_tsv.write(line.strip())
                    out_tsv.write("\tPolymorphicAnnotation\tPolyATail\tTSD\tInsertMotif\n")
                    count = 1
                    continue
                row = line.strip().split('\t')
                if row[0] not in ref_contigs:
                    ref_contigs[row[0]] = ref_file.fetch(row[0])
                sequence = row[4]
                read_positions = defaultdict(list)
                #print(row[0:3])
                for read_insert in row[5].split(','):
                    if read_insert == "NA":
                        continue
                    read = read_insert.split(':')[1]
                    orientation = read_insert.split(':')[2]
                    insert_start = int(read_insert.split(':')[3].split('-')[0])
                    insert_end = int(read_insert.split(':')[3].split('-')[1])
                    read_positions[read] = [insert_start, insert_end, orientation]
                q_list[count%threads].put([row, read_positions])
                count += 1
            print("Read in Data")
            thread_list = [None] *threads
            for i in range(threads):
                thread_list[i] = myThread2(i, q_list[i], seen, control_sample, centromere_positions, telomere_positions, bam_read_positions, ref_aligner, t_lock, min_mapq, bam_ref_positions, chrom_list, min_reads, softclips, ref_contigs, out_tsv)
                thread_list[i].start()
            print("Setup Threads")
            for i in range(threads):
                thread_list[i].join()
            print("Done")

               

