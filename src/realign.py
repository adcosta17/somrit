import pysam
import argparse
import sys
import re
import parasail
from intervaltree import Interval, IntervalTree
from collections import defaultdict

scoring_matrix = parasail.matrix_create("ACGT", 5, -1)


def get_id(header, name):
    i = 0
    while i < len(header['SQ']):
        if header['SQ'][i]['SN'] == name:
            return i
        else:
            i += 1
    return 0


def get_hap_set(chrom, start, window, sample_to_fastq, sample_to_bam, ref_fa, insert_list, ref_window):
    haplotype_sequences = defaultdict(list)
    reads_per_sample = defaultdict(list)
    read_sequences = defaultdict(str)
    best_haplotype_per_read= defaultdict(list)
    insert_reads = {}
    for sample in sample_to_bam:
        in_bam = pysam.AlignmentFile(sample_to_bam[sample])
        sam_reader = in_bam.fetch(region=chrom+":"+str(start)+"-"+str(start+window))
        fq = pysam.FastaFile(filename = sample_to_fastq[sample])
        for record in sam_reader:
            # Grab read sequences from fastq
            read_sequence = fq.fetch(record.query_name)
            read_sequences[record.query_name+"\t"+sample] = read_sequence
    # Evaluate the base haplotype 
    base_ref_seq = ref_fa.fetch(region=chrom+":"+str(start-ref_window)+"-"+str(start+ref_window))
    haplotype_sequences["base"] = [base_ref_seq, 0, ""]
    for insert in insert_list:
        # Compute where the insert sequence should be spliced in
        avg_insert_pos = int((int(insert.data[0][1]) + int(insert.data[0][2]))/2)
        insert_pos = ((avg_insert_pos - start) + ref_window)
        haplotype_sequences[insert.data[0][3]] =  [base_ref_seq[:insert_pos]+insert.data[0][7]+base_ref_seq[insert_pos:], insert_pos - ref_window, insert.data[0][7]]
        insert_reads[insert.data[0][3]] = 1
    # Once we have haplotype set, align each read to each haplotype and mark wich one was best scoring
    haplotype_profiles = {}
    for name in haplotype_sequences:
        haplotype_profiles[name] = parasail.profile_create_16(haplotype_sequences[name][0], scoring_matrix)
    for item in read_sequences:
        read = item.split('\t')[0]
        if read not in insert_reads:
            continue
        sample = item.split('\t')[1]
        #print(read)
        read_sequence = read_sequences[read+"\t"+sample]
        for name in haplotype_sequences:
            result = parasail.sw_striped_profile_16(haplotype_profiles[name], read_sequence, 5, 4)
            best_haplotype_per_read[read].append([result.score, name, read, sample])
            #print([result.score, name, read, sample])
    return (best_haplotype_per_read, haplotype_sequences, haplotype_profiles, read_sequences)

def get_realign_window(haplotype_profiles, read_sequences):
    # Get all reads in region per bam
    best_haplotype_per_read = defaultdict(list)
    for item in read_sequences:
        read = item.split('\t')[0]
        sample = item.split('\t')[1]
        # Grab read sequences from fastq
        #print("realign "+item)
        read_sequence = read_sequences[item]
        # Align each read to best haplotypes
        best_score = 0
        best_hap = ""
        for hap in haplotype_profiles:
            result = parasail.sw_striped_profile_16(haplotype_profiles[hap], read_sequence, 5, 4)
            if result.score > best_score:
                best_score = result.score
                best_hap = hap
        best_haplotype_per_read[item] = [best_hap, result]
    return best_haplotype_per_read


def get_insert_pos(result, insert_ref_start_pos, insert_length):
    insert_ref_end_pos = insert_ref_start_pos + insert_length
    cigar = result.cigar.decode.decode('utf-8')
    # Ref seq is the haplotype, query is the read
    cgs = re.findall('[0-9]*[A-Z=]', cigar)
    ref_hap_count = 0
    read_count = 0
    start_pos = 0
    end_pos = 0
    for cg in cgs:
        if cg.endswith("M"):
            ref_hap_count += int(cg[:cg.find("M")])
            read_count += int(cg[:cg.find("M")])
        elif cg.endswith("I"):
            read_count += int(cg[:cg.find("I")])
        elif cg.endswith("D"):
            ref_hap_count += int(cg[:cg.find("D")])
        elif cg.endswith("N"):
            ref_hap_count += int(cg[:cg.find("N")])
        elif cg.endswith("S"):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith("="):
            ref_hap_count += int(cg[:cg.find("=")])
            read_count += int(cg[:cg.find("=")])
        elif cg.endswith("X"):
            ref_hap_count += int(cg[:cg.find("X")])
            read_count += int(cg[:cg.find("X")])
        if ref_hap_count <= insert_ref_start_pos:
            start_pos = read_count
        if ref_hap_count >= insert_ref_end_pos and end_pos == 0:
            end_pos = read_count
    return(start_pos, end_pos)

def get_start_end_pos(result):
    cigar = result.cigar.decode.decode('utf-8')
    # Ref seq is the haplotype, query is the read
    cgs = re.findall('[0-9]*[A-Z=]', cigar)
    ref_hap_count = 0
    read_count = 0
    start_pos = 0
    end_pos = 0
    for cg in cgs:
        if cg.endswith("M"):
            ref_hap_count += int(cg[:cg.find("M")])
            read_count += int(cg[:cg.find("M")])
        elif cg.endswith("I"):
            read_count += int(cg[:cg.find("I")])
        elif cg.endswith("D"):
            ref_hap_count += int(cg[:cg.find("D")])
        elif cg.endswith("N"):
            ref_hap_count += int(cg[:cg.find("N")])
        elif cg.endswith("S"):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith("H"):
            read_count += int(cg[:cg.find("H")])
        elif cg.endswith("="):
            ref_hap_count += int(cg[:cg.find("=")])
            read_count += int(cg[:cg.find("=")])
        elif cg.endswith("X"):
            ref_hap_count += int(cg[:cg.find("X")])
            read_count += int(cg[:cg.find("X")])
        if read_count == 0:
            start_pos = ref_hap_count
        end_pos = ref_hap_count
    return(start_pos, end_pos)

def get_start_pos(result):
    cigar = result.cigar.decode.decode('utf-8')
    # Ref seq is the haplotype, query is the read
    cgs = re.findall('[0-9]*[A-Z=]', cigar)
    ref_hap_count = 0
    read_count = 0
    start_pos = 0
    count = 0
    for cg in cgs:
        if cg.endswith("M"):
            ref_hap_count += int(cg[:cg.find("M")])
            read_count += int(cg[:cg.find("M")])
            count += 1
        elif cg.endswith("I"):
            read_count += int(cg[:cg.find("I")])
            count += 1
        elif cg.endswith("D"):
            ref_hap_count += int(cg[:cg.find("D")])
            if count == 0:
                # Start on a deletion
                start_pos = ref_hap_count
            count += 1
        elif cg.endswith("N"):
            ref_hap_count += int(cg[:cg.find("N")])
            count += 1
        elif cg.endswith("S"):
            read_count += int(cg[:cg.find("S")])
            count += 1
        elif cg.endswith("H"):
            read_count += int(cg[:cg.find("H")])
            count += 1
        elif cg.endswith("="):
            ref_hap_count += int(cg[:cg.find("=")])
            read_count += int(cg[:cg.find("=")])
            count += 1
        elif cg.endswith("X"):
            ref_hap_count += int(cg[:cg.find("X")])
            read_count += int(cg[:cg.find("X")])
            count += 1
    cleaned_cigar = "".join(cgs[1:len(cgs) - 1])
    return(start_pos, cleaned_cigar)

def get_best_haplotypes(best_haplotype_per_read, haplotype_profiles, haplotype_sequences):
    best_haplotype_profiles = {}
    best_haplotype_profiles["base"] = haplotype_profiles["base"]
    # For each read get set of haps that are within 5% of the best to account for small scoring variations
    top_set = defaultdict(list)
    for read in best_haplotype_per_read:
        best_haplotype_per_read[read] = sorted(best_haplotype_per_read[read], key=lambda x: x[0], reverse=True)
        max_score = -1
        for item in best_haplotype_per_read[read]:
            if item[1] == "base":
                continue
            if max_score == -1:
                max_score = item[0]
                continue
            if float(item[0])/max_score > 0.95:
                top_set[read].append(item[1])
    # Top set should have all the best haplotypes for each read
    reads = list(top_set.keys())
    seen = {}
    groups = defaultdict(list)
    count = 0
    for i in range(len(reads)):
        if i not in seen:
            seen[i] = 1
            groups[i].append(reads[i])
        for j in range(len(reads)):
            if i >= j:
                continue
            if j in seen:
                continue
            if top_set[reads[i]] == top_set[reads[j]]:
                seen[j] = 1
                groups[i].append(reads[j])
    # Have the groups 
    for i in groups:
        best_haps = top_set[groups[i][0]]
        if len(best_haps) == 0:
            # The base was the best hap only, Its already been added to the set
            pass
        else:
            # Select the longest hap from the best haps list
            max_len = 0
            max_hap = ""
            for hap in best_haps:
                if len(haplotype_sequences[hap]) > max_len:
                    max_len = len(haplotype_sequences[hap])
                    max_hap = hap
            best_haplotype_profiles[max_hap] = haplotype_profiles[max_hap]
    return best_haplotype_profiles

# Assumes that the bam_list and tsv_list line up, will throw an error if they don't
def realign_candidate_insertions(bam_list, tsv_list, fastq_list, output_tsv, output_bam, window, ref, ref_window):
    # Parse bam and tsv lists. Check if they are same size
    bams = bam_list.split(',')
    tsvs = tsv_list.split(',')
    fastqs = fastq_list.split(',')
    assert len(bams) == len(tsvs), "The length of the tsv list and bam list should match"
    assert len(bams) == len(fastqs), "The length of the fastq list and bam list should match"
    sample_to_fastq = defaultdict(str)
    sample_to_tsv = defaultdict(str)
    sample_to_bam = defaultdict(str)
    for i in range(len(bams)):
        path = bams[i].split('/')
        sample = path[len(path)-1].split(".bam")[0]
        sample_to_tsv[sample] = tsvs[i]
        sample_to_bam[sample] = bams[i]
        sample_to_fastq[sample] = fastqs[i]
    header = pysam.AlignmentFile(bams[0]).header
    # Read in insertions. Only consider passing inserts
    inserts_per_chrom = defaultdict(IntervalTree)
    min_pos_per_chrom = defaultdict(int)
    max_pos_per_chrom = defaultdict(int)
    for sample in sample_to_tsv:
        tsv = sample_to_tsv[sample]
        with open(tsv, 'r') as in_tsv:
            count = 0
            for line in in_tsv:
                if count == 0:
                    count = 1
                    continue
                row = line.split('\t')
                if row[8] != "PASS":
                    continue
                chrom = row[0]
                inserts_per_chrom[chrom][int(row[1]):int(row[2])] = (row, sample)
                if min_pos_per_chrom[chrom] == 0:
                    min_pos_per_chrom[chrom] = int(row[1])
                if int(row[1]) < min_pos_per_chrom[chrom]:
                    min_pos_per_chrom[chrom] = int(row[1])
                if int(row[1]) > max_pos_per_chrom[chrom]:
                    max_pos_per_chrom[chrom] = int(row[1])
    # Re-align each insertion on its own first. If the insert does not map end to end then it may be off
    # We will align these insertions to other passinng inserts that do align end to end in a later step
    with open(output_tsv, 'w') as out_tsv:
        out_tsv.write("Chromosome\tStart\tEnd\tAltHaplotypeRead\tInsertSequence\tSupportingReads\n")
        with pysam.AlignmentFile(output_bam, "wb", header=header) as outf:
            with pysam.FastaFile(filename=ref) as ref_fa:
                # Go over each window starting at min pos
                for chrom in min_pos_per_chrom:
                    chrom_seq = ref_fa.fetch(chrom)
                    for i in range(min_pos_per_chrom[chrom], max_pos_per_chrom[chrom]+window, window):
                        insert_list = inserts_per_chrom[chrom][i:i+window]
                        insert_reads = {}
                        if len(insert_list) == 0:
                            #print("No inserts at "+ str(i))
                            continue
                        print(i, file=sys.stderr)
                        for item in insert_list:
                            insert_reads[item.data[0][3]+"\t"+item.data[1]] = 1
                        # Generate set of alternative haplotypes by splicing in insert sequence
                        best_haplotype_per_read, haplotype_sequences, haplotype_profiles, read_sequences = get_hap_set(chrom, i, window, sample_to_fastq, sample_to_bam, ref_fa, insert_list, ref_window)
                        # Compute best mapping haplotype set
                        best_haplotype_profiles = get_best_haplotypes(best_haplotype_per_read, haplotype_profiles, haplotype_sequences)
                        # Map all reads to the set of best haplotype sequences
                        reads_in_window_to_haplotypes = get_realign_window(best_haplotype_profiles, read_sequences)
                        reads_per_haplotype = defaultdict(list)
                        # Get list of reads that support each haplotype in reads_in_window_to_halpotypes
                        for item in reads_in_window_to_haplotypes:
                            read = item.split('\t')[0]
                            sample = item.split('\t')[1]
                            hap = reads_in_window_to_haplotypes[item][0]
                            reads_per_haplotype[hap].append([read,sample])
                        # For each hap compute the TSV record. Include insertion sequence and list of supporting reads
                        for hap in reads_per_haplotype:
                            if hap == "base":
                                continue
                            # First get the positions of the main entry 
                            tsv_row = [chrom, str(i+haplotype_sequences[hap][1]), str(i+haplotype_sequences[hap][1]+1), hap+'_insert_hap', haplotype_sequences[hap][2]]
                            # Then get the positions for each supporting read baed on the re-alignment for that read
                            reads = []
                            for item in reads_per_haplotype[hap]:
                                # Get the re-alignmnet record and compute the position of the insertion based on parsing the cigar string
                                #result = reads_in_window_to_haplotypes[item[0]+'\t'+item[1]][1]
                                #read_sequence = read_sequences[item[0]+'\t'+item[1]]
                                #upstream_result = parasail.sw_trace_striped_profile_16(upstream_profile, read_sequence, 5, 4)
                                #downstream_result = parasail.sw_trace_striped_profile_16(downstream_profile, read_sequence, 5, 4)
                                #start_pos, end_pos = get_start_end_pos(result, haplotype_sequences[hap][1] + ref_window, len(haplotype_sequences[hap][2]))
                                reads.append(item[0]+'_'+item[1])
                            tsv_row.append(','.join(reads))
                            out_tsv.write('\t'.join(tsv_row)+'\n')
                            # Can then generate a BAM record for each haplotype sequence mapped to the reference
                            # Need to Re-Map haplotype sequence to the reference to get cigar string
                            result = parasail.sg_trace_striped_16(haplotype_sequences[hap][0], chrom_seq[i-len(haplotype_sequences[hap][0]):i+len(haplotype_sequences[hap][0])], 5, 4, scoring_matrix)
                            start_pos, cleaned_cigar = get_start_pos(result)
                            ###################################################################################
                            # Use Different CIGAR string as parasial may not align insert as single insertion #
                            ###################################################################################
                            cleaned_cigar = str(ref_window)+"M"+str(len(haplotype_sequences[hap][2]))+"I"+str(ref_window)+"M"
                            # Create alignment record for the haplotype based on CIGAR string
                            start_pos = i-len(haplotype_sequences[hap][0]) + start_pos
                            a = pysam.AlignedSegment(header=header)
                            a.query_name = hap+'_insert_hap'
                            a.flag = 0
                            a.reference_id = get_id(header, chrom)
                            a.reference_start = start_pos
                            a.mapping_quality = 60
                            a.cigarstring = cleaned_cigar
                            a.query_sequence = haplotype_sequences[hap][0]
                            outf.write(a)



