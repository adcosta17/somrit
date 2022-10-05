import argparse
import sys
import pysam
from ctypes import cdll
from ctypes import c_char_p, c_int
import os
import os.path
from collections import defaultdict

htslib = cdll.LoadLibrary('./htslib/libhts.so')
lib = cdll.LoadLibrary('./librealign.so')

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_chr_list():
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"


def merge_samples(args, bams, tsvs):
    samples = {}
    sample_input_bam = {}
    sample_output_bam = {}
    seen = {}
    input_bams = args.bam_list.split(',')
    output_bams = args.output_bam_list.split(',')
    assert len(input_bams) == len(output_bams), "The length of the input bam list and output bam list should match, if the output bam list is provided"
    assert len(input_bams) > 0, "Must have at least one input bam"
    # Setup the naming dictionaries so we can link input and output bams to the same sample
    filedata = {}
    header = pysam.AlignmentFile(input_bams[0]).header
    for i in range(len(input_bams)):
        path = input_bams[i].split('/')
        sample = path[len(path)-1].split(".bam")[0]
        print(sample)
        samples[sample] = defaultdict(int)
        sample_input_bam[sample] = input_bams[i]
        sample_output_bam[sample] = output_bams[i]
        filedata[sample] = pysam.AlignmentFile(output_bams[i], "wb", header=header)
        seen[sample] = {}
    print(1)
    print(bams)
    for bam in bams:
        print(bam)
        reader = pysam.AlignmentFile(bam)
        for record in reader.fetch():
            # Get the read group tag
            tags = record.get_tags()
            for t in tags:
                if t[0] == "RG":
                    sample = t[1]
                    if sample in filedata:
                        filedata[sample].write(record)
                        seen[sample][record.query_name] = 1
    # Write out the rest of the reads for each sample that weren't re-aligned 
    print(2)
    for sample in samples:
        print(sample)
        reader = pysam.AlignmentFile(sample_input_bam[sample])
        for record in reader.fetch():
            if record.query_name not in seen[sample]:
                filedata[sample].write(record)
    # merge and clean up the files
    #for bam in bams:
    #    os.remove(bam)
    for file in filedata.values():
        file.close()
    print("Done Bams")
    #with open(args.output_dir+"/"+args.tsv_prefix+".tsv", "w") as myfile:
    #    for tsv in tsvs:
    #        with open(tsv, 'r') as in_tsv:
    #            count = 0
    #            for line in in_tsv:
    #                if count == 0:
    #                    if len(bams) == 0:
    #                        myfile.write(line.strip()+"\n")
    #                    count += 1
    #                    continue
    #                myfile.write(line)
    #for tsv in tsvs:
    #    os.remove(tsv)


def merge_all(args, bams, tsvs):
    pysam.merge(args.output_dir+"/"+args.bam_prefix+".tmp.bam", *bams)
    pysam.index(args.output_dir+"/"+args.bam_prefix+".tmp.bam")
    for bam in bams:
        os.remove(bam)
    bams = [args.output_dir+"/"+args.bam_prefix+".tmp.bam"]
    # Read in the reads in the tmp bam. Get their names, and then write out the remaining records
    seen = {}
    header = pysam.AlignmentFile(args.bam_list.split(',')[0]).header
    reader = pysam.AlignmentFile(args.output_dir+"/"+args.bam_prefix+".tmp.bam")
    for record in reader.fetch():
        seen[record.query_name] = 1
    with pysam.AlignmentFile(args.output_dir+"/"+args.bam_prefix+".rest.bam", "wb", header=header) as outf:
        for bam in args.bam_list.split(','):
            reader = pysam.AlignmentFile(bam)
            for record in reader.fetch():
                if record.query_name not in seen:
                    outf.write(record)
    bams.append(args.output_dir+"/"+args.bam_prefix+".rest.bam")
    pysam.merge(args.output_dir+"/"+args.bam_prefix+".bam", *bams)
    os.remove(args.output_dir+"/"+args.bam_prefix+".tmp.bam")
    os.remove(args.output_dir+"/"+args.bam_prefix+".tmp.bam.bai")
    os.remove(args.output_dir+"/"+args.bam_prefix+".rest.bam")
    with open(args.output_dir+"/"+args.tsv_prefix+".tsv", "w") as myfile:
        for tsv in tsvs:
            with open(tsv, 'r') as in_tsv:
                count = 0
                for line in in_tsv:
                    if count == 0:
                        if len(bams) == 0:
                            myfile.write(line.strip()+"\n")
                        count += 1
                        continue
                    myfile.write(line)
    for tsv in tsvs:
        os.remove(tsv)

class ReAlign(object):
  
    # constructor
    def __init__(self):
        # attribute
        self.obj = lib.ReAlign_new()
  
    # define method
    def realign_all(self, bam_list, tsv_list, fastq_list, samples_list, output_tsv, output_bam, cluster_window, reference_genome, gap_open, gap_extend, min_mapq, chr_list, pass_only, depth_filter, threads, max_mem, high_mem, include_haps, max_window_size, max_insert_size):
        lib.ReAlign_realign_all(self.obj, bam_list, tsv_list, fastq_list, samples_list, output_tsv, output_bam, cluster_window, reference_genome, gap_open, gap_extend, min_mapq, chr_list, pass_only, depth_filter, threads, max_mem, high_mem, include_haps, max_window_size, max_insert_size)

    def haps_only(self, bam_list, tsv_list, fastq_list, samples_list, output_fa, cluster_window, reference_genome, gap_open, gap_extend, min_mapq, chr_list, pass_only, depth_filter, threads, max_mem, high_mem, include_haps, max_window_size, max_insert_size):
        lib.ReAlign_haps_only(self.obj, bam_list, tsv_list, fastq_list, samples_list, output_fa, cluster_window, reference_genome, gap_open, gap_extend, min_mapq, chr_list, pass_only, depth_filter, threads, max_mem, high_mem, include_haps, max_window_size, max_insert_size)

parser = argparse.ArgumentParser(description="somrit - Somatic Retrotransposon Insertion Toolkit")
subparsers = parser.add_subparsers(help='subprograms')
extract_parser = subparsers.add_parser("extract")
realign_parser = subparsers.add_parser("realign")
classify_parser = subparsers.add_parser("classify")
filter_parser = subparsers.add_parser("filter")
merge_parser = subparsers.add_parser("merge")

extract_parser.add_argument('--bam', type=str, required=True)
extract_parser.add_argument('--output-tsv', type=str, required=True)
extract_parser.add_argument('--output-merged', type=str, required=True)
extract_parser.add_argument('--min-insertion-length', type=int, default=100)
extract_parser.add_argument('--min-detected-inclusion-length', type=int, default=50)
extract_parser.add_argument('--min-flank-size', required=False, default=100)
extract_parser.add_argument('--min-mapq', required=False, type=int, default=20)
extract_parser.add_argument('--reference-gap-minimum', type=int, default=100)
extract_parser.add_argument('--fastq-file', type=str, required=True)
extract_parser.add_argument('--threads', type=int, required=False, default=1)

realign_parser.add_argument('--bam-list', type=str, required=True)
realign_parser.add_argument('--tsv-list', type=str, required=True)
realign_parser.add_argument('--fastq-list', type=str, required=True)
realign_parser.add_argument('--output-dir', type=str, required=True)
realign_parser.add_argument('--tsv-prefix', type=str, required=True)
realign_parser.add_argument('--bam-prefix', type=str, required=True)
realign_parser.add_argument('--reference-genome', type=str, required=True)
realign_parser.add_argument('--cluster-window', type=int, required=False, default=1000)
realign_parser.add_argument('--gap-open', type=int, required=False, default=11)
realign_parser.add_argument('--gap-extend', type=int, required=False, default=1)
realign_parser.add_argument('--threads', type=int, required=False, default=1)
realign_parser.add_argument('--min-mapq', type=int, required=False, default=20)
realign_parser.add_argument('--max-mem', type=int, required=False, default=5000000000)
realign_parser.add_argument('--chromosome-list', type=str, required=False, default="all_main")
realign_parser.add_argument("--filter-pass", type=str2bool, nargs='?',const=True, default=False)
realign_parser.add_argument('--max-depth', type=int, required=False, default=1000)
realign_parser.add_argument("--filter-depth", type=str2bool, nargs='?',const=True, default=False)
realign_parser.add_argument("--include-haps", type=str2bool, nargs='?',const=True, default=False)
realign_parser.add_argument("--high-mem", type=str2bool, nargs='?',const=True, default=False)
realign_parser.add_argument("--only-haps", type=str2bool, nargs='?',const=True, default=False)
realign_parser.add_argument("--only-realign", type=str2bool, nargs='?',const=True, default=False)
realign_parser.add_argument('--max-window-size', type=int, required=False, default=10000)
realign_parser.add_argument('--max-insert-size', type=int, required=False, default=25000)

merge_parser.add_argument('--bam-prefix', type=str, required=True)
merge_parser.add_argument('--output-dir', type=str, required=True)
merge_parser.add_argument('--output-bam-list', type=str, required=False)
merge_parser.add_argument('--tsv-prefix', type=str, required=True)
merge_parser.add_argument('--chromosome-list', type=str, required=False, default="all_main")
merge_parser.add_argument('--bam-list', type=str, required=True)

classify_parser.add_argument('--bam-list', type=str, required=True)
classify_parser.add_argument('--sample-list', type=str, required=True)
classify_parser.add_argument('--realign-tsv', type=str, required=False)
classify_parser.add_argument('--tsv-list', type=str, required=True)
classify_parser.add_argument('--output-tsv', type=str, required=True)
classify_parser.add_argument('--annotation-file', type=str, required=True)
classify_parser.add_argument('--fastq-list', type=str, required=True)

filter_parser.add_argument('--input-tsv', type=str, required=True)
filter_parser.add_argument('--output-tsv', type=str, required=True)
filter_parser.add_argument('--centromeres', type=str, required=True)
filter_parser.add_argument('--telomeres', type=str, required=True)
filter_parser.add_argument('--bam', type=str, required=True)
filter_parser.add_argument('--fastq-list', type=str, required=True)
filter_parser.add_argument('--control-sample', type=str, required=False, default="NA")
filter_parser.add_argument('--min-mapq', type=int, required=False, default=20)
filter_parser.add_argument('--min-reads', type=int, required=False, default=1)
filter_parser.add_argument('--reference-genome', type=str, required=True)
filter_parser.add_argument('--cluster-window', type=int, required=False, default=500)
filter_parser.add_argument('--chromosome-list', type=str, required=False, default="all_main")
filter_parser.add_argument('--threads', type=int, required=False, default=1)

args = parser.parse_args()

if len(sys.argv) < 2:
    print("\nPlease select an option from: [extract, realign, classify, filter]\n")
elif sys.argv[1] == "extract":
    print("Extracting Insertions")
    from extract import extract_candidate_insertions
    extract_candidate_insertions(args.bam, args.output_tsv, args.output_merged, args.min_insertion_length, args.min_detected_inclusion_length, args.min_flank_size, args.min_mapq, args.reference_gap_minimum, args.fastq_file, args.threads)
elif sys.argv[1] == "realign":
    print("Realigning Insertions")
    #from src.realign import realign_candidate_insertions
    #realign_candidate_insertions(args.bam_list, args.tsv_list, args.fastq_list, args.output_tsv, args.output_bam, args.cluster_window, args.reference_genome, args.reference_window, args.threads)
    bams = args.bam_list.split(',')
    tsvs = args.tsv_list.split(',')
    fastqs = args.fastq_list.split(',')
    assert len(bams) == len(tsvs), "The length of the tsv list and bam list should match"
    assert len(bams) == len(fastqs), "The length of the fastq list and bam list should match"
    samples = []
    for i in range(len(bams)):
        path = bams[i].split('/')
        sample = path[len(path)-1].split(".bam")[0]
        samples.append(sample)
    # Get the set of chromososmes to use as csv
    chr_to_use = ""
    if args.chromosome_list == "all_main":
        chr_to_use = get_chr_list()
    elif args.chromosome_list != "all":
        chr_to_use = args.chromosome_list
    pass_only = 0
    if args.filter_pass:
        pass_only = 1
    depth_filter = 0
    if args.filter_depth:
        depth_filter = int(args.max_depth)
    include_haps = 0
    if args.include_haps:
        include_haps = 1
    only_haps = 0
    if args.only_haps:
        only_haps = 1
    high_mem = 0
    if args.high_mem:
        high_mem = 1
    if os.path.isfile(args.output_dir+"/"+args.bam_prefix+".bam"):
        # ReAlignment already done here
        exit(0)
    # Use the C++ library for realign
    if args.only_haps:
        # Get just a fasta of alternative haplotypes based on insertions. Don't care about re-aligning reads to the haplotypes
        for chrom in chr_to_use.split(','):
            r = ReAlign()
            out_fasta = args.output_dir+"/"+args.tsv_prefix+"."+chrom+".fa"
            r.haps_only(c_char_p(args.bam_list.encode('utf-8')), c_char_p(args.tsv_list.encode('utf-8')), c_char_p(args.fastq_list.encode('utf-8')), c_char_p(','.join(samples).encode('utf-8')), c_char_p(out_fasta.encode('utf-8')), c_int(int(args.cluster_window)), c_char_p(args.reference_genome.encode('utf-8')), c_int(int(args.gap_open)), c_int(int(args.gap_extend)), c_int(int(args.min_mapq)), c_char_p(chrom.encode('utf-8')), c_int(pass_only), c_int(depth_filter), c_int(int(args.threads)), c_int(int(args.max_mem)), c_int(high_mem), c_int(include_haps), c_int(args.max_window_size), c_int(args.max_insert_size))
            with open(args.output_dir+"/"+args.tsv_prefix+".fa", "a") as myfile:
                with open(out_fasta, 'r') as in_fa:
                    for line in in_fa:
                        myfile.write(line)
            os.remove(out_fasta)
    elif args.only_realign:
        # only do the re-alignment step, don't merge
        for chrom in chr_to_use.split(','):
            r = ReAlign()
            out_tsv = args.output_dir+"/"+chrom+"."+args.tsv_prefix+".tsv"
            out_sam = args.output_dir+"/"+chrom+"."+args.bam_prefix+".sam"
            if os.path.isfile(args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam"):
                bams.append(args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam")
                tsvs.append(out_tsv)
                print("found "+chrom)
                continue
            r.realign_all(c_char_p(args.bam_list.encode('utf-8')), c_char_p(args.tsv_list.encode('utf-8')), c_char_p(args.fastq_list.encode('utf-8')), c_char_p(','.join(samples).encode('utf-8')), c_char_p(out_tsv.encode('utf-8')), c_char_p(out_sam.encode('utf-8')), c_int(int(args.cluster_window)), c_char_p(args.reference_genome.encode('utf-8')), c_int(int(args.gap_open)), c_int(int(args.gap_extend)), c_int(int(args.min_mapq)), c_char_p(chrom.encode('utf-8')), c_int(pass_only), c_int(depth_filter), c_int(int(args.threads)), c_int(int(args.max_mem)), c_int(high_mem), c_int(include_haps), c_int(args.max_window_size), c_int(args.max_insert_size))
            pysam.sort("-o", args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam", args.output_dir+"/"+chrom+"."+args.bam_prefix+".sam")
            os.remove(out_sam)
    else:
        # Do each chrom one by one
        bams = []
        tsvs = []
        for chrom in chr_to_use.split(','):
            r = ReAlign()
            out_tsv = args.output_dir+"/"+chrom+"."+args.tsv_prefix+".tsv"
            out_sam = args.output_dir+"/"+chrom+"."+args.bam_prefix+".sam"
            if os.path.isfile(args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam"):
                bams.append(args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam")
                tsvs.append(out_tsv)
                print("found "+chrom)
                continue
            r.realign_all(c_char_p(args.bam_list.encode('utf-8')), c_char_p(args.tsv_list.encode('utf-8')), c_char_p(args.fastq_list.encode('utf-8')), c_char_p(','.join(samples).encode('utf-8')), c_char_p(out_tsv.encode('utf-8')), c_char_p(out_sam.encode('utf-8')), c_int(int(args.cluster_window)), c_char_p(args.reference_genome.encode('utf-8')), c_int(int(args.gap_open)), c_int(int(args.gap_extend)), c_int(int(args.min_mapq)), c_char_p(chrom.encode('utf-8')), c_int(pass_only), c_int(depth_filter), c_int(int(args.threads)), c_int(int(args.max_mem)), c_int(high_mem), c_int(include_haps), c_int(args.max_window_size), c_int(args.max_insert_size))
            pysam.sort("-o", args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam", args.output_dir+"/"+chrom+"."+args.bam_prefix+".sam")
            bams.append(args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam")
            os.remove(out_sam)
            tsvs.append(out_tsv)
        merge_all(args, bams, tsvs)
elif sys.argv[1] == "merge":
    print("Merging ReAligned Insertions")
    chr_to_use = get_chr_list()
    if args.chromosome_list == "all_main":
        chr_to_use = get_chr_list()
    elif args.chromosome_list != "all":
        chr_to_use = args.chromosome_list
    print(chr_to_use)
    bams = []
    tsvs = []
    for chrom in chr_to_use.split(','):
        out_tsv = args.output_dir+"/"+chrom+"."+args.tsv_prefix+".tsv"
        if os.path.isfile(args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam"):
            print(args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam")
            bams.append(args.output_dir+"/"+chrom+"."+args.bam_prefix+".bam")
            tsvs.append(out_tsv)
    if args.output_bam_list:
        # get a list of reads from each input bam, then generate an output bam for each sample
        merge_samples(args, bams, tsvs)
    else:
        merge_all(args, bams, tsvs)
elif sys.argv[1] == "classify":
    print("Classifing Insertions")
    from classify import classify_insertions
    rt = None
    if args.realign_tsv:
        rt = args.realign_tsv
    samples = args.sample_list.split(',')
    classify_insertions(args.tsv_list, args.bam_list, args.annotation_file, args.output_tsv, args.fastq_list, rt, samples)
elif sys.argv[1] == "filter":
    print("Filtering Insertions")
    from filter import filter_insertions
    centromeres = None
    if args.centromeres:
        centromeres = args.centromeres
    telomeres = None
    if args.telomeres:
        telomeres = args.telomeres
    chr_to_use = ""
    if args.chromosome_list == "all_main":
        chr_to_use = get_chr_list()
    else:
        chr_to_use = args.chromosome_list
    threads = args.threads - 1
    if threads <= 0:
        threads = 1
    filter_insertions(args.input_tsv, args.output_tsv, args.bam, args.fastq_list, centromeres, telomeres, args.control_sample, args.min_mapq, args.reference_genome, args.cluster_window, chr_to_use, args.min_reads, threads)
    pass

