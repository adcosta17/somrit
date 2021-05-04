import argparse
import sys

parser = argparse.ArgumentParser(description="somrit - Somatic Retrotransposon Insertion Toolkit")
subparsers = parser.add_subparsers(help='subprograms')
extract_parser = subparsers.add_parser("extract")
realign_parser = subparsers.add_parser("realign")
classify_parser = subparsers.add_parser("classify")
filter_parser = subparsers.add_parser("filter")

extract_parser.add_argument('--bam', type=str, required=True)
extract_parser.add_argument('--output-tsv', type=str, required=True)
extract_parser.add_argument('--min-insertion-length', type=int, default=100)
extract_parser.add_argument('--min-detected-inclusion-length', type=int, default=50)
extract_parser.add_argument('--min-flank-size', required=False, default=100)
extract_parser.add_argument('--min-mapq', required=False, type=int, default=20)
extract_parser.add_argument('--reference-gap-minimum', type=int, default=100)
extract_parser.add_argument('--fastq-file', type=str, required=True)

realign_parser.add_argument('--bam-list', type=str, required=True)
realign_parser.add_argument('--tsv-list', type=str, required=True)
realign_parser.add_argument('--fastq-list', type=str, required=True)
realign_parser.add_argument('--output-tsv', type=str, required=True)
realign_parser.add_argument('--reference-genome', type=str, required=True)
realign_parser.add_argument('--output-bam', type=str, required=False, default="")
realign_parser.add_argument('--cluster-window', type=int, required=False, default=1000)
realign_parser.add_argument('--reference-window', type=int, required=False, default=5000)

args = parser.parse_args()

if len(sys.argv) < 2:
    print("\nPlease select an option from: [extract, realign, classify, filter]\n")
elif sys.argv[1] == "extract":
    print("Extracting Insertions")
    from src.extract import extract_candidate_insertions
    extract_candidate_insertions(args.bam, args.output_tsv, args.min_insertion_length, args.min_detected_inclusion_length, args.min_flank_size, args.min_mapq, args.reference_gap_minimum, args.fastq_file)
elif sys.argv[1] == "realign":
    #print("Realigning Insertions")
    from src.realign import realign_candidate_insertions
    realign_candidate_insertions(args.bam_list, args.tsv_list, args.fastq_list, args.output_tsv, args.output_bam, args.cluster_window, args.reference_genome, args.reference_window)
elif sys.argv[1] == "classify":
    pass
elif sys.argv[1] == "filter":
    pass

