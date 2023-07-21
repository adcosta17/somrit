# Somatic Retrotransposon Insertion Toolkit (somrit)

Somrit is a toolkit designed to identify and analyze Somatic Retrotransposon Insertions from Long Reads. It consists of 5 modules designed to be run on their own or together in sequence as part of a larger analysis.

## Dependencies

- pysam
- minimap2 & mappy
- IntervalTree
- abPOA
- htslib

## Setup

somrit can be setup in a few short commands. After recursively cloning the repo, a simple make call is needed to generate the realignment library. 

```sh
# Clone the somrit repo
git clone --recursive https://github.com/adcosta17/somrit.git

# Ensure submodule dependencies are also cloned
git submodule update --init --recursive

# Make minimap2 & htslib libraries, can use system htslib version instead if in path.
make htslib/libhts.so
make minimap2/libminimap2.a

# Make the realignment library
make 
```

## Modules

Somrit consists of 4 modules that can be run on their own given the correct input or as part of a larger analysis. Each module is described in more detail below. Each module has its own CLI and is run using the following format: 
```sh
python somrit.py <module_name>
```
### Extract
The extract module takes in an indexed BAM of reads aligned to a reference genome as input. It identifies insertions in the reads relative to the reference and generates a tsv listing the insertions it has identified. By default it identifies insertions at least 50bp in length, with at least 100bp of sequence on the read on either side and a minimum mapping quality of at least 20. 

Some insertions may be missed by the aligner, and may result in a read having two aligned segments that are adjacent to each other on the reference, but split by an insertion on the read. Somrit extract is able to recover these insertions and outputs a list of reads that have had insertions recovered this way. 

```sh
# somrit extract 
python somrit.py extract --bam <input.bam> --output-tsv <output.tsv> --output-merged <merged.txt> --fastq-file {input.fastq> --threads {num_threads}
```
#### Parameters

**--bam** Input BAM File

**--output-tsv** Output TSV of insertions

**--output-merged** Output list of reads whose alignments were merged to recover missed insertions. 

**--fastq-file** Read fastq file

#### Optional Parameters

**--min-insertion-length** Sets the minimum passing insertion size [100bp]

**--min-detection-inclusion-length** Sets the minimum size where we detect an insertion. Insertions larger than this size but smaller than the --min-insertion-length are included in the output but annotated as being to short. Done to include insertions for later filtering [50bp]

**--min-mapq** Sets the minimum mapping quality a read alignment must have [20]

**--min-read-len** Sets the minimum read length to identify insertions from [1000bp]

**--min-flank-size** Sets the minimum flank required on each side of the insertion [500bp]

**--threads** Number of threads to use [1]

### Realign
The realignment module takes as input the output tsv of somrit extract, the same input bam, as well the read fastq and reference genome the reads were aligned to. It supports multiple bams being realigned together, with a comma separated list of tsvs, bams and fastqs passed as arguments. Note for multiple samples the order of each csv must match, so that the first tsv is generated from the first bam with reads from the first fastq. 

By default Realign outputs a single BAM file containing the realigned alignments along with the original alignment records for reads that were not realigned, along with a tsv of realigned insertions. Realign optionally takes in a list of chromosomes if only certain chromosomes require realignment. Run `python somrit.py realign -h` for a full breakdown of the arguments and options available for realignment.

```sh
# somrit realign 
python somrit.py realign --bam-list <input.bam> --tsv-list <input.tsv> --fastq-list <input.fastq> --output-dir <path/to/output/folder> --tsv-prefix <output_tsv_prefix> --bam-prefix <bam_output_prefix> --reference-genome <ref.fa> -- threads <num_threads> 
```
#### Parameters

**--bam-list** Input list of BAM files. Can be a single BAM or multiple

**--tsv-list** Input list of TSV files from somrit extract. Must be matched to the --bam-list

**--fastq-list** Input list of fastq files. Must be matched to the --bam-list

**--output--dir** Location of written output

**--tsv-prefix** Prefix of output tsv

**--bam-prefix** Prefix of output bam

**--reference-genome** Reference genome used for alignments in BAM file

#### Optional Parameters

**--cluster-window** The insertion clustering window. All inserts within a window up to this size will be realigned together [1000bp]

**--min-mapq** The minimum mapping quality used for insertion identification [20]

**--max-haps** The maximum number of alternative haplotypes to consider when computing a multi-haplotype consensus [3]

**--chromosome-list** If a subset of chromosomes is being used rather than all chromosomes, a csv of chromosome names to use is passed

**--filter-pass** Flag to only realign passing insertion from the extract step [False]

**--max-depth** The maximum number of insertions considered in a cluster window [1000]

**--filter-depth** Flag to skip windows where the maximum depth has been exceeded [False]

**--high-mem** Flag to allow a higher memory approach for consensus sequence generation with abPOA. By default somrit realign breaks up generation of the consensus sequence into 500bp windows. [False]

**--only-haps** Flag to output only consensus haplotypes in the tsv rather than performing realignment of the insertions [False]

**--only-realign** Flag to only perform local realignment but not to replace the old BAM records and TSV records. By default somrit realign will perform realignment and then output an updated BAM and TSV with realigned read alignments and calls. If this flag is passed the output files will only contain the realigned records not any of the non-realigned BAM or TSV records. [False]

**--max-window-size** Maximum size to grow the cluster window to if adjacent windows have insertions [10000bp]

**--max-insert-size** Maximum size of an insertion to be realigned [25000]

**--threads** Number of threads to use [1]

#### Realignment Methodology

Somrit Realign first clusters insertions found in the reads relative to the reference, within 1000 bp on the reference as default, merging adjacent clusters together up to a maximum genomic window size of 10000 bp. It then computes a consensus within the cluster window over all the reads that support an insertion using abPOA, aligning the consensus sequence to the reference with minimap2 to identify the exact insertion position(s) and sequence(s). The consensus insertions are then applied to the reference at the identified genomic positions to generate an alternative haplotype that contains the insertion sequence. Reads that align to the genomic region spanned by the cluster window are then aligned to both the reference and the alternative haplotype. If a read aligns to the alternative haplotype with a higher alignment score than the reference, Somrit Realign uses the alignment of the read to the alternative haplotype and the alignment of the alternative haplotype to the  reference to guide and project the alignment of the read to the reference, correctly placing the insertion in the read. This process allows for a reduction in the ambiguity of the insertion position between reads that support the same insertion event, and allows for the insertions to be detected in reads that support the insertion event but for which an insertion has not been called by the aligner.  

#### Merge

Somrit merge is a helper module designed to be used only in the event that the realignment module is called on multiple chromosomes individually. If the realign module is run in its default settings then the merge module is implicity called from realign to merge the individual chromosome BAM together into a single output BAM. However if individual calls to realign are made for each chromosome, then the somrit merge module combines the realigned calls together into a single BAM, removing the old alignments for reads that were realigned. If multiple bams are passed into the --output-bam-list is must match the --bam-list input in length. As somrit realign can realign reads across multiple samples the output BAM and TSV files for each chromosome run contain reads and calls from multiple samples. The merge helper module will, if multiple bams are passed in the --output-bam-list argument generate an updated BAM and TSV per sample. 

```sh
# somrit merge 
# The --bam-list, --output-dir, --tsv-prefix and --bam-prefix are input variables and should be consistent between the realign runs being merged
python somrit.py merge --bam-list <input.bam> --output-dir <path/to/output/folder> --tsv-prefix <realign_tsv_prefix> --bam-prefix <realign_bam_prefix> --output-bam-list <merged_bam.bam> 
```

### Classify

Somrit classify takes in the tsv output of somrit extract, and optionally the tsv output of somrit realign if realign has been run. Additionally it takes in the fastq and bam for the sample and an annotation fasta file, listing the retrotransposon consensus sequences it aims to classify insertions as. Somrit Classify aligns insertions sequences to the list of RT consensus sequences using minimap2 and outputs a single tsv of insertions, annotated with the retrotransposon repeat family and sub family they best align to if possible. Note that classify can combine multiple samples together like realign, taking in an paired order sensitive list of bams, fastqs and tsvs. 

```sh
# somrit classify 
python somrit.py classify --bam-list <input.bam> --tsv-list <original.tsv> --realign-tsv <realigned.tsv> --annotation-file <data/repeats.fa> --fastq-list <input.fastq> --output-tsv <output.tsv>
```
#### Parameters

**--bam-list** Input BAM(s) used for realignment

**--realign-tsv** The tsv output from somrit realign 

**--tsv-list** The tsv(s) output from somrit extract. These are the original tsvs for each sample  

**--annotation-file** The annotation file of retrotransposon consensus sequences used to classify insertions to repeat families. An example file is provided in the utils folder.

**-fastq-list** The sample fastq(s)

**--output-tsv** The output-tsv


### Filter

Somrit filter takes in the tsv output of somrit classify along with the BAMs and Fastqs of the samples used to generate the tsv and applies a series of filtering steps to identify and remove mapping artifacts and polymorphic insertions. 

```sh
# somrit filter 
python somrit.py filter --input-tsv <input.tsv> --bam <input.bam> --reference-genome <ref.fa> --centromeres <data/centromeres.bed> --telomeres <data/telomeres.bed> --output-tsv <output.tsv>
```
#### Parameters

**--input-tsv** The input tsv, this tsv file is generated somrit classify 

**--output-tsv** The filtered output tsv file

**--bam** The Realigned BAM for this sample 

**--reference-genome** The Reference genome used for realignment 

#### Optional Parameters

**--centromeres** A BED file of centromere positions [None]

**--telomeres** A BED file of telomere positions [None]

**--control-sample** The sample name of a control sample. If multiple samples are realigned together and one represents a control [None]

**--min-mapq** Minimum Mapping quality of read alignments [20]

**--min-reads** Minimum number of reads supporting an insertion event [1]

**--cluster-window** Window to group insertions for polymorphic classification [500]

**--chromosome-list** List of chromosomes if only some are being considered. Comma separated

**--threads** Number of threads [1]




