#include <iostream>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <cmath> 
#include "ReAlign.hpp"
#include "TsvRecord.hpp"
#include <stdexcept>
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "parasail.h"

struct AlignmentResult {
    int score;
    std::string probe_name;
    std::string query;
    std::string ref;
    std::string comp;
    std::string oligo;
    char orientation;
    char* cigar;
    unsigned char* alignment;
    int* endLocation;
    int alignmentLength;
};

std::string dna_reverse_complement(std::string seq) {
    reverse(seq.begin(),seq.end());
    for (std::size_t i = 0; i < seq.length(); ++i){
        switch (seq[i]){
        case 'A':
            seq[i] = 'T';
            break;
        case 'C':
            seq[i] = 'G';
            break;
        case 'G':
            seq[i] = 'C';
            break;
        case 'T':
            seq[i] = 'A';
            break;
        }
    }
    return seq;
}

parasail_matrix_t *matrix = parasail_matrix_create("ACGT", 5, -1);

inline AlignmentResult align_parasail(std::string reads, std::string control, parasail_matrix_t *matrix) {
    AlignmentResult complete_result;
    std::string reverse = dna_reverse_complement(reads);
    
    parasail_result_t* result = parasail_sg_trace_scan_16(reads.c_str(),reads.length(),control.c_str(),control.length(),11,1,matrix);
    parasail_traceback_t* traceback = parasail_result_get_traceback(result,reads.c_str(), reads.length(), control.c_str(), control.length(),matrix,'|','*','*');
    parasail_result_t* result_reverse = parasail_sg_trace_scan_16(reverse.c_str(),reverse.length(), control.c_str(), control.length(),11,1,matrix);
    parasail_traceback_t* traceback_reverse = parasail_result_get_traceback(result_reverse,reverse.c_str(), reverse.length(), control.c_str(), control.length(),matrix,'|','*','*');
    parasail_cigar_t* cigar = result->score > result_reverse->score ? parasail_result_get_cigar(result, reads.c_str(), reads.length(), control.c_str(), control.length(), matrix) : parasail_result_get_cigar(result_reverse, reverse.c_str(), reverse.length(), control.c_str(), control.length(), matrix);
    char* cigar_decoded = parasail_cigar_decode(cigar);
    
    if(result->score > result_reverse->score) {
	
        complete_result.score = result->score;
        complete_result.ref = traceback ->ref;
        complete_result.comp = traceback->comp;
        complete_result.query = traceback->query;
        complete_result.orientation = '+';
        complete_result.cigar = cigar_decoded;
    }
    else {
        complete_result.score = result_reverse->score;
        complete_result.ref = traceback_reverse ->ref;
        complete_result.comp = traceback_reverse->comp;
        complete_result.query = traceback_reverse->query;
        complete_result.orientation = '-';
        complete_result.cigar = cigar_decoded;
    }
    parasail_traceback_free(traceback);
    parasail_result_free(result);
    parasail_traceback_free(traceback_reverse);
    parasail_result_free(result_reverse);

    return complete_result;
}

// Get a map of all positions for reads in bam, for each window get the read IDs and store id in the sample dict
std::unordered_map<std::string, std::map<int, std::vector<std::pair<std::string,std::string>>>> get_reads_in_window(std::unordered_map<std::string, std::string>& samples_to_bam, std::unordered_map<std::string, std::map<int,int>>& positions, int window){
	using namespace std;
	unordered_map<string, map<int, vector<pair<string,string>>>> positions_to_reads;
	for(unordered_map<string, string>::iterator it = samples_to_bam.begin(); it != samples_to_bam.end(); ++it){
		string sample = it.first;
		string bam = it.second;
		// Open up the bam
		auto bam_file = sam_open(bam.c_str(), file_type.c_str());
		auto bam_header = sam_hdr_read(bam_file);
		auto bam_record = bam_init1();
		while (sam_read1(bam_file, bam_header, bam_record) >= 0){
			// Compute the reference positions
			int start = (int) bam_record->core.pos - window;
			string chrom = string(bam_header->target_name[aln->core.tid]);
			auto cigar = bam_get_cigar(bam_record);
			string q_name = string(bam_get_qname(bam_record));
			int end = start;
			for (int k = 0; k < bam_record->core.n_cigar; k++){
				int op = bam_cigar_op(cigar[k]);
				int ol = bam_cigar_oplen(cigar[k]);
				if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF){
					end += ol;
				}
			}
			end += window;
			// Need to now check if there are any entries within positions between start and end
			auto lower = positions[chrom].lower_bound(start);
    		auto upper = positions[chrom].upper_bound(end);
    		if(lower != upper){
    			// have a hit
    			for(auto it2 = lower; it2 != upper; ++it2){
    				int window = it2.first;
    				positions_to_reads[chrom][window].push_back(make_pair(q_name, sample));
    			}
    		}
		}
	}
	return positions_to_reads;
}


// For now void, will return some number of maps later
std::unordered_map<std::string, std::string> get_hap_set(std::string& chrom, int start, int window, std::vector<std::pair<std::string, std::string>>& reads_in_window, std::unordered_map<std::string, std::string>& samples_to_fastq, std::string& ref_sequence, std::unordered_map<std::string, std::vector<TsvRecord>>& inserts_per_sample, int ref_window){
	// For each insert in window, splice insert into ref sequence to create haplotype
	using namespace std;
	unordered_map<string, string> hap_set; 
	for(unordered_map<string, vector<TsvRecord>>::iterator it = inserts_per_sample.begin(); it != inserts_per_sample.end(); ++it){
		string sample = it.first;
		for(size_t i = 0; i < it.second.size(); i++){
			// Compute insertion position on ref seq and splice in insertion
			int ref_start = start - ref_window;
			int insert_start = it.second[i].ref_start - ref_start;
			int insert_end = it.second[i].ref_end - ref_start;
			// Insert the sequence at insert_start
			string hap = ref_sequence.substr(0, insert_start) + it.second[i].seq + ref_sequence.substr(insert_end);
			hap_set[it.second[i].read_id] = hap;
		}
	}
	return hap_set;
}

std::string get_best_hap(std::unordered_map<std::string, std::vector<TsvRecord>>& inserts_per_sample, std::unordered_map<std::string, std::string>& hap_set){
	// For each insert in inserts_per_sample, align it to each hap
	using namespace std;
	unordered_map<string, int> scores;
	for(unordered_map<string, vector<TsvRecord>>::iterator it = inserts_per_sample.begin(); it != inserts_per_sample.end(); ++it){
		string sample = it.first;
		auto read_index = fai_load(samples_to_fastq[sample].c_str());
		for(size_t i = 0; i < it.second.size(); i++){
			// Compute insertion position on ref seq and splice in insertion
			char* tmp_read = fai_fetch(read_index,it.second[i].read_id.c_str());
			string read_sequence = string(tmp_read);
			for(unordered_map<string,string>::iterator it2 = hap_set.begin(); it2 != hap_set.end(); it2++){
				string hap_sequence = it2.second;
				AlignmentResult result = align_parasail(read_sequence, hap_sequence, matrix);
				scores[it2.first] += result.score;
			}
			free(tmp_read);
		}
	}
	// Compute and return the highest scoring hap
	int best_score = 0;
	string best_hap = "";
	for(unordered_map<string, int>::iterator it = scores.begin(); it != scores.end(); it++){
		if(it.second > best_score){
			best_score = it.second;
			best_hap = it.first;
		}
	}
	return best_hap;
}

void ReAlign::realign_all(std::string& bam_list, std::string& tsv_list, std::string& fastq_list, std::string& sample_list, std::string& output_tsv, std::string& output_bam, int window, std::string& ref, int ref_window){
	using namespace std;
	// Read in the lists
	vector<string> samples;
	vector<string> bams;
	vector<string> tsvs;
	vector<string> fastqs;
	stringstream sss(sample_list);
	for (int i; sss >> i;) {
		samples.push_back(i);
		if (sss.peek() == ',')
			sss.ignore();
	}
	unordered_map<string, string> samples_to_bam;
	stringstream ssb(bam_list);
	int count = 0;
	for (int i; ssb >> i;) {
		bams.push_back(i);
		samples_to_bam[samples[count]] = i;
		count++;
		if (ssb.peek() == ',')
			ssb.ignore();
	}
	unordered_map<string, string> samples_to_fastq;
	stringstream ssf(fastq_list);
	int count = 0;
	for (int i; ssf >> i;) {
		fastqs.push_back(i);
		samples_to_fastq[samples[count]] = i;
		count++;
		if (ssf.peek() == ',')
			ssf.ignore();
	}
	unordered_map<string, string> samples_to_tsv;
	stringstream sst(tsv_list);
	int count = 0;
	for (int i; sst >> i;) {
		tsvs.push_back(i);
		samples_to_tsv[samples[count]] = i;
		count++;
		if (sst.peek() == ',')
			sst.ignore();
	}
	// Open each tsv and read in data
	all_tsv_records = unordered_map<string, map<int, vector<TsvRecord>>>;
	positions_per_chrom = unordered_map<string, pair<int,int>>;
	for(unordered_map<string, string>::iterator it = samples_to_tsv.begin(); it != samples_to_tsv.end(); ++it){
		string sample = it.first;
		string tsv_file = it.second;
		// Open and read tsv file
		ifstream infile(it.second);
		string chrom, read, seq, annotation;
		int ref_start, ref_end, read_start, read_end, merged
		// First line is header
		getline(infile, chrom);
		while(infile >> chrom >> ref_start >> ref_end >> read >> read_start >> read_end >> merged >> seq >> annotation){
			// For each line create a TSV record class that has everything set correctly if it is passing
			if(annotation != "PASS"){
				continue;
			}
			// Setup tsv record
			r = TsvRecord(read, read_start, read_end, chrom, ref_start, ref_end, merged, seq, sample);
			all_tsv_records[chrom][ref_start].push_back(r);
			if(positions_per_chrom.count(chrom) == 0){
				positions_per_chrom[chrom] = make_pair(numeric_limits<int>::max(), 0);
			}
			if(ref_start < positions_per_chrom[chrom].first){
				positions_per_chrom[chrom].first = ref_start;
			}
			if(ref_end > positions_per_chrom[chrom].second){
				positions_per_chrom[chrom].second = ref_end;
			}

		}
	}
	// Iterate through each window and get a vector of all tsv_records in that window
	unordered_map<string, map<int,int>> positions;
	for(unordered_map<string, map<int, vector<TsvRecord>>>::iterator it = all_tsv_records.begin(); it != all_tsv_records.end(); ++it){
		string chrom = it.first;
		map<int,int> tmp;
		positions[chrom] = tmp;
		// TODO : Get Reference sequence from the ref
		ref_sequence = "";
		int count = positions_per_chrom[chrom].first;
		while(count < positions_per_chrom[chrom].second){
			auto lower = it.second.lower_bound(count);
    		auto upper = it.second.upper_bound(count+window);
    		if(lower != upper){
    			positions[chrom][count] = 1;
    		}
    		count += window;
    	}
    }

    // Setup faidx for reference
    auto ref_index = fai_load(ref.c_str());

    // Pass these positions to get a list of reads/samples per window
    unordered_map<string, map<int, vector<pair<string,string>>>> reads_in_window = get_reads_in_window(samples_to_bam, positions, window);
    // For each window now that we have all reads, get the haplotypes, then align haplotypes to insert supporting reads, and then 
    for(unordered_map<string, map<int, vector<TsvRecord>>>::iterator it = all_tsv_records.begin(); it != all_tsv_records.end(); ++it){
		string chrom = it.first;
		map<int,int> tmp;
		positions[chrom] = tmp;
		// TODO : Get Reference sequence from the ref
		char* tmp_ref = fai_fetch(ref_index,chrom.c_str());
		string ref_sequence = string(tmp_ref);
		int count = positions_per_chrom[chrom].first;
		while(count < positions_per_chrom[chrom].second){
			auto lower = it.second.lower_bound(count);
    		auto upper = it.second.upper_bound(count+window);
    		if(lower != upper){
    			// Have at least one insert in this range
    			// Setup a list per sample
    			inserts_per_sample = unordered_map<string, vector<TsvRecord>>;
    			for(auto it2 = lower; it2 != upper; ++it2){
    				for(int j = 0; j < it2.second.size(); j++){
    					inserts_per_sample[it2.second[j].sample].push_back(it2.second[j]);
    				}
    			}
    			// Run re-align on the window
    			unordered_map<string, string> hap_set = get_hap_set(chrom, count, window, reads_in_window[chrom][count], samples_to_fastq, ref_sequence, inserts_per_sample, ref_window);
    			// Take hap sets and align it to each read with an insertion, Select best haplotype
    			string best_hap= get_best_haps(reads_in_window[chrom][count], hap_set);
    			// Align reads to best haplotype and ref. If read aligns better to haplotype than ref, update record
    			unorderd_map<string, vector<string>> reads_per_sample;
    			for(size_t i = 0; i < reads_in_window[chrom][count].size(); i++){
    				reads_per_sample[reads_in_window[chrom][count][i].first].push_back(reads_in_window[chrom][count][i].second);
    			}
    			AlignmentResult hap_to_ref = align_parasail(best_hap, ref_sequence, matrix);
    			for(unordered_map<string, vector<string>>::iterator it2 = reads_per_sample.begin(); it2 != reads_per_sample.end(); ++it2){
    				auto read_index = fai_load(samples_to_fastq[it2.first].c_str());
    				for(size_t i = 0; i < it2.second.size(); i++){
    					char* tmp_read = fai_fetch(read_index, it2.second[i].c_str());
						string read_sequence = string(tmp_read);
						AlignmentResult ref_result = align_parasail(read_sequence, ref_sequence, matrix);
						AlignmentResult alt_result = align_parasail(read_sequence, best_hap, matrix);
						free(tmp_read);
						if(alt_result.score > ref_result.score){
							// Need to output a tsv record here
						}
    				}
    			}
    		}
    		count += window;
    	}
    	free(tmp_ref);
    }
}



















