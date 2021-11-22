#include <iostream>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <tuple> 
#include <algorithm>
#include <limits>
#include <cmath>
#include <thread>
#include <chrono>
#include <queue>
#include <functional>
#include <mutex>
#include "ReAlign.hpp"
#include <stdexcept>
#include "htslib/htslib/vcf.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include "htslib/htslib/hts.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include "spoa/spoa.hpp"

extern "C" {

#include "minimap.h"
#include "include/abpoa.h"

std::mutex bam_lock;
std::mutex tsv_lock;
std::mutex out_lock;
std::mutex queue_lock;
std::queue<std::tuple<std::string,int,int>> positions_for_threads;
std::mutex seen_lock;
std::unordered_set<std::string> seen_full_reads;
std::mutex hap_lock;
std::mutex faidx_lock;
std::queue<std::pair<std::string, std::string>> alt_hap_queue;
std::unordered_map<std::string, std::unordered_map<std::string, std::string>> alt_hap_per_read_per_sample;
std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<int, int, std::string, int, int, int, std::string>>> consensus_tuples_per_pos;
std::mutex count_lock;
long align_time_count = 0;
long map_time_count = 0;
long project_time_count = 0;
long hap_time_count = 0;

// AaCcGgTtNn ... ==> 0,1,2,3,4 ...
// BbDdEeFf   ... ==> 5,6,7,8 ...
unsigned char _char26_table[256] = {
	 0,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15, 
	16, 17, 18, 19,  20, 21, 22, 23,  24, 25, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26,  0,  5,  1,   6,  7,  8,  2,   9, 10, 11, 12,  13, 14,  4, 15, 
	16, 17, 18, 19,   3, 20, 21, 22,  23, 24, 25, 26,  26, 26, 26, 26, 
	26,  0,  5,  1,   6,  7,  8,  2,   9, 10, 11, 12,  13, 14,  4, 15, 
	16, 17, 18, 19,   3, 20, 21, 22,  23, 24, 25, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26
};

// 0/1/2/3/4=>ACGTN
// 5/6/7/8=>BDEF ...
const char _char256_table[256] = {
	'A', 'C', 'G', 'T',  'N', 'B', 'D', 'E',  'F', 'H', 'I',  'J', 'K', 'L', 'M', 'O',
	'P', 'Q', 'R', 'S',  'U', 'V', 'W', 'X',  'Y', 'Z', '*', '-',  '*', '*', '*', '*',
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', 'A', 'B', 'C',  'D', 'E', 'F', 'G',  'H', 'I', 'J', 'K',  'L', 'M', 'N', 'O', 
	'P', 'Q', 'R', 'S',  'T', 'U', 'V', 'W',  'X', 'Y', 'Z', '*',  '*', '*', '*', '*', 
	'*', 'A', 'B', 'C',  'D', 'E', 'F', 'G',  'H', 'I', 'J', 'K',  'L', 'M', 'N', 'O', 
	'P', 'Q', 'R', 'S',  'T', 'U', 'V', 'W',  'X', 'Y', 'Z', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*'
};


unsigned char _nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


// struct for storing the result of Parasail Alignment
struct AlignmentResult {
	int score;
	std::string probe_name;
	std::string query;
	std::string ref;
	std::string comp;
	std::string oligo;
	char orientation;
	std::string cigar;
	unsigned char* alignment;
	int endLocation;
	int startLocation;
	int alignmentLength;
	int read_start;
	int read_end;
};

// struct to store TSV records read in from file
struct TsvRecord {
	std::string read_id;
	int read_start;
	int read_end;
	std::string chrom;
	int ref_start;
	int ref_end;
	int mapped;
	std::string seq;
	std::string sample;
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

void print_debug(std::string string_to_print){
	out_lock.lock();
	std::cout << string_to_print << std::endl;
	out_lock.unlock();
}

AlignmentResult align_minimap(const std::string& read, mm_idx_t *mi, mm_mapopt_t* mopt, const std::string& output_file, const std::string& read_name, const std::string& ref_name, const bool parse_cigar, mm_tbuf_t *tbuf) {
	AlignmentResult complete_result; 
	mm_reg1_t *reg;
	int j, i, n_reg;
	int best_score = 0;
	char orientation = '+';
	std::string cigar = "";
	int startLocation = 0;
	int endLocation = 0;
	int read_start = 0;
	int read_end = 0;
	if(parse_cigar){
		mopt->flag |= MM_F_CIGAR;
	}
	reg = mm_map(mi, read.length()+1, read.c_str(), &n_reg, tbuf, mopt, 0); // get all hits for the query
	for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
		mm_reg1_t *r = &reg[j];
		if(output_file != ""){
			assert(r->p); // with MM_F_CIGAR, this should not be NULL
			std::ofstream output;
			output.open(output_file, std::ios_base::app);
			output << read_name << "\t" << read.length() << "\t" << r->qs << "\t" << r->qe << "\t" << "+-"[r->rev] << "\t" << ref_name << "\t" << mi->seq[r->rid].len << "\t" << r->rs << "\t" << r->re << "\t" << r->mlen << "\t" << r->blen << "\t" <<r->mapq << "\t" << "cg:Z:";
			for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
				output << (r->p->cigar[i]>>4) << MM_CIGAR_STR[r->p->cigar[i]&0xf];
			output << "\n";
			output.close();
			best_score = 1;
			free(r->p);
		}
		else if(r->score > best_score){
			// have a better scoring hit, use it
			best_score = r->score;
			orientation = "+-"[r->rev];
			startLocation = r->rs;
			endLocation = r->re;
			read_start = r->qs;
			read_end = r->qe;
			cigar = "";
			if(parse_cigar){
				assert(r->p); // with MM_F_CIGAR, this should not be NULL
				for (i = 0; i < r->p->n_cigar; ++i){ // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
					cigar += std::to_string((int)r->p->cigar[i]>>4);
					cigar.push_back(MM_CIGAR_STR[r->p->cigar[i]&0xf]);
				}
				free(r->p);
			}

		}
	}
	free(reg);
	complete_result.score = best_score;
	complete_result.orientation = orientation;
	complete_result.cigar = cigar;
	complete_result.startLocation = startLocation;
	complete_result.endLocation = endLocation;
	complete_result.read_start = read_start;
	complete_result.read_end = read_end;
	return complete_result;
}


AlignmentResult align_minimap_single(const std::string& read, const std::string& ref, const bool parse_cigar, const std::string& output_file, const std::string& read_name, const std::string& ref_name, mm_tbuf_t *tbuf){ 
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mm_set_opt(0, &iopt, &mopt);
	if(output_file != "" || parse_cigar){
		mopt.flag |= MM_F_CIGAR; // perform alignment
	}
	const char* seq = ref.c_str();
	mm_idx_t *mi = mm_idx_str(10,15,0,16,1, &seq, NULL);
	mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
	mopt.cap_kalloc = 100000000;
	mm_set_opt("ava_ont", &iopt, &mopt);
	AlignmentResult complete_result = align_minimap(read, mi, &mopt, output_file, read_name, ref_name, parse_cigar, tbuf);
	mm_idx_destroy(mi);
	return complete_result;
}

std::string align_minimap_multiple(const std::string& read, std::vector<mm_idx_t*>& mi_vec, mm_mapopt_t* mopt, mm_tbuf_t *tbuf, const bool parse_cigar, const bool print){
	if(parse_cigar){
		mopt->flag |= MM_F_CIGAR; // perform alignment
	}
	mm_reg1_t *reg;
	int j, n_reg;
	int best_score = 0;
	std::string best_hap = "";
	for(size_t i = 0; i < mi_vec.size(); i++){
		mm_idx_t* mi = mi_vec[i];
		reg = mm_map(mi, read.length()+1, read.c_str(), &n_reg, tbuf, mopt, 0); // get all hits for the query
		for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
			mm_reg1_t *r = &reg[j];
			if(print){
				printf("%s\t%d\t%d\t%d\t%c\t", "read", read.size(), r->qs, r->qe, "+-"[r->rev]);
				printf("%s\t%d\t%d\t%d\t%d\t%d\t%d", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
				putchar('\n');
			}
			if(r->score > best_score){
				// have a better scoring hit, use it
				best_score = r->score;
				best_hap = std::string(mi->seq[r->rid].name);
				if(parse_cigar){
					assert(r->p); // with MM_F_CIGAR, this should not be NULL
					free(r->p);
				}
			}
		}
		free(reg);
	}
	return best_hap;
}

// Get a map of all positions for reads in bam, for each window get the read IDs and store id in the sample dict 
std::vector<std::pair<bam1_t*,std::string>> get_reads_in_window(const std::unordered_map<std::string, std::string>& samples_to_bam, const int start, const int window, const std::string& chrom){
	using namespace std;
	vector<pair<bam1_t*,string>> positions_to_reads;
	for(auto& it: samples_to_bam){
		string sample = it.first;
		string bam = it.second;
		// Open up the bam
		string file_type = "r";
		if(bam.substr(bam.find_last_of(".") + 1) == "bam") {
				file_type = "rb";
		}

		auto bam_file = sam_open(bam.c_str(), file_type.c_str());
		auto bam_header = sam_hdr_read(bam_file);
		auto idx = sam_index_load(bam_file, bam.c_str());
		if(idx==NULL){
			out_lock.lock();
			cerr << "Error opening BAM index" << endl;
			out_lock.unlock();
			continue;
		}
		string region = chrom+":"+to_string(start)+"-"+to_string(start+window);
		auto iter = sam_itr_querys(idx, bam_header, region.c_str());
		if(iter==NULL){
			out_lock.lock();
			cerr << "Error creating BAM iterator" << endl;
			out_lock.unlock();
			continue;
		}
		auto bam_record = bam_init1();
		while (sam_itr_next(bam_file, iter, bam_record) >= 0) {
			bam1_t* tmp = bam_init1();
			bam1_t* ret = bam_copy1(tmp, bam_record);
			if(ret != tmp){
				out_lock.lock();
				cerr << "Error copying bam record" << endl;
				out_lock.unlock();
				bam_destroy1(tmp);
				continue;
			}
			positions_to_reads.push_back(make_pair(tmp, sample));
		}
		bam_destroy1(bam_record);
		hts_idx_destroy(idx);
		sam_itr_destroy(iter);
		sam_hdr_destroy(bam_header);
		sam_close(bam_file);
	}
	return positions_to_reads;
}


// Gets a set of all possible haplotype HIns within the window
std::unordered_map<std::string, std::tuple<int, int, std::string, int, int, int, std::string>> get_hap_set(const int start, const std::string& ref_sequence, const std::unordered_map<std::string, std::vector<TsvRecord>>& inserts_per_sample, const int ref_window, const int window, int t){
	// For each insert in window, splice insert into ref sequence to create haplotype
	using namespace std;
	unordered_map<string, tuple<int, int, string, int, int, int, string>> hap_set; 
	for(auto it = inserts_per_sample.begin(); it != inserts_per_sample.end(); ++it){
		string sample = it->first;
		for(size_t i = 0; i < it->second.size(); i++){
			// Compute insertion position on ref seq and splice in insertion
			int ref_start = start;
			if(ref_start < 0 ){
				ref_start = 0;
			}
			int window_to_use = ref_window;
			if(ref_start + ref_window > (int) ref_sequence.length()){
				window_to_use = ref_sequence.length() - ref_start - 1;
			} 
			int insert_start = it->second[i].ref_start - ref_start;
			// Insert the sequence at insert_start
			string hap = ref_sequence.substr(ref_start, insert_start) + it->second[i].seq + ref_sequence.substr(ref_start+insert_start, window_to_use);
			hap_set[it->second[i].read_id] = make_tuple(insert_start, it->second[i].seq.length(), it->second[i].seq, it->second[i].ref_start, it->second[i].ref_end,window,it->second[i].read_id);
		}
	}
	return hap_set;
}

// Gets the best haplotype based on the alignment of the reads to the haplotypes 
std::string get_best_hap(std::unordered_map<std::string, std::vector<TsvRecord>>& inserts_per_sample, std::unordered_map<std::string, std::tuple<int, int, std::string, int, int, int, std::string>>& hap_set, std::unordered_map<std::string, std::unordered_map<std::string,std::string>>& read_sequences_per_sample, int gap_open, int gap_extend, int t, const int max_mem, bool edlib, mm_tbuf_t *tbuf){
	// For each insert in inserts_per_sample, align it to each hap
	using namespace std;
	unordered_map<string, int> scores;
	for(auto it2 = hap_set.begin(); it2 != hap_set.end(); it2++){
		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		mm_set_opt(0, &iopt, &mopt);
		const char* seq = get<2>(it2->second).c_str();
		mm_idx_t *mi = mm_idx_str(10,15,0,16,1, &seq, NULL);
		mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
		mopt.cap_kalloc = 100000000;
		mm_set_opt("ava_ont", &iopt, &mopt);
		for(auto it = inserts_per_sample.begin(); it != inserts_per_sample.end(); ++it){
			string sample = it->first;
			for(size_t i = 0; i < it->second.size(); i++){
				string read_seq = read_sequences_per_sample.at(sample).at(it->second[i].read_id);
				AlignmentResult result = align_minimap(read_seq, mi, &mopt, "","","",false, tbuf);
				scores[it2->first] += result.score;
			}
		}
		mm_idx_destroy(mi);
	}
	// Compute and return the highest scoring hap
	int best_score = std::numeric_limits<int>::min();
	if(edlib){
		best_score = std::numeric_limits<int>::max();
	}
	string best_hap = "";
	for(auto it = scores.begin(); it != scores.end(); it++){
		if(edlib){
			if(it->second < best_score){
				best_score = it->second;
				best_hap = it->first;
			}
		} else {
			if(it->second > best_score){
				best_score = it->second;
				best_hap = it->first;
			}
		}
	}
	return best_hap;
}

int parse_cigar_string(const std::string& cigar){
	using namespace std;
	int hap_count = 0;
	int read_count = 0;
	istringstream parser(cigar);
	char op;
	uint32_t ol;
	while(parser >> ol >> op) {
		if (op == 'M' || op == '=' || op == 'X'){
			//Parse both
			hap_count += ol;
			read_count += ol;
		} else if(op == 'I'){
			read_count += ol;
		} else if(op == 'D'){
			hap_count += ol;
		}
	}
	return read_count;
}

std::unordered_map<std::string, std::string> populate_dict(const std::string& item_list, std::vector<std::string>& samples){
	std::unordered_map<std::string, std::string> samples_to_item;
	std::stringstream ssb(item_list);
	int count = 0;
	while( ssb.good() ){
		std::string substr;
		std::getline( ssb, substr, ',' );
		samples_to_item[samples[count]] = substr;
		count++;
	}
	return samples_to_item;
}

bool validate_sequence(std::string& seq){
	using namespace std;
	if(seq.find('b') != std::string::npos || seq.find('d') != std::string::npos || seq.find('e') != std::string::npos || seq.find('f') != std::string::npos || seq.find('h') != std::string::npos || seq.find('i') != std::string::npos || seq.find('j') != std::string::npos || seq.find('k') != std::string::npos || seq.find('l') != std::string::npos || seq.find('m') != std::string::npos || seq.find('z') != std::string::npos || seq.find('o') != std::string::npos || seq.find('p') != std::string::npos || seq.find('q') != std::string::npos || seq.find('r') != std::string::npos || seq.find('s') != std::string::npos || seq.find('u') != std::string::npos || seq.find('v') != std::string::npos || seq.find('w') != std::string::npos || seq.find('x') != std::string::npos || seq.find('y') != std::string::npos){
		return false;
	}
	return true;
}

void get_all_tsv_positions(const std::unordered_map<std::string, std::string>& samples_to_tsv, const std::unordered_map<std::string, std::string>& samples_to_bam, bam_hdr_t*& header, std::unordered_map<std::string, std::map<int, std::vector<TsvRecord>>>& all_tsv_records, std::unordered_map<std::string, std::pair<int,int>>& start_end_positions_per_chrom, htsFile* outfile, bool filter_pass){
	using namespace std;
	int seen_samples = 0;
	for(auto it = samples_to_tsv.begin(); it != samples_to_tsv.end(); ++it){
		string sample = it->first;
		if(seen_samples == 0){
			htsFile *infile = hts_open(samples_to_bam.at(sample).c_str(),"rb");
			header = sam_hdr_read(infile);
			int ret_val = sam_hdr_write(outfile, header);
			if ( ret_val < 0 ) {
				out_lock.lock();
				cerr << "ReAlign: error copying header from input bamfile to new bamfile" << endl;
				out_lock.unlock();
			}
			hts_close(infile);
		}
		seen_samples += 1;
		string tsv_file = it->second;
		// Open and read tsv file
		ifstream infile(it->second);
		string chrom, read, seq, annotation, line;
		int ref_start, ref_end, read_start, read_end;
		float merged;
		int count = 0;
		while(getline(infile, line)){
			if(count == 0){
				count = 1;
				continue;
			}
			// For each line create a TSV record class that has everything set correctly if it is passing
			istringstream iss(line);
			iss >> chrom >> ref_start >> ref_end >> read >> read_start >> read_end >> merged >> seq >> annotation;
			if(filter_pass){
				// check that the annotation contains PASS, if not skip this record
				if(annotation != "PASS"){
					continue;
				}
			}
			// Setup tsv record
			TsvRecord r;
			r.read_id = read;
			r.read_start = read_start;
			r.read_end = read_end;
			r.chrom = chrom;
			r.ref_start = ref_start;
			r.ref_end = ref_end;
			r.mapped = merged;
			if(validate_sequence(seq)){
				r.seq = seq;
			} else {
				r.seq = "N";
			}
			r.sample = sample;
			all_tsv_records[chrom][ref_start].push_back(r);
			if(start_end_positions_per_chrom.count(chrom) == 0){
				start_end_positions_per_chrom[chrom] = make_pair(numeric_limits<int>::max(), 0);
			}
			if(ref_start < start_end_positions_per_chrom[chrom].first){
				start_end_positions_per_chrom[chrom].first = ref_start;
			}
			if(ref_end > start_end_positions_per_chrom[chrom].second){
				start_end_positions_per_chrom[chrom].second = ref_end;
			}
		}
	}
}

std::unordered_map<std::string, std::set<std::pair<int,int>>> get_positions(const std::unordered_map<std::string, std::map<int, std::vector<TsvRecord>>>& all_tsv_records, const std::unordered_map<std::string, std::pair<int,int>>& start_end_positions_per_chrom, const int window){
	using namespace std;
	unordered_map<string, set<pair<int,int>>> positions; // -- Use a type def to a 
	for(auto it = all_tsv_records.begin(); it != all_tsv_records.end(); ++it){
		string chrom = it->first;
		int count = start_end_positions_per_chrom.at(chrom).first;
		int prev = 0;
		int prev_window = 0;
		while(count < start_end_positions_per_chrom.at(chrom).second){
			auto lower = it->second.lower_bound(count);
			auto upper = it->second.upper_bound(count+window);
			if(lower != upper){
				// Get the actual start position of the insert in the window
				int pos = lower->first;
				if(prev == 0){
					prev = pos-window;
					prev_window = 2*window;
				} else {
					//Just saw an insertion in the adjacent window, add it to this one
					prev_window += window;
				}
				count = pos+window;
			} else {
				if(prev != 0){
					prev_window += window;
					positions[chrom].insert(make_pair(prev, prev_window));
					count += window;
				}
				count += window;
				prev = 0;
				prev_window = 0;
			}
		}
		if(prev != 0){
			positions[chrom].insert(make_pair(prev, prev_window));
		}
	}
	return positions;
}

void get_windows(int pos, const std::string chrom, const int window, const std::unordered_map<std::string, std::map<int, std::vector<TsvRecord>>>& all_tsv_records, std::unordered_map<std::string, std::vector<TsvRecord>>& inserts_per_sample, std::unordered_map<std::string, std::vector<bam1_t*>>& alignments_per_read, std::unordered_map<std::string, std::unordered_map<std::string,std::string>>& read_sequences_per_sample, size_t& max_read_len, const std::unordered_map<std::string, faidx_t*>& samples_to_faidx, const std::vector<std::pair<bam1_t*,std::string>>& reads_in_window){
	// Setup a list per sample
	using namespace std;
	auto lower = all_tsv_records.at(chrom).lower_bound(pos);
	auto upper = all_tsv_records.at(chrom).upper_bound(pos+window);
	for(auto it = lower; it != upper; ++it){
		for(size_t j = 0; j < it->second.size(); j++){
			inserts_per_sample[it->second[j].sample].push_back(it->second[j]);
		}
	}
	unordered_map<string, vector<string>> reads_per_sample;
	for(size_t i = 0; i < reads_in_window.size(); i++){
		reads_per_sample[reads_in_window[i].second].push_back(bam_get_qname(reads_in_window[i].first));
		alignments_per_read[bam_get_qname(reads_in_window[i].first)].push_back(reads_in_window[i].first);
	}
	// Compute the maximum read length within the window
	for(auto it2 = reads_per_sample.begin(); it2 != reads_per_sample.end(); ++it2){
		auto read_index = samples_to_faidx.at(it2->first);
		for(size_t i = 0; i < it2->second.size(); i++){
			int size = 0;
			faidx_lock.lock();
			char* tmp_read = fai_fetch(read_index, it2->second[i].c_str(), &size);
			faidx_lock.unlock();
			string read_sequence = string(tmp_read);
			read_sequences_per_sample[it2->first][it2->second[i]] =  read_sequence;
			free(tmp_read);
			if(read_sequence.length() > max_read_len){
				max_read_len = read_sequence.length();
			}
		}
	}
}

int write_bam_record(const std::string& cigar, const std::string& name, const int start_pos, const int orientation, const std::string& chrom, const int op_number, const std::string& seq, htsFile *outfile, bam_hdr_t *header, const int t){
	using namespace std;
	uint32_t *a_cigar = (uint32_t*) malloc(cigar.size()+1);
	size_t a_mem = 0;
	sam_parse_cigar(cigar.c_str(), NULL, &a_cigar, &a_mem);
	bam1_t *new_hap_record = bam_init1();
	int tid = sam_hdr_name2tid(header, chrom.c_str());
	int ret = bam_set1(new_hap_record, name.length(), name.c_str(), orientation, tid, start_pos, 60, op_number, a_cigar, tid, start_pos, 0, seq.length(), seq.c_str(), NULL, 0);
	if(ret < 0){
		out_lock.lock();
		cerr << t << " Unable to generate bam for "<< name << endl;
		out_lock.unlock();
	}
	bam_lock.lock();
	ret = sam_write1(outfile, header, new_hap_record);
	bam_lock.unlock();
	if(ret == -1){
		out_lock.lock();
		cerr << t << " Write to bam failed for " << name << endl;
		out_lock.unlock();
	}
	bam_destroy1(new_hap_record);
	return ret;
}

void get_cigar_length(const std::string& cigar, int insertion_size){
	using namespace std;
	istringstream parser(cigar);
	int hap_count = 0;
	int read_count = 0;
	char op;
	uint32_t ol;
	bool at_insert = false;
	int before_insert = 0;
	int after_insert = 0;
	while(parser >> ol >> op) {
		if (op == 'M' || op == '=' || op == 'X'){
			//Parse both
			hap_count += ol;
			read_count += ol;
			if(at_insert){
				after_insert += ol;
			} else{
				before_insert += ol;
			}
		} else if(op == 'I' || op == 'S'){
			read_count += ol;
			if(ol == insertion_size && op == 'I'){
				at_insert = true;
				continue;
			}
			if(at_insert){
				after_insert += ol;
			} else{
				before_insert += ol;
			}
		} else if(op == 'D'){
			hap_count += ol;
		}
		else {
		}
	}
}

std::pair<int, char> get_projected_op(uint32_t ol, char read_to_hap_op, char hap_to_ref_op){
	using namespace std;
	char op_to_project = ' ';
	if(hap_to_ref_op == 'M'){
		if(read_to_hap_op == 'M'){
			op_to_project = 'M';
		}
		else if(read_to_hap_op == 'I'){
			op_to_project = 'I';
		}
		else if(read_to_hap_op == 'D'){
			op_to_project = 'D';
		}
	} else { // Insert
		if(read_to_hap_op == 'M'){
			op_to_project = 'I';
		}
		else if(read_to_hap_op == 'I'){
			op_to_project = 'I';
		}
		else if(read_to_hap_op == 'D'){
			// Don't do anything here
		}
	}
	pair<int, int> item = make_pair(0, 'M');
	if(op_to_project != ' '){
		item = make_pair(ol, op_to_project);
	}
	return item;
}

std::pair<std::string, int> get_final_cig(std::vector<std::pair<int, char>>& final_cig){
	using namespace std;
	if(final_cig.size() == 0){
		return make_pair("", 0);
	}
	int op_count = 0;
	string final_cig_str = "";
	int prev_count = final_cig[0].first;
	char prev_op = final_cig[0].second;
	for(size_t i = 1; i < final_cig.size(); i++){
		if(final_cig[i].second == prev_op){
			// append to the op
			prev_count += final_cig[i].first;
		} else {
			// Different op, don't need to merge
			final_cig_str += to_string(prev_count)+prev_op;
			op_count += 1;
			prev_count = final_cig[i].first;
			prev_op = final_cig[i].second;
		}
	}
	final_cig_str += to_string(prev_count)+prev_op;
	op_count += 1;
	return make_pair(final_cig_str, op_count);
}

std::vector<std::pair<int,int>> project_alignment(int hap_start, std::vector<std::pair<int,int>> hap_insert_positions, const std::string& alt_cigar, const std::string& hap_name, const std::string read_name, std::string& read_sequence, const char read_orientation, const std::string chrom, htsFile *outfile, bam_hdr_t *header, const int startLocation, const int endLocation, const int readStart, const int readEnd, const int t){
	using namespace std;
	vector<pair<int,int>> ret;

	// -- ReWork to do two cigar string projection with ref, hap and read counts on base level. At each base check if Match/Insert on ref to hap vs M/I/D on hap to read CIGAR.
	// Convert hap_insert_positions to a CIGAR String
	int hap_to_ref_count = startLocation;
	hap_to_ref_count += 1;
	string hap_to_ref_cigar = "";
	int count_cigar_ops = 0;
	for(size_t i = 0; i < hap_insert_positions.size(); i++){
		// First get the match and then the insert
		int insert_start = hap_insert_positions[i].first;
		int insert_end = hap_insert_positions[i].second;
		int matches = insert_start - hap_to_ref_count;
		if(matches <= 0){
			return ret;
		}
		hap_to_ref_cigar += to_string(matches);
		hap_to_ref_cigar += "M";
		hap_to_ref_count += matches;
		count_cigar_ops++;
		int insert = insert_end - insert_start;
		hap_to_ref_cigar += to_string(insert);
		if(insert <= 0){
			return ret;
		}
		hap_to_ref_cigar += "I";
		hap_to_ref_count += insert;
		count_cigar_ops++;
	}
	// Now add last bit of hap
	int matches = endLocation - hap_to_ref_count;
	hap_to_ref_cigar += to_string(matches);
	hap_to_ref_cigar += "M";
	if(matches <= 0){
		return ret;
	}	
	count_cigar_ops++;
	bool print = false;
	if(read_name == "19e39e5c-3ee8-4a37-a129-df3930dfe12f"){
		print = true;
		cout << hap_name << "\t" << hap_to_ref_cigar << endl;
		cout << startLocation << "\t" << hap_start << endl;
	}
	// now iterate base by base
	hap_start += startLocation;
	int hap_count = -1;
	int read_count = 0;
	vector<pair<int, char>> final_cig;
	if(readStart > 0 && read_orientation == '+'){
		final_cig.push_back(make_pair(readStart, 'S'));
		read_count = readStart;
		if(print){
			print_debug(read_name+" Added "+to_string(readStart)+"S");
		}
	} else if (read_sequence.length() - readEnd > 0 && read_orientation == '-'){
		final_cig.push_back(make_pair(read_sequence.length() - readEnd, 'S'));
		read_count = read_sequence.length() - readEnd;
		if(print){
			print_debug(read_name+" Added "+to_string(read_sequence.length() - readEnd)+"S");
		}
	}
	if(read_count > read_sequence.length()){
		// Problem with alignment;
		return ret;
	}
	read_count = 0;
	int hap_length = parse_cigar_string(hap_to_ref_cigar);
	// Convert CIGAR strings to vectors of ops with offsets 
	// Start parsing
	istringstream hap_to_ref_parser(hap_to_ref_cigar);
	istringstream read_to_hap_parser(alt_cigar);
	char hr_op;
	uint32_t hr_ol;
	char op;
	uint32_t ol;
	hap_to_ref_parser >> hr_ol >> hr_op;
	hap_to_ref_count = hr_ol;
	int seen_ops = 1;
	while(read_to_hap_parser >> ol >> op){
		if(hap_count > hap_length){
			// Should not iterate past the end of the hap
			break;
		}
		if(hap_count > hap_to_ref_count){
			// parse another op
			if(hap_to_ref_parser >> hr_ol >> hr_op){
				hap_to_ref_count += hr_ol;
				seen_ops++;
			} else {
				// Should not run out of ops to parse before the end of the alignment
				out_lock.lock();
				cerr << "Not enough CIGAR Ops for 1: " << read_name << endl;
				out_lock.unlock();
				return ret;
			}
		}
		if (op == 'M' || op == '=' || op == 'X'){
			//Parse both
			hap_count += ol;
			read_count += ol;
		} else if(op == 'I'){
			read_count += ol;
		} else if(op == 'D'){
			hap_count += ol;
		}
		// Check if our read to hap op is fully covered by the current hap to ref op 
		if(hap_to_ref_count - hap_count >= 0){
			// Covered can project as needed to our output CIGAR based on hr_op and op values
			pair<int, char> item = get_projected_op(ol, op, hr_op);
			if(item.first > 0){
				final_cig.push_back(item);
			}
		} else {
			// op spans a break between hr_op values. Need to go base by base
			int prev_count = hap_count - ol;
			int diff = hap_to_ref_count - prev_count;
			// Output diff number of operations
			pair<int, char> item = get_projected_op(diff, op, hr_op);
			if(item.first > 0){
				final_cig.push_back(item);
			}
			// Now update the hap_to_ref_count by parsing the next op
			if(seen_ops < count_cigar_ops){
				hap_to_ref_parser >> hr_ol >> hr_op;
				hap_to_ref_count += hr_ol;
				seen_ops++;
			}
			// Compute the remainder of the ol to now be projected to this new op
			int rest = ol - diff;
			// Project rest bases against hap with new hr_op value
			if(rest > 0){
				item = get_projected_op(rest, op, hr_op);
				if(item.first > 0){
					final_cig.push_back(item);
				}
			}
		}
	}
	// Add the end soft clip
	if(read_sequence.length()-readEnd > 0 && read_orientation == '+'){
		if(read_sequence.length()-readEnd > read_sequence.length()){
			return ret;
		}
		final_cig.push_back(make_pair(read_sequence.length()-readEnd, 'S'));
	} else if(readStart > 0 && read_orientation == '-'){
		if(readStart > read_sequence.length()){
			return ret;
		}
		final_cig.push_back(make_pair(readStart, 'S'));
	}
	// Collapse the final_cig vector into a CIGAR string
	pair<string, int> final_cig_data = get_final_cig(final_cig);
	// Output bam record for updated alignment
	int orientation = 0;
	if (read_orientation == '-'){
		orientation = 16;
		read_sequence = dna_reverse_complement(read_sequence);
	}
	if(print){
		cout << read_sequence.length() << "  " << parse_cigar_string(final_cig_data.first);
	}
	write_bam_record(final_cig_data.first, read_name, hap_start, orientation, chrom, final_cig_data.second, read_sequence, outfile, header, t);
	return ret;
}

void write_original(htsFile *outfile, bam_hdr_t *header, std::vector<bam1_t*>& read_alignments, const int t){
	using namespace std;
	for(size_t i = 0; i < read_alignments.size(); i++){
		string name = string(bam_get_qname(read_alignments[i]));
		seen_lock.lock();
		if(seen_full_reads.count(name) > 0){
			seen_lock.unlock();
			continue;
		} else {
			seen_full_reads.insert(name);
			seen_lock.unlock();
		}
		bam_lock.lock();
		int ret = sam_write1(outfile, header, read_alignments[i]);
		bam_lock.unlock();
		if(ret == -1){
			out_lock.lock();
			cerr << t << " Write to bam failed for " << name << endl;
			out_lock.unlock();
		}
	}
}

std::unordered_map<std::string,std::unordered_map<std::string,std::string>> align_all_reads(std::unordered_map<std::string, std::unordered_map<std::string,std::string>>& read_sequences_per_sample, const std::string& ref_sequence, const std::vector<std::tuple<int, int, std::string, int, int, int,std::string>>& haplotypes, const int gap_open, const int gap_extend, const int min_mapq, std::string& output, const std::string& chrom, const int pos, const int max_read_len, htsFile *outfile, bam_hdr_t *header, std::unordered_set<std::string>& seen, const std::unordered_map<std::string, std::vector<bam1_t*>>& alignments_per_read, const int t, const int max_mem, mm_tbuf_t *tbuf){
	using namespace std;
	unordered_map<string,unordered_map<string, string>> return_map;
	int ref_start = pos - max_read_len;
	int ref_end = ref_start + 2*max_read_len;
	if(ref_start < 0){
		ref_start = 0;
	}
	if(ref_end > (int) ref_sequence.length()){
		ref_end = (int) ref_sequence.length() - 1;
	}
	string ref_seq = ref_sequence.substr(ref_start, ref_end - ref_start);
	// Get the ref seq and the rest of the sequences 
	unordered_map<string, string> seq_map;
	seq_map["ref"] = ref_seq;
	unordered_map<string,string> output_map;
	unordered_map<string,int> supporting_counts;
	unordered_map<string, mm_idx_t*> mi_list; 
	unordered_map<string, pair<int,int>> mapping_params;
	int n_seqs = 1;
	for(size_t i =0; i < haplotypes.size(); i++){
		// Align reads to best haplotype and ref. If read aligns better to haplotype than ref, update record
		// Output an alignment record for the haplotype
		int insert_size = get<1>(haplotypes[i]);
		int offset = get<0>(haplotypes[i]);
		string best_hap_sequence = get<2>(haplotypes[i]);
		int hap_ref_start = get<3>(haplotypes[i]);
		int hap_ref_end = get<4>(haplotypes[i]);
		string hap_name = get<6>(haplotypes[i]);
		int left_flank = pos+offset-max_read_len;
		int right_flank = pos+offset;
		int right_flank_length = max_read_len;
		int left_flank_length = max_read_len;
		if(left_flank < 0){
			left_flank = 0;
			left_flank_length = max_read_len + left_flank;
		}
		if((right_flank + max_read_len) > (ref_sequence.length() -1)){
			right_flank_length = right_flank + max_read_len - ref_sequence.length() -1;
		}
		string new_hap_seq = ref_sequence.substr(left_flank,left_flank_length)+best_hap_sequence+ref_sequence.substr(right_flank,right_flank_length);
		string hap_cigar = to_string(left_flank_length)+"="+to_string(insert_size)+"I"+to_string(right_flank_length)+"=";
		//print_debug(to_string(t)+" "+to_string(left_flank_length)+" "+to_string(insert_size)+" "+to_string(right_flank_length));
		if(insert_size < 10) {
			continue;
		}
		int ret = write_bam_record(hap_cigar, "hap_"+hap_name, left_flank, 0, chrom, 3, new_hap_seq, outfile, header, t);
		output_map[hap_name] = chrom+"\t"+to_string(hap_ref_start)+"\t"+to_string(hap_ref_end)+"\t"+hap_name+"\t"+best_hap_sequence+"\t";
		supporting_counts[hap_name] = 0;
		seq_map[hap_name] = new_hap_seq;
		n_seqs++;
		// Generate the mi for each alt hap and store params for re-alignment
		mapping_params[hap_name] = make_pair(left_flank, insert_size);
	}
	// Get and index an array of seqs
	vector<mm_idx_t*> mi_vec;
	mm_idxopt_t iopt_all;
	mm_mapopt_t mopt_all;
	mm_set_opt(0, &iopt_all, &mopt_all);
	mopt_all.cap_kalloc = 100000000;
	mm_set_opt("ava-ont", &iopt_all, &mopt_all);
	mopt_all.best_n = 20;
	for(auto it = seq_map.begin(); it != seq_map.end(); ++it){
		const char* s = it->second.c_str();
		const char* n = it->first.c_str();
		auto start1 = chrono::high_resolution_clock::now();
		mm_idx_t *mi_all = mm_idx_str(10,15,0,16,1, &s, &n);
		mm_mapopt_update(&mopt_all, mi_all);
		mi_vec.push_back(mi_all);
		auto end1 = chrono::high_resolution_clock::now();
		auto duration1 = chrono::duration_cast<chrono::microseconds>(end1 - start1);
		count_lock.lock();
		map_time_count += duration1.count();
		count_lock.unlock();
	}
	// Align each read against the set of seqs
	for(auto it = read_sequences_per_sample.begin(); it != read_sequences_per_sample.end(); ++it){
		string sample = it->first;
		for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2){
			bool print = false;
			if(it2->first == "b7a9e8ec-cb85-45cb-aff4-b42929eee391" || it2->first == "0a1f0897-e6f0-4752-98f0-1163309593cb"){
				cout << it2->first << endl;
				print = true;
			}
			string read_sequence = it2->second;
			if(seen.count(it2->first) > 0 ) {
				continue;
			}
			seen.insert(it2->first);
			vector<bam1_t*> read_alignments = alignments_per_read.at(it2->first);
			bool low_mapq = false;
			for(size_t j = 0; j < read_alignments.size(); j++){
				if (read_alignments[j]->core.qual < min_mapq){
					low_mapq = true;
				}
			}
			if(low_mapq){
				// Write out original alignment 
				//write_original(outfile, header, read_alignments, t);
				continue;
			}
			// Align to all seqs and get the best hit
			if(print){
				cout << n_seqs << endl;
			}
			auto start3 = chrono::high_resolution_clock::now();
			string best_hap = align_minimap_multiple(read_sequence, mi_vec, &mopt_all, tbuf, false, print);
			if(print){
				cout << it2->first << " to " << best_hap << endl;
			}
			auto end3 = chrono::high_resolution_clock::now();
			auto duration3 = chrono::duration_cast<chrono::microseconds>(end3 - start3);
			count_lock.lock();
			align_time_count += duration3.count();
			count_lock.unlock();
			if(best_hap != "ref"){
				return_map[sample][it2->first] = best_hap;
				/*
				// Need to project and output
				int insert_size = mapping_params[best_hap].second;
				int left_flank = mapping_params[best_hap].first;
				pair<int,int> ret = project_alignment(pos - max_read_len, left_flank, left_flank+insert_size-2, alt_results[best_hap].cigar, it2->first, read_sequence, alt_results[best_hap].orientation, chrom, outfile, header, alt_results[best_hap].startLocation, alt_results[best_hap].read_start, alt_results[best_hap].read_end, t);
				if(ret.first < 0){
					write_original(outfile, header,read_alignments,t);
				} else {
					if(supporting_counts[best_hap] != 0){
						output_map[best_hap] += ",";
					}
					supporting_counts[best_hap]++;
					output_map[best_hap] += sample+":"+it2->first+":"+alt_results[best_hap].orientation+":";
					output_map[best_hap] += to_string(ret.first);
					output_map[best_hap] += "-";
					output_map[best_hap] += to_string(ret.second);
				}
				*/
			}
		}
	}
	for(size_t i = 0; i < mi_vec.size(); i++){
		mm_idx_destroy(mi_vec[i]);
	}
	/*
	p = 0;
	for(auto it = seq_map.begin(); it != seq_map.end(); ++it){
		delete [] seqs[p];
		delete [] names[p];
		p++;
	}
	delete [] seqs;
	delete [] names;
	*/
	return return_map;
}

std::vector<std::tuple<int,int,int>> get_insert_position(std::string& cigar, int seq_start, int seq_end, int ref_start, int ref_end, char read_orientation, std::string& seq){
	using namespace std;
	vector<tuple<int,int,int>> larger_inserts;
	istringstream cigar_parser(cigar);
	char op;
	uint32_t ol;
	cigar_parser >> ol >> op;
	int seq_count = 0;
	int ref_count = ref_start;
	if(seq_start > 0 && read_orientation == '+'){
		seq_count = seq_start;
	} else if (seq.length() - seq_end > 0 && read_orientation == '-'){
		seq_count = seq.length() - seq_end;
	}
	while(cigar_parser >> ol >> op){
		if (op == 'M' || op == '=' || op == 'X'){
			//Parse both
			ref_count += ol;
			seq_count += ol;
		} else if(op == 'I'){
			// Check if insert length is long enough to be considered 
			if(ol >= 100){
				// Add to the list of larger_inserts
				larger_inserts.push_back(make_tuple(seq_count,ol,ref_count));
			}
			seq_count += ol;
		} else if(op == 'D'){
			ref_count += ol;
		}
	}
	return larger_inserts;
}

std::vector<std::tuple<int, int, std::string, int, int, int, std::string>> get_consensus_hap(const int start, const std::string& ref_sequence, const std::unordered_map<std::string, std::vector<TsvRecord>>& inserts_per_sample, const int ref_window, std::unordered_map<std::string, std::unordered_map<std::string,std::string>>& read_sequences_per_sample, int gap_open, int gap_extend, int t, const int max_mem, const int window, const std::string& racon_path, mm_tbuf_t *tbuf){
	// Extract the read sequences of any insert supporting reads in this window.
	// Compute a consensus for them, select the largest one
	using namespace std;
	int ref_start = start;
	string hap_name = "";
	vector<tuple<int, int, string, int, int, string, string>> read_haps;
	unordered_map<string,int> hap_to_pos;
	vector<string> haps_to_align; 
	int count = 0;
	for(auto it = inserts_per_sample.begin(); it != inserts_per_sample.end(); ++it){
		string sample = it->first;
		for(size_t i = 0; i < it->second.size(); i++){
			// Compute insertion position on ref seq and splice in insertion 
			int insert_start = it->second[i].ref_start - ref_start+1;
			string read_sequence = read_sequences_per_sample[sample][it->second[i].read_id];
			int read_insert_start = it->second[i].read_start - window+1;
			if(read_insert_start < 0){
				read_insert_start = 0;
			}
			int read_insert_end = it->second[i].read_end+window;
			if(read_insert_end > read_sequence.length()-1){
				read_insert_end = read_sequence.length() -1;
			}
			hap_name = it->second[i].read_id+"_"+to_string(start);
			haps_to_align.push_back(read_sequence.substr(read_insert_start, read_insert_end-read_insert_start));
			read_haps.push_back(make_tuple(insert_start, it->second[i].seq.length(), it->second[i].seq, it->second[i].ref_start, it->second[i].ref_end, it->second[i].read_id, read_sequence.substr(read_insert_start, read_insert_end-read_insert_start)));
			hap_to_pos[it->second[i].read_id] = count;
			count++;
		}
	}
	// Pairwise alignment of sequences
	//unordered_map<string, unordered_set<string>> overlap_set;
	//for(size_t i = 0; i < read_haps.size(); i++){
	//	int pos_i = get<0>(read_haps[i]);
	//	for(size_t j = 0; j < read_haps.size(); j++){
	//		int pos_j = get<0>(read_haps[j]);
	//		if(abs(pos_i - pos_j) < 10){
	//			overlap_set[get<5>(read_haps[i])].insert(get<5>(read_haps[j]));
	//			overlap_set[get<5>(read_haps[j])].insert(get<5>(read_haps[i]));
	//		}
	//	}
	//}
	// generate concensus sequence
	int n_seqs = haps_to_align.size();
	sort(haps_to_align.begin(), haps_to_align.end(), []
    (const std::string& first, const std::string& second){
        return first.size() > second.size();
    });
	string hap_sequence = "";
	vector<tuple<int, int, string, int, int, int,string>> corrected_seqs;
	if(n_seqs > 1){
		/*
		auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
		spoa::Graph graph{};
		for (const auto& it : haps_to_align) {
			auto alignment = alignment_engine->Align(it, graph);
			graph.AddAlignment(alignment, it);
		}
		auto consensus = graph.GenerateConsensus();
		hap_sequence = (string) consensus;
		*/
		// initialize variables
		abpoa_t *ab = abpoa_init();
		abpoa_para_t *abpt = abpoa_init_para();
		abpt->out_msa = 0; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
		abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
		abpt->w = 6, abpt->k = 15; abpt->min_w = 10; // minimizer-based seeding and partition
		//abpt->progressive_poa = 1;
		abpt->disable_seeding = 0;
		abpoa_post_set_para(abpt);
		if(hap_name == "0d2c1c54-8f9b-4b53-a4c2-e4c9a6f979b4_25791350"){
			ofstream MyFile("test_failure_reads.fa");
			for(size_t i = 0; i < read_haps.size(); i++){
				string seq_i = get<6>(read_haps[i]);
				MyFile << ">" << get<5>(read_haps[i]) << "\n" << seq_i << "\n";
			}
			MyFile.close();
		}
		// Get the lengths of each hap sequence and convert from ACGT to 0123
		//cout << n_seqs << endl;
		int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
		uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
	    for(size_t i = 0; i < haps_to_align.size(); i++){
			//cout << seq_i.length() << endl;
			//MyFile << ">" << get<5>(read_haps[i]) << "\n" << seq_i << "\n";
	        seq_lens[i] = haps_to_align[i].length();
	        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
	        for (size_t j = 0; j < seq_lens[i]; ++j){
	            bseqs[i][j] = _nt4_table[(int)haps_to_align[i][j]];
	        }
	    }
	    //MyFile.close();
	    //abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	    // 2. variables to store result
	    uint8_t **cons_seq; int **cons_cov, *cons_l, cons_n=0;
	    uint8_t **msa_seq; int msa_l=0;
	    // perform abpoa-msa
	    ab->abs->n_seq = 0; // To re-use ab, n_seq needs to be set as 0
	    int val = abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);
	    if(val == 1){
	    	// MSA failed
	    	cout << "Failure in ABPOA" << endl;
	    	if (cons_n) {
		        for (size_t i = 0; i < cons_n; ++i) {
		            free(cons_seq[i]);
		            free(cons_cov[i]);
		        } 
		        free(cons_seq);
		        free(cons_cov);
		        free(cons_l);
		    }
		    for (size_t i = 0; i < n_seqs; ++i){
		    	free(bseqs[i]);
		    }
		    free(bseqs);
		    free(seq_lens);
			abpoa_free(ab);
			abpoa_free_para(abpt);
			return corrected_seqs;
	    }
	    if(cons_n > 0){
	    	string tmp(cons_l[0], 'N');
	    	int i = 0;
	        for (size_t j = 0; j < cons_l[0]; ++j){
	            tmp[i] = (char)_char256_table[cons_seq[0][j]];
	            i++;
	        }
	        hap_sequence = tmp;
	    }
		if (cons_n) {
	        for (size_t i = 0; i < cons_n; ++i) {
	            free(cons_seq[i]);
	            free(cons_cov[i]);
	        } 
	        free(cons_seq);
	        free(cons_cov);
	        free(cons_l);
	    }
	    for (size_t i = 0; i < n_seqs; ++i){
	    	free(bseqs[i]);
	    }
	    free(bseqs);
	    free(seq_lens);
		abpoa_free(ab);
		abpoa_free_para(abpt);
	} else {
		// Only have one sequence
		for(size_t i = 0; i < read_haps.size(); i++){
			hap_sequence = get<6>(read_haps[i]);
		}
	}
	// Align hap_seq to ref with minimap2 and parse the CIGAR to get the position of the insert on the hap and ref
	// Extract the ref pos as the position where the corrected hap seq should go and add it to the ref there and return
	string ref_seq = ref_sequence.substr(start-window, 2*window);
	if(hap_name == "0d2c1c54-8f9b-4b53-a4c2-e4c9a6f979b4_25791350"){
		ofstream MyFile("test_failure_ref.fa");
		MyFile << ">ref_seq\n" << ref_seq << "\n";
		MyFile.close();
	}
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mm_set_opt(0, &iopt, &mopt);
	const char* seq = ref_seq.c_str();
	mm_idx_t *mi = mm_idx_str(10,15,0,16,1, &seq, NULL);
	mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
	mopt.cap_kalloc = 100000000;
	mm_set_opt("map-ont", &iopt, &mopt);
	if(hap_name == "0d2c1c54-8f9b-4b53-a4c2-e4c9a6f979b4_25791350"){
		ofstream MyFile("test_con.fa");
		MyFile << ">hap_seq\n" << hap_sequence << "\n";
		MyFile.close();
	}
	AlignmentResult complete_result = align_minimap(hap_sequence, mi, &mopt, "", "", "", true, tbuf);
	if(hap_name == "0d2c1c54-8f9b-4b53-a4c2-e4c9a6f979b4_25791350"){
		cout << complete_result.cigar << endl;
	}
	vector<tuple<int,int,int>> positions_on_ref = get_insert_position(complete_result.cigar, complete_result.read_start, complete_result.read_end, complete_result.startLocation, complete_result.endLocation, complete_result.orientation, hap_sequence);
	if(positions_on_ref.size() > 1 || positions_on_ref.size() == 0){
		// More than one insert detected when we compute the concensus 
		return corrected_seqs;
	}
	// Only have one seq
	tuple<int,int,int> concensus_pos = positions_on_ref[0];
	int insert_start = (start - window + get<2>(concensus_pos)) - ref_start+1;
	string insert_seq = hap_sequence.substr(get<0>(concensus_pos), get<1>(concensus_pos));
	if(hap_name == "0d2c1c54-8f9b-4b53-a4c2-e4c9a6f979b4_25791350"){
		cout << hap_name << " " << insert_start << " " << insert_seq.length() << endl;
	}
	corrected_seqs.push_back(make_tuple(insert_start, insert_seq.length(), insert_seq, start - window + get<2>(concensus_pos), start - window + get<2>(concensus_pos) + 1, window, "consensus_"+hap_name));
	return corrected_seqs;
}

void realign_reads(const int t, const std::unordered_map<std::string, faidx_t*>& samples_to_faidx, const std::unordered_map<std::string, std::string>& ref_sequences, htsFile *outfile, bam_hdr_t *header, const std::string& tsv_output){
	using namespace std;
	mm_tbuf_t *tbuf = mm_tbuf_init();
	while(true){
		queue_lock.lock();
		if(alt_hap_queue.empty()){
			queue_lock.unlock();
			mm_tbuf_destroy(tbuf);
			print_debug("Done Thread: "+to_string(t));
			return;
		}
		pair<string, string> read_to_project = alt_hap_queue.front();
		alt_hap_queue.pop();
		queue_lock.unlock();
		string sample = get<0>(read_to_project);
		string read_name = get<1>(read_to_project);
		unordered_map<string, vector<pair<int,string>>> positions_for_chrom;
		int total_positions = 0;
		unordered_map<string, string> val = alt_hap_per_read_per_sample[sample+":"+read_name];
		for(auto it = val.begin(); it != val.end(); it++){
			string chrom; 
			string position;
			// -- Use a data struct to store as seperate fields. 
			stringstream  combined(it->first);
			getline(combined, chrom, ':');
			getline(combined, position, ':');
			int pos = stoi(position);
			string alt_hap_name = it->second;
			positions_for_chrom[chrom].push_back(make_pair(pos, alt_hap_name));
			total_positions++;
		}
		print_debug(to_string(t)+" "+sample+" "+read_name+" "+to_string(positions_for_chrom.size())+" "+to_string(total_positions));
		auto read_index = samples_to_faidx.at(sample);
		int size = 0;
		faidx_lock.lock();
		char* tmp_read = fai_fetch(read_index, read_name.c_str(), &size);
		faidx_lock.unlock();
		string read_sequence = string(tmp_read);
		free(tmp_read);
		for(auto it = positions_for_chrom.begin(); it != positions_for_chrom.end(); it++){
			// For each chrom sort the positions and then apply them to the reference genome to create an alt hap specific to this reads insertions
			vector<pair<int,string>> insert_positions = it->second;
			if(insert_positions.size() > 1){
				sort(insert_positions.begin(), insert_positions.end(), [](const pair<int,string> &left, const pair<int,string> &right) {return left.first < right.first;});
			}
			int window_to_use = 0;
			for(size_t i = 0; i < insert_positions.size(); i++){
				int pos = insert_positions[i].first;
				string hap_name = insert_positions[i].second;
				tuple<int, int, string, int, int, int, string> hap_tuple = consensus_tuples_per_pos[it->first+":"+to_string(pos)][hap_name];
				int hap_window = get<5>(hap_tuple);
				if(hap_window > window_to_use){
					window_to_use = hap_window;
				}
			}
			string ref_sequence = ref_sequences.at(it->first);
			int min_pos = insert_positions[0].first;
			min_pos = min_pos - (read_sequence.length()+window_to_use);
			if(min_pos < 0){
				min_pos = 0;
			}
			int max_pos = insert_positions[insert_positions.size()-1].first;
			max_pos = max_pos + read_sequence.length()+window_to_use;
			if(max_pos > ref_sequence.length()-1){
				max_pos = ref_sequence.length()-1;
			}
			string alt_hap_seq = ref_sequence.substr(min_pos, max_pos-min_pos);
			// Add hapsequences in one by one to base alt hap seq
			int added_bases = 0;
			vector<pair<int,int>> added_insert_positions;
			string hap_names = "";
			bool print = false;
			if(read_name == "b7a9e8ec-cb85-45cb-aff4-b42929eee391" || read_name == "0a1f0897-e6f0-4752-98f0-1163309593cb"){
				print = true;
				cout << read_name << endl;
			}
			for(size_t i = 0; i < insert_positions.size(); i++){
				int pos = insert_positions[i].first;
				string hap_name = insert_positions[i].second;
				if(hap_names != ""){
					hap_names = hap_names+","+hap_name;
				} else {
					hap_names = hap_name;
				}
				tuple<int, int, string, int, int, int, string> hap_tuple = consensus_tuples_per_pos[it->first+":"+to_string(pos)][hap_name];
				int hap_position = (pos + get<0>(hap_tuple)) - (min_pos - added_bases);
				if(print){
					cout << "Added " << hap_name << endl;
				}
				alt_hap_seq = alt_hap_seq.substr(0, hap_position)+get<2>(hap_tuple)+alt_hap_seq.substr(hap_position);
				added_bases += get<1>(hap_tuple);
				added_insert_positions.push_back(make_pair(hap_position, hap_position+get<1>(hap_tuple)));
			}
			// Align the read to this larger alt hap and then project the alignment to generate a new record.
			mm_idxopt_t iopt;
			mm_mapopt_t mopt;
			mm_set_opt(0, &iopt, &mopt);
			const char* seq = alt_hap_seq.c_str();
			mm_idx_t *mi = mm_idx_str(10,15,0,16,1, &seq, NULL);
			mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
			mopt.cap_kalloc = 100000000;
			mm_set_opt("ava_ont", &iopt, &mopt);
			AlignmentResult complete_result = align_minimap(read_sequence, mi, &mopt, "", "", "", true, tbuf);
			//if(read_name == "b7a9e8ec-cb85-45cb-aff4-b42929eee391"){
			//	// write the hap to file 
			//	ofstream seq_output("hap_seq.fa");
			//	seq_output << ">" << hap_names << "\n" << alt_hap_seq << "\n";
			//	seq_output.close();
			//}
			vector<pair<int,int>> ret = project_alignment(min_pos, added_insert_positions, complete_result.cigar, hap_names, read_name, read_sequence, complete_result.orientation, it->first, outfile, header, complete_result.startLocation, complete_result.endLocation, complete_result.read_start, complete_result.read_end, t);
			// Output as needed here
			mm_idx_destroy(mi);
		}
	}
}


void realign_chrom(const int t, const std::unordered_map<std::string, std::map<int, std::vector<TsvRecord>>>& all_tsv_records, const std::unordered_map<std::string, faidx_t*>& samples_to_faidx, const std::unordered_map<std::string, std::string>& samples_to_bam, const std::unordered_map<std::string, std::string>& ref_sequences, const int gap_open, const int gap_extend, const int min_mapq, const int depth_filter, htsFile *outfile, bam_hdr_t *header, const std::string& tsv_output, const int max_mem, const std::string& racon_path){
	using namespace std;
	mm_tbuf_t *tbuf = mm_tbuf_init();
	while(true){
		queue_lock.lock();
		if(positions_for_threads.empty()){
			queue_lock.unlock();
			mm_tbuf_destroy(tbuf);
			print_debug("Done Thread: "+to_string(t));
			return;
		}
		tuple<string,int,int> q_pos = positions_for_threads.front();
		positions_for_threads.pop();
		queue_lock.unlock();
		string chrom = get<0>(q_pos);
		int pos = get<1>(q_pos);
		int window = get<2>(q_pos);
		string output = ""; 
		string ref_sequence = ref_sequences.at(chrom);
		unordered_map<string, vector<TsvRecord>> inserts_per_sample;
		unordered_map<string, vector<bam1_t*>> alignments_per_read;
		unordered_map<string, unordered_map<string,string>> read_sequences_per_sample;
		unordered_set<string> seen;
		size_t max_read_len = 0;
		// Get all the reads in the window, insert set and all reads
		vector<pair<bam1_t*,string>> reads_in_window = get_reads_in_window(samples_to_bam, pos, window, chrom);
		if((int) reads_in_window.size() > depth_filter && depth_filter > 0){
			// More reads than our filter, skip this window
			for(size_t j = 0; j < reads_in_window.size(); j++){
				vector<bam1_t*> tmp; 
				tmp.push_back(reads_in_window[j].first);
				write_original(outfile, header,tmp,t);
				bam_destroy1(reads_in_window[j].first);
			}
			continue;
		}
		get_windows(pos, chrom, window, all_tsv_records, inserts_per_sample, alignments_per_read, read_sequences_per_sample, max_read_len, samples_to_faidx, reads_in_window);
		// Get consensus alternative haplotype sequence for re-alignment

		// -- Make a parameters struct
		// -- Return a struct 
		auto start1 = chrono::high_resolution_clock::now();
		vector<tuple<int, int, string, int, int, int, string>> consensus_tuples = get_consensus_hap(pos, ref_sequence, inserts_per_sample, max_read_len, read_sequences_per_sample, gap_open, gap_extend, t, max_mem, window, racon_path, tbuf);
		// Get a set of alternative haplotypes for the window
		if(consensus_tuples.size() == 0){
			unordered_map<string, tuple<int, int, string, int, int, int, string>> hap_set = get_hap_set(pos, ref_sequence, inserts_per_sample, max_read_len, window, t);
			// Take hap sets and align it to each read with an insertion, Select best haplotype
			string best_hap= get_best_hap(inserts_per_sample, hap_set, read_sequences_per_sample, gap_open, gap_extend, t, max_mem, false, tbuf);
			consensus_tuples.push_back(hap_set[best_hap]);
		}
		auto end1 = chrono::high_resolution_clock::now();
		auto duration1 = chrono::duration_cast<chrono::microseconds>(end1 - start1);
		count_lock.lock();
		hap_time_count += duration1.count();
		count_lock.unlock();
		print_debug(to_string(t)+" "+chrom+":"+to_string(pos)+" "+to_string(reads_in_window.size())+" "+to_string(max_read_len)+" "+to_string(window)+" "+to_string(consensus_tuples.size()));
		unordered_map<string,unordered_map<string, string>> alt_hap_map = align_all_reads(read_sequences_per_sample, ref_sequence, consensus_tuples, gap_open, gap_extend, min_mapq, output, chrom, pos, max_read_len, outfile, header, seen, alignments_per_read, t, max_mem, tbuf);
		//Update global list of reads per sample to map each read to its respective alternative haplotype
		hap_lock.lock();
		for(auto it = alt_hap_map.begin(); it !=alt_hap_map.end(); it++){
			for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++){
				alt_hap_per_read_per_sample[it->first+":"+it2->first][chrom+":"+to_string(pos)] = it2->second;
			}
		}
		for(size_t i = 0; i < consensus_tuples.size(); i++){
			consensus_tuples_per_pos[chrom+":"+to_string(pos)][get<6>(consensus_tuples[i])] = consensus_tuples[i];
		}
		hap_lock.unlock();
		for(size_t j = 0; j < reads_in_window.size(); j++){
			bam_destroy1(reads_in_window[j].first);
		}
		// Clear all maps/sets
		for(auto it = inserts_per_sample.begin(); it != inserts_per_sample.end(); it++){
			it->second.erase(it->second.begin(), it->second.end());
		}
		inserts_per_sample.erase(inserts_per_sample.begin(), inserts_per_sample.end());
		for(auto it = read_sequences_per_sample.begin(); it != read_sequences_per_sample.end(); it++){
			it->second.erase(it->second.begin(), it->second.end());
		}
		read_sequences_per_sample.erase(read_sequences_per_sample.begin(), read_sequences_per_sample.end());
		/*
		ofstream output_file;
		tsv_lock.lock();
		output_file.open(tsv_output, ios_base::app);
		for(auto it = output_map.begin(); it != output_map.end(); it++){
			output_file << it->second;
		}
		output_file.close();
		tsv_lock.unlock();
		*/
	}
}	

// Main wrapper function that gets called from somrit interface
void ReAlign::realign_all(const std::string& bam_list, const std::string& tsv_list, const std::string& fastq_list, const std::string& sample_list, const std::string& output_tsv, const std::string& output_bam, const int window, const std::string& ref, const int gap_open, const int gap_extend, const int min_mapq, const std::string& chrom_list, const int pass_only, const int depth_filter, const int threads, const int max_mem, const std::string& racon_path){
	using namespace std;
	// Read in the lists of bams, fastqs and tsvs
	htsFile *outfile= hts_open(output_bam.c_str(), "w");
	ofstream tsv_output;
	tsv_output.open(output_tsv);
	tsv_output << "Chromosome\tStart\tEnd\tAltHaplotypeRead\tInsertSequence\tSupportingReads\n";
	tsv_output.close();
	vector<string> samples;
	stringstream sss(sample_list);
	while( sss.good() ){
		string substr;
		getline( sss, substr, ',' );
		samples.push_back( substr );
	}
	const unordered_map<string, string> samples_to_bam = populate_dict(bam_list, samples);
	const unordered_map<string, string> samples_to_fastq = populate_dict(fastq_list, samples);
	const unordered_map<string, string> samples_to_tsv = populate_dict(tsv_list, samples);
	unordered_map<string, faidx_t*> samples_to_faidx;
	for(auto it = samples_to_fastq.begin(); it != samples_to_fastq.end(); it++){
		samples_to_faidx[it->first] = fai_load(it->second.c_str());
	}
	unordered_set<string> chromosomes;
	stringstream scl(chrom_list);
	while( scl.good() ){
		string substr;
		getline( scl, substr, ',' );
		chromosomes.insert( substr );
	}
	// Open each tsv and read in data
	unordered_map<string, map<int, vector<TsvRecord>>> all_tsv_records;
	unordered_map<string, pair<int,int>> start_end_positions_per_chrom;
	bam_hdr_t *header = NULL;
	bool filter_pass = false;
	if(pass_only > 0){
		filter_pass = true;
	}
	print_debug("Get all positions");
	get_all_tsv_positions(samples_to_tsv, samples_to_bam, header, all_tsv_records, start_end_positions_per_chrom, outfile, filter_pass);
	// Iterate through each window and get a vector of all tsv_records in that window
	// -- Use a type def
	const unordered_map<string, set<pair<int,int>>> positions = get_positions(all_tsv_records, start_end_positions_per_chrom, window);
	// Setup faidx for reference
	auto ref_index = fai_load(ref.c_str());
	// For each window now that we have all reads, get the haplotypes, then align haplotypes to insert supporting reads
	unordered_map<string, string> ref_sequences;
	print_debug("Get all ref sequences");
	for(auto it = all_tsv_records.begin(); it != all_tsv_records.end(); ++it){
		string chrom = it->first;
		if(chromosomes.count(chrom) == 0){
			//cout << "Cannot find " << chrom << endl;
			continue;
		}
		int size = 0;
		char* tmp_ref = fai_fetch(ref_index,chrom.c_str(), &size);
		string ref_sequence = string(tmp_ref);
		ref_sequences[chrom] = ref_sequence;
		free(tmp_ref);
		// Setup a map of position per thread
		set<pair<int,int>> positions_for_chrom = positions.at(chrom);
		//cout << chrom << " " << positions_for_chrom.size() << endl;
		for(auto pos = positions_for_chrom.begin(); pos != positions_for_chrom.end(); pos++){
			positions_for_threads.push(make_tuple(chrom,pos->first,pos->second));
		}
	}
	fai_destroy(ref_index);
	print_debug("Identify and Align reads to Alternative Haplotypes");
	vector<thread> thread_pool;
	for (int i = 0; i < threads; i++) {
		thread_pool.push_back(thread(realign_chrom, i, cref(all_tsv_records), cref(samples_to_faidx), cref(samples_to_bam), cref(ref_sequences), gap_open, gap_extend, min_mapq, (int) depth_filter*samples.size(), outfile, header, cref(output_tsv), max_mem, racon_path));
	}
	for (auto &th : thread_pool) {
		th.join();
	}
	print_debug("ReAlign reads to Alternative Haplotypes");
	// For each read in alt_hap_per_read_per_sample, re-algin it
	for(auto it = alt_hap_per_read_per_sample.begin(); it != alt_hap_per_read_per_sample.end(); it++){
		string sample; 
		string read;
		stringstream  combined(it->first);
		getline(combined, sample, ':');
		getline(combined, read, ':');
		alt_hap_queue.push(make_pair(sample, read));
	}
	vector<thread> thread_pool2;
	for (int i = 0; i < threads; i++) {
		thread_pool2.push_back(thread(realign_reads, i, cref(samples_to_faidx), cref(ref_sequences), outfile, header, cref(output_tsv)));
	}
	for (auto &th : thread_pool2) {
		th.join();
	}
	for(auto it = samples_to_faidx.begin(); it != samples_to_faidx.end(); it++){
		fai_destroy(it->second);;
	}
	sam_hdr_destroy(header);
	hts_close(outfile);
	print_debug(to_string(hap_time_count)+" "+to_string(map_time_count)+" "+to_string(align_time_count));
}
// close extern C
}


int main(){
	std::cout << "TEST" << std::endl;
	ReAlign* tmp = new ReAlign();
	tmp->realign_all("/.mounts/labs/simpsonlab/users/adcosta/projects/HPRC_Validation/validation/HG00438/Sample1/mapped/Sample1.sorted.bam", "/.mounts/labs/simpsonlab/users/adcosta/projects/HPRC_Validation/validation/HG00438/Sample1/read_analysis/Sample1.all.read_insertions.tsv", "/.mounts/labs/simpsonlab/users/adcosta/projects/HPRC_Validation/validation/HG00438/Sample1/fastq/Sample1.fastq.gz", "Sample1", "/.mounts/labs/simpsonlab/users/adcosta/projects/HPRC_Validation/validation/HG00438/Sample1.test.tsv", "/.mounts/labs/simpsonlab/users/adcosta/projects/HPRC_Validation/validation/HG00438/Sample1.test.bam", 1000, "/.mounts/labs/simpsonlab/data/references/GRCh38_no_alt_analysis_set.GCA_000001405.15.fna", 11,1, 10, "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY", 1, 1000, 20, 5000000000, "/u/adcosta/SimpsonLab/racon/build/bin/racon");
	delete tmp;
	std::cout << "END TEST" << std::endl;
	return 0;
}



