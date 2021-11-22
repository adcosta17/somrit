#pragma once
#include <unordered_map>
#include <map>
#include <set>
#include <unordered_set>
#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <tuple>

class ReAlign 
{
	public:
	void realign_all(const std::string& bam_list, const std::string& tsv_list, const std::string& fastq_list, const std::string& sample_list, const std::string& output_tsv, const std::string& output_bam, const int window, const std::string& ref, const int gap_open, const int gap_extend, const int min_mapq, const std::string& chrom_list, const int pass_only, const int depth_filter, const int threads, const int max_mem, const std::string& racon_path);
};

extern "C" {
	ReAlign* ReAlign_new(){ return new ReAlign(); }
	void ReAlign_realign_all(ReAlign* r, const char* bam_list, const char* tsv_list, const char* fastq_list, const char* sample_list, const char* output_tsv, const char* output_bam, const int window, const char* ref, const int gap_open, const int gap_extend, const int min_mapq, const char* chrom_list, const int pass_only, const int depth_filter, const int threads, const int max_mem, const char* racon_path){ 
		std::cout << " realign " << std::endl;
		std::string bl(bam_list);
		std::string tl(tsv_list);
		std::string fl(fastq_list);
		std::string sl(sample_list);
		std::string ot(output_tsv);
		std::string ob(output_bam);
		std::string rs(ref);
		std::string cl(chrom_list);
		std::string rp(racon_path);
		r->realign_all(bl,tl,fl,sl,ot,ob,window,rs,gap_open,gap_extend,min_mapq,cl,pass_only,depth_filter,threads,max_mem,rp);
		std::cout << " realign completed" << std::endl;
	}
}
