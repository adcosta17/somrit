//---------------------------------------------------------
// Copyright 2023 Ontario Institute for Cancer Research
// Written by Alister D'Costa (alister.d'costa@oicr.on.ca)
//---------------------------------------------------------
//
// ReAlign.hpp

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
	void realign_all(const std::string& bam_list, const std::string& tsv_list, const std::string& fastq_list, const std::string& sample_list, const std::string& output_tsv, const std::string& output_bam, const int window, const std::string& ref, const int gap_open, const int gap_extend, const int min_mapq, const std::string& chrom_list, const int pass_only, const int depth_filter, const int threads, const int max_mem, const int high_mem, const int include_haps, const int max_window_size, const int max_insert_size, const int max_con);
	void haps_only(const std::string& bam_list, const std::string& tsv_list, const std::string& fastq_list, const std::string& sample_list, const std::string& output_fasta, const int window, const std::string& ref, const int gap_open, const int gap_extend, const int min_mapq, const std::string& chrom_list, const int pass_only, const int depth_filter, const int threads, const int max_mem, const int high_mem, const int include_haps, const int max_window_size, const int max_insert_size, const int max_con);
};

extern "C" {
	ReAlign* ReAlign_new(){ return new ReAlign(); }
	void ReAlign_realign_all(ReAlign* r, const char* bam_list, const char* tsv_list, const char* fastq_list, const char* sample_list, const char* output_tsv, const char* output_bam, const int window, const char* ref, const int gap_open, const int gap_extend, const int min_mapq, const char* chrom_list, const int pass_only, const int depth_filter, const int threads, const int max_mem, const int high_mem, const int include_haps, const int max_window_size, const int max_insert_size, const int max_con){ 
		std::cerr << " realign " << std::endl;
		std::string bl(bam_list);
		std::string tl(tsv_list);
		std::string fl(fastq_list);
		std::string sl(sample_list);
		std::string ot(output_tsv);
		std::string ob(output_bam);
		std::string rs(ref);
		std::string cl(chrom_list);
		r->realign_all(bl,tl,fl,sl,ot,ob,window,rs,gap_open,gap_extend,min_mapq,cl,pass_only,depth_filter,threads,max_mem,high_mem,include_haps,max_window_size,max_insert_size, max_con);
		std::cerr << " realign completed" << std::endl;
	}
	void ReAlign_haps_only(ReAlign* r, const char* bam_list, const char* tsv_list, const char* fastq_list, const char* sample_list, const char* output_fa, const int window, const char* ref, const int gap_open, const int gap_extend, const int min_mapq, const char* chrom_list, const int pass_only, const int depth_filter, const int threads, const int max_mem, const int high_mem, const int include_haps, const int max_window_size, const int max_insert_size, const int max_con){ 
		std::cerr << " haps_only " << std::endl;
		std::string bl(bam_list);
		std::string tl(tsv_list);
		std::string fl(fastq_list);
		std::string sl(sample_list);
		std::string of(output_fa);
		std::string rs(ref);
		std::string cl(chrom_list);
		r->haps_only(bl,tl,fl,sl,of,window,rs,gap_open,gap_extend,min_mapq,cl,pass_only,depth_filter,threads,max_mem,high_mem,include_haps,max_window_size,max_insert_size, max_con);
		std::cerr << " haps_only completed" << std::endl;
	}
}
