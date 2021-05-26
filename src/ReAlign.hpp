#pragma once
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <tuple>

class ParseSVs 
{
	public:
	static std::unordered_map<std::string, std::vector<std::tuple<int,int,std::string,int,int>>> parse_inserts_tsv(std::string& read_positions_file, int window);
	static std::unordered_map<int, std::tuple<std::string, int, std::string, int, std::unordered_set<std::string>>> parse_vcf(std::string& read_positions_file, std::string& sv_type_filter);
	static std::unordered_map<int, std::tuple<std::string, int, std::string, int, std::unordered_set<std::string>>> parse_bedpe(std::string& read_positions_file, std::string& sv_type_filter);
	static std::unordered_map<int, std::tuple<std::string, int, std::string, int, std::unordered_map<std::string, std::pair<int,int>>>> get_read_bam_positions(std::unordered_map<int, std::tuple<std::string, int, std::string, int, std::unordered_set<std::string>>>& reads_to_check, std::string& input_bam);
};