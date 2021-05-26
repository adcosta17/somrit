#pragma once
#include <string>

class TsvRecord 
{
public:
	std::string read_id;
	int read_start;
	int read_end;
	std::string chrom;
	int ref_start;
	int ref_end;
	int mapped;
	std::string seq;
	std::string sample;

	TsvRecord(std::string i, int rs, int re, std::string c, int fs, int fe, int m, std::string sq, std::string sp);

	};