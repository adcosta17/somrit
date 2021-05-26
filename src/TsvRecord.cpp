#include "TsvRecord.hpp"

TsvRecord::TsvRecord(std::string i, int rs, int re, std::string c, int fs, int fe, int m, std::string sq, std::string sp) :
		read_id(i), read_start(rs),
		read_end(re), chrom(c), ref_start(fs),
		ref_end(fe), mapped(m), seq(sq), 
		sample(sp) {}


