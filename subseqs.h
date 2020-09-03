#ifndef _SUBSEQS_H
#define _SUBSEQS_H

#include <vector>
#include <string>
#include "chain.h"






struct Subseq_context{
	vector <std::string> data_file_names;
	vector <std::string> res_file_names;
	vector <uint32_t> aligned_idx;
	vector <uint32_t> seg_len;

};



struct Pkg_context{

	vector <std::string> files_names;
	vector <size_t> files_size;
};

// template<typename Context>
struct Thread_p{
	unsigned thread_id;
	Pkg_context pkg_context;
	// Context context;
};


void thread_seq_align(Pkg_context context, uint16_t thread_id);
void *pool_worker(void *p);
void launch_thread_pool(Pkg_context &context, unsigned threads);




void write_raw(const char *ptr, size_t count, FILE *_f);
void writesubseqs(_chain *chain);
void subseqfile(int fileID, uint32_t *start, uint32_t *end);



std::string extract_dir( std::string s);
void set_path();
std::string set_filename(uint32_t fileID);
// void sub_file(Subseq_context *context, uint32_t fileID, uint32_t range_first, 
// 	uint32_t range_second, _chain *chain);
void sub_file(Subseq_context *context, uint32_t fileID, uint32_t range_first, 
	uint32_t range_second, std::vector<_chain> chain);
// void write_sub_file(_chain *chain);
char *getline(FILE *_fp );
void write_sub_file( std::vector<_chain> chain);


#endif