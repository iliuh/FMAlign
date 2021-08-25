#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <pthread.h>
#include <thread>
#include <mutex> 
#include <string.h>
#include "common_seed.h"
#include "FMtree/bwt.h"
#include "load_seqs.h"
#include "Auxiliary.h"
#include "subseqs.h"

using std::vector;
using std::string;


const char dir_separator = '/';
std::string tmpdir;
std::mutex mtx;

// int threads;

// ***********************thread***********************

void thread_seq_align(Pkg_context context, uint16_t thread_id)
{
	
	std::string cmnd;
	// pkg = MAFFT;//MAFFT
	switch (pkg)
	{
		case HALIGN:
			cmnd.append("java -cp ")
				.append("packages/HAlign:")
				.append("packages/HAlign/HAlign2.1.jar halignWrapper \"");
			break;
		case MAFFT:
			cmnd.append("packages/MAFFT/mafftWrapper \"");
			break;
	}

	// std::cout<<context.files_names[thread_id]<<std::endl;
	cmnd.append(context.files_names[thread_id].c_str()).append("\" ")
		.append(std::to_string(context.files_size[thread_id])).append(" ")
		.append(std::to_string(thread_id))
		.append(" ").append(std::to_string(0))
		.append(" ").append(std::to_string(1));

	cmnd.append(" > /dev/null");
	system(cmnd.c_str());
}


void *pool_worker(void *p)
{
	//这里就调用seq_align，seq_align的两个参数来自P,context是pkgcontext
    struct Thread_p *_context;
    _context = (struct Thread_p *)p;

    thread_seq_align(_context->pkg_context, _context->thread_id );
	pthread_exit(NULL);

}

void launch_thread_pool(Pkg_context &context, unsigned threads)
{

	pthread_t t[threads];
	Thread_p p[threads];


	for(unsigned int i = 0; i < threads; ++i) {

		p[i].thread_id = i;
		p[i].pkg_context = context;

		// fprintf(stdout, "thread id %d\n", i);
		int ret = pthread_create(&t[i], NULL, pool_worker, (void*)&p[i]);
		if (ret !=0)
		{
           fprintf(stderr, "pthread_create error: error_code = %d\n", ret);
           exit (-1);
		}

	}

	for (unsigned int i = 0; i < threads; ++i)
	{
		pthread_join(t[i], NULL); 
	}
	

}

// ***********************thread***********************



void removefile(char * filename)
{
	if ((remove(filename)) != 0)
	{
		fprintf(stderr, "Failed to delete file...\n");
	}
}

void write_raw(const char *ptr, size_t count, FILE *_f)
{
	size_t n;
	if ((n = fwrite((const void*)ptr, 1, count, _f)) != count)
	{
		fprintf(stderr, "write file err\n");
	}
}

std::string extract_dir( std::string s)
{
	//如果没找到dir_separator就是置空，找到就拷贝root
	return s.find_last_of(dir_separator) == std::string::npos ? "" : s.substr(0, s.find_last_of(dir_separator));
}
void set_path()
{

	// char tmpdir == "/dev/shm";
	// char tmpdir[PATH_MAX_LEN];
	// std::string tmpdir;
	char path[PATH_MAX_LEN];
	// char root[PATH_MAX_LEN];
	if (readlink("/proc/self/exe", path, sizeof(path) - 1) == -1)
		throw std::string("readlink() failed");

	std::string root = extract_dir(path);
	if (tmpdir.length() == 0)
	{
		tmpdir = root;
		tmpdir.append("/temp");
	}

	// std::cout<<tmpdir<<std::endl;
	// std::string file_name_;
	// std::stringstream ss;
	// ss.setf(std::ios::hex, std::ios::basefield);
	// if (tmpdir != "")
	// 	ss << tmpdir << dir_separator;
	// // ss << "chain-" << hash_key << "-" << to_string(0) << ".tmp";
	// ss << "subfile"  << "-" << std::to_string(0) << ".tmp";
	// ss >> file_name_;
	// std::cout<<file_name_<<std::endl;

	// if (output_file.length() == 0)
	// {
	// 	output_file = input_ref_file;
	// 	output_file.append(".malign");
	// }
}


std::string set_filename(uint32_t fileID)
{
	timeval count;
	gettimeofday(&count, NULL);
	long long hash_key = count.tv_sec + count.tv_usec + getpid();


	std::string file_name_;
	std::stringstream ss;
	ss.setf(std::ios::hex, std::ios::basefield);
	if (tmpdir != "")
		ss << tmpdir << dir_separator;
	// ss << "chain-" << hash_key << "-" << to_string(0) << ".tmp";
	ss << "subfile-"  << hash_key << "-" << std::to_string(fileID) << ".tmp";
	ss >> file_name_;
	// std::cout<<file_name_<<std::endl;
	ss.str("");
	return file_name_;
}


void sub_file(Subseq_context *context, uint32_t threadID, uint32_t range_first, 
	uint32_t range_second, std::vector<_chain> chain)
{
	mtx.lock();
	std::vector <uint32_t> s_pos(njob,0);//记录子文件中序列开始的位置，然后根据长度得到子序列
	// uint32_t s_pos[njob] = {0};
	uint32_t tmp_range;

	if (range_first > 0)
	{
		_chain chain_now = chain[std::min(range_first, chain_size_) - 1];
		for (uint32_t i = 0; i < njob; ++i)
		{
			s_pos[i] = chain_now.pos[i] + chain_now.wide;
			// fprintf(stderr, "s_pos  %u\n", chain_now.wide);
		}
	}

	for (uint32_t i = range_first; i < range_second; i++)
	{
		std::string result_file;
		std::string data_file;
		result_file = set_filename(i);
		data_file = set_filename(i);

		FILE* _result_fp = fopen(result_file.c_str(), "wb+");

		if (_result_fp == NULL)
		{
			fprintf(stdout, "Failed to make %s! will exit ...\n", result_file.c_str());	
		}

		FILE* _data_fp = fopen( data_file.c_str(), "wb+");
		if (_data_fp == NULL)
		{
			fprintf(stdout, "Failed to make %s! will exit ...\n", data_file.c_str());	
		}


		uint32_t max_len = 0;
		uint32_t len = 0;
		char *str;
		str = (char *)malloc(sizeof(char) * SEQ_MAX_LENGTH);

		for (uint32_t k = 0; k < njob; ++k)//每次写入一行
		{
			len = (i < chain_size_ ? chain[i].pos[k] : seqset[k].seq_size) - s_pos[k];
			if (len > max_len) max_len = len;
			// fprintf(stdout, "..%u\n", len);
			// fprintf(stdout, "%u \n", chain[i].pos[k]);
			// fprintf(stdout, "%u : %u - %u\n", chain[i].pos[k],  seqset[k].seq_size, s_pos[k]);
			memcpy(str, seqset[k].seq + s_pos[k], len);
			str[len] = '\0';

			write_raw(seqset[k].name, seqset[k].name_size, _data_fp);
			write_raw("\n", 1, _data_fp);
		
			write_raw(str, len, _data_fp);
			write_raw("\n", 1, _data_fp);
			// free(str);

			if(i < chain_size_) s_pos[k] = chain[i].pos[k] + chain[i].wide;
		}
		// fprintf(stdout, "\n", len);
		fclose(_result_fp);
		fclose(_data_fp);

		context->seg_len.push_back( max_len ) ;
		context->aligned_idx.push_back( i ) ;

		context->res_file_names.push_back( result_file);
		context->data_file_names.push_back( data_file );
		// free(str);
	}
	
	mtx.unlock();
}


char *getline(FILE *_fp )
{
	char *out_str;
	out_str = (char *)malloc(sizeof(char) * 1);

	uint32_t n = 0;

	char *str;
	str = (char *)malloc(sizeof(char) * PER_READ_LEN);
	
	fread(str, 1, PER_READ_LEN, _fp);
	
	
	const char *p = (const char*)memchr(str, '\n', PER_READ_LEN);
	if (p == 0)
	{
		out_str = (char *)realloc(out_str, sizeof(char) * (n + PER_READ_LEN));
		memcpy(out_str + n, str, PER_READ_LEN);
		// strncpy(out_str + n, str, PER_READ_LEN);
		free(str);

		n += PER_READ_LEN;
		fread(str, 1, PER_READ_LEN, _fp);
	}
	else
	{

		const int num = p - str;



		out_str = (char *)realloc(out_str, sizeof(char) * (n + num + 1));
		memcpy(out_str + n, str, num);
		fprintf(stdout, "%s  %d\n", out_str, num);

		
		free(str);

		n += num; 
	}

	int offset;
	offset = &str[PER_READ_LEN - 1]- p ;
	fseek( _fp, -offset, SEEK_CUR );
	fprintf(stdout, "offset  %d\n", offset);

	
	return out_str;

}
 
void to_lower(char *s){
    int len=strlen(s);
    for(int i=0;i<len;i++){

        s[i]=tolower(s[i]);
    }
}




void write_sub_file( std::vector<_chain> chain)    //If the number of commonseeds is not 0
{

	// threads = std::min(threads, (int)chain_size_ + 1);
	// fprintf(stdout, "Treads %d\n", threads);

	vector<uint16_t> p_range(threads + 1);//计算每个线程应该计算的子文件数量，
	for (size_t i = 1, r = (chain_size_ + 1) % threads; r > 0; i++, r--) p_range[i]++;
	for (size_t i = 1, m = (chain_size_ + 1) / threads; i <= threads; i++) p_range[i] += p_range[i - 1] + m;
	
	vector <Subseq_context> context_(threads); 
	std::thread t[threads];
	for (uint16_t i = 0; i < threads; ++i)
	{

		t[i] = std::thread(sub_file, &context_[i], i, p_range[i], p_range[i+1], chain);
	}
	for (uint16_t i = 0; i < threads; ++i)
	{
		t[i].join();
	}
	




	//需要对subfile的信息进行整理，好进行并行计算
	struct Seg_len_index {
		size_t len, index;
		bool operator < (const Seg_len_index &right)
		{
			return right.len < len;
		}
	};
	vector <string> res_file_names(chain_size_ + 1);//这个就是需要存储n个文件名
	vector <string> data_file_names(chain_size_ + 1);
	vector <uint32_t> aligned_idx(chain_size_ + 1);
	vector <Seg_len_index>  seg_len_index(chain_size_ + 1);

	
	for (size_t i = 0, n = 0; i < threads; ++i)//这里就把线程里的信息全都存储在文件这个文件里面
	{
		Subseq_context *cntxt = &context_[i];//第i个线程的文件
		for (size_t k = 0; k < cntxt->aligned_idx.size(); k++, n++)
		{
			res_file_names[n] = cntxt->res_file_names[k];
			data_file_names[n] = cntxt->data_file_names[k];
			aligned_idx[n] = cntxt->aligned_idx[k];
			seg_len_index[n] = { cntxt->seg_len[k], n };//每个子文件里面最长序列的长度，n就是index
		}
	}
	context_.clear();

	// std::sort(&seg_len_index[0], &seg_len_index[seg_len_index.size() - 1] + 1);
	vector <vector<uint32_t>> files_names_index(threads);
	vector <uint32_t> files_size(threads);//每个线程总文件数
	vector <size_t> tmp_files_len(threads);//每个线程中，序列的总长
	for (size_t i = 0; i < seg_len_index.size(); i++)
	{
		size_t idx = 0;
		for (size_t k = 1; k < threads; k++)
			if (tmp_files_len[k] < tmp_files_len[idx]) idx = k;
		tmp_files_len[idx] += seg_len_index[i].len;

		files_size[idx]++;
		files_names_index[idx].push_back(seg_len_index[i].index);
	}



	vector <std::string> files_names(threads);

	for (uint32_t i = 0; i < threads; ++i)
	{
		std::string tmp_file;
		tmp_file = set_filename(i);

		FILE* _tmp_fp = fopen(tmp_file.c_str(), "wb+");//用于保存子文件和比对结果文件的文件名
		if (_tmp_fp == NULL)
		{
			fprintf(stdout, "Failed to make %s! will exit ...\n", tmp_file.c_str());	
		}

		files_names[i] = tmp_file;

		for (uint32_t k = 0; k < files_size[i]; k++)
		{
			write_raw(data_file_names[files_names_index[i][k]].c_str(), data_file_names[files_names_index[i][k]].size(), _tmp_fp);
			write_raw("\n", 1, _tmp_fp);
			write_raw(res_file_names[files_names_index[i][k]].c_str(), res_file_names[files_names_index[i][k]].size(), _tmp_fp);
			write_raw("\n", 1, _tmp_fp);
		}

		fclose(_tmp_fp);
	}


	//删除vector
	struct timeval start_align, end_align;
	double timeconsumed;
	gettimeofday(&start_align, NULL);

	Pkg_context pkg;
	pkg.files_names.assign(files_names.begin(), files_names.end());
	pkg.files_size.assign(files_size.begin(), files_size.end());
	launch_thread_pool(pkg, threads);

	gettimeofday(&end_align, NULL);
	timeconsumed = end_align.tv_sec-start_align.tv_sec +(end_align.tv_usec-start_align.tv_usec)/1000000.0;
	fprintf(stdout, "Align squences time: %5.3f seconds\n",  timeconsumed);



	//对结果进行拼接

	
	char **chain_seq;
	chain_seq = (char **)malloc(sizeof(char *) * chain_size_);
	for (int i = 0; i < chain_size_; ++i)
	{

		chain_seq[i] = (char *)malloc(sizeof(char) * (chain[i].wide + 2));
		strncpy(chain_seq[i], seqset[0].seq + chain[i].pos[0], chain[i].wide);
		to_lower(chain_seq[i]);
		chain_seq[i][chain[i].wide] = '\0';
		// puts(chain_seq[i]);	
	}


	vector <vector <string>> aligned;
	aligned.resize(chain_size_ + 1);//这是要存下所有比对的字符串吗？
	for (size_t i = 0; i <= chain_size_; i++)
		aligned[i].resize(njob);

	uint32_t aligned_size[njob] ={0};
	uint32_t tmp_size = 0;


	char *tmp_aligned;
	tmp_aligned = (char *)malloc(sizeof(char) * SEQ_MAX_LENGTH);

	uint32_t seq_index = 0;
	uint32_t size = 0;

	for (uint32_t i = 0; i <= chain_size_ ; ++i)
	{
		FILE* _aligned_fp = fopen(res_file_names[i].c_str(), "r");
		if (_aligned_fp == NULL)
		{
			fprintf(stdout, "result Failed to make %s! will exit ...\n", res_file_names[i].c_str());	
		}

		fgets(tmp_aligned ,SEQ_MAX_LENGTH,_aligned_fp);
		for (int j = 0; j < njob; ++j)
		{
			fgets(tmp_aligned ,SEQ_MAX_LENGTH,_aligned_fp);

			while(!feof(_aligned_fp) && tmp_aligned[0] != '>')
			{
				tmp_size = strlen(tmp_aligned);
				tmp_aligned[tmp_size - 1] = '\0';
				// fprintf(stdout, "row %u \n", tmp_size);
				// puts(tmp_aligned);
				// aligned[i].resize(tmp_size + 1);
				aligned[i][j] += tmp_aligned;

				fgets(tmp_aligned ,SEQ_MAX_LENGTH,_aligned_fp);
			} 

		}
			
		fclose(_aligned_fp);	
	}
	free(tmp_aligned);


	if (outputfile == NULL)
	{
		outputfile = (char *)malloc(sizeof(char) * PATH_MAX_LEN);
		strcpy(outputfile, inputfile);
		strcpy(outputfile + strlen(inputfile), ".fmalign");
		puts(outputfile);
	}



// ---------------------------------*----------------------------------

	FILE* _out_fp = fopen(outputfile, "wb+");
	if (_out_fp == NULL)
	{
		fprintf(stdout, "outputfile Failed to make %s! will exit ...\n", outputfile);	
	}


	for (uint32_t i = 0; i < njob; ++i)//每次写入一行
	{
		write_raw(seqset[i].name, seqset[i].name_size, _out_fp);
		write_raw("\n", 1, _out_fp);
		for (uint32_t k = 0; k < chain_size_; k++)
		{
			write_raw(aligned[k][i].c_str(), aligned[k][i].size(), _out_fp);//写每个文件的第i行
			write_raw(chain_seq[k], strlen(chain_seq[k]), _out_fp);
		}
		write_raw(aligned[chain_size_][i].c_str(), aligned[chain_size_][i].size(), _out_fp);
		
		write_raw("\n", 1, _out_fp);
	}

	fclose(_out_fp);



	std::string cmd = "rm -f ";
	cmd.append("subfile*.*");
	system(cmd.c_str());

}

void align_direct()                                //If the number of commonseeds is 0 without writting sub-files
{
	std::string cmnd;
	std::string file_name;			
	file_name = set_filename(0);
	if (outputfile == NULL)
	{
		outputfile = (char *)malloc(sizeof(char) * PATH_MAX_LEN);
		strcpy(outputfile, inputfile);
		strcpy(outputfile + strlen(inputfile), ".fmalign");
		puts(outputfile);
	}

	FILE* _fp = fopen(file_name.c_str(), "wb+");//用于保存子文件和比对结果文件的文件名
	if (_fp == NULL)
	{
		fprintf(stdout, "Failed to make %s! will exit ...\n", file_name.c_str());	
	}


	write_raw(inputfile, strlen(inputfile), _fp);
	write_raw("\n", 1, _fp);
	write_raw(outputfile, strlen(outputfile), _fp);
	write_raw("\n", 1, _fp);
	
	fclose(_fp);


	switch (pkg)
	{
		case HALIGN:
			cmnd.append("java -cp ")
				.append("packages/HAlign:")
				.append("packages/HAlign/HAlign2.1.jar halignWrapper \"");
			break;
		case MAFFT:
			cmnd.append("packages/MAFFT/mafftWrapper \"");
			break;
	}

	
	cmnd.append(file_name).append("\" ")
		.append(std::to_string(1).append(" "))
		.append(std::to_string(0))
		.append(" ").append(std::to_string(0))
		.append(" ").append(std::to_string(threads));

	cmnd.append(" > /dev/null");
	system(cmnd.c_str());

	std::string cmd = "rm -f ";
	cmd.append("subfile*.*");
	system(cmd.c_str());
}





