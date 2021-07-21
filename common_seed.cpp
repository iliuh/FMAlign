#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <cstring>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <ctype.h>
#include <algorithm>
#include "common_seed.h"
#include "FMtree/bwt.h"
#include "load_seqs.h"
#include "Auxiliary.h"
#include "chain.h"
#include "minimizer.h"

unsigned int comseedsize;
int K;
uint32_t threshold;


bool even(int num)
{
	bool flag = flag;
	if ( ( num & 1 ) == 0) return true;

	return flag;
}


void setkvalue()
{

	int readlen = seqsrank[1];
	int tmp;
	tmp = ceil(log(readlen)/log(2));
	K = std::min(tmp , KMER_MAX_SIZE);
	fprintf(stderr, "readlen is %u kvalue is %d\n", readlen, K);

}

void setkvalue_sqrt()
{

	int readlen = seqset[minseqID].seq_size;
	int tmp;
	tmp = ceil(log2(readlen));
	// if (even(tmp)) tmp += 1;
	// else 
	// {
	// 	tmp += 2;
	// }

	
 	K = std::min(tmp*2+1 , KMER_MAX_SIZE);


	fprintf(stderr, "read_length=%u  kvalue=%d\n", readlen, K);

}



int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

//这里是传进来一个seed的信息
bool getcommonseed(unsigned int hits, unsigned int *locates, _seqseed *comseed)
{
	
	int _com_seqcount = 0;//拥有共同种子的序列
	_seqseed * noise_ss;     //这是最初始得到的一些种子
	int seqsn = -1;        //记录序列数

	
	uint32_t tmpid = 0;
	int tmpfreq;

	std::sort(locates, locates + hits);
	

	
	uint32_t *seqsid;
	seqsid = (uint32_t *)malloc(sizeof(uint32_t) * hits);
	for (int i = 0; i < hits; ++i)
	{
		seqsid[i] = getseqid(locates[i]);
	}

	
	for (int i = 0; i < hits; ++i)
	{
		if ((seqsid[i] > 0) && (seqsid[i] != tmpid))//此时建立一条新的seqseed
		{

			tmpid = seqsid[i];//tmpid都大于1
			// fprintf(stdout, "id %u\n", tmpid);
			tmpfreq = 0;
			++ seqsn;

			comseed[tmpid-1].locLIST = (unsigned int*)malloc(sizeof(unsigned int) * (tmpfreq + 1));

			// comseed[tmpid-1].freqonSEQ = tmpfreq + 1;
			comseed[tmpid-1].locLIST[tmpfreq] = locates[i] - seqsrank[seqsid[i] - 1];
			comseed[tmpid-1].freqonSEQ = ++ tmpfreq;
			
			// tmpfreq++;

		}
		else if ((seqsid[i] > 0) && (seqsid[i] == tmpid))
		{

			comseed[tmpid-1].locLIST = (unsigned int*)realloc(comseed[tmpid-1].locLIST, 
				sizeof(unsigned int) * (tmpfreq + 1));
			// comseed[tmpid-1].freqonSEQ = tmpfreq + 1;
			comseed[tmpid-1].locLIST[tmpfreq] = locates[i] - seqsrank[seqsid[i] - 1];
			comseed[tmpid-1].freqonSEQ = ++tmpfreq;
			
			// tmpfreq++;
		}
		else//这个位置是在序列拼接的间隙上
		{
			continue;
		}
	}

	free(seqsid);

	bool flag = true;
	if (seqsn < njob - 1) flag = false;//这个种子不存在所有序列上，应该过滤





	return flag;


}

uint32_t fileter_threshold(uint32_t *freq, uint32_t seednum)
{

	struct freq_count
	{
		uint32_t tmp_req;//频率
		uint32_t tmp_count;//对应频率出现次数
	};

	freq_count *tmp_freq_count;
	tmp_freq_count = (freq_count*)malloc(sizeof(freq_count) * seednum);
	uint32_t counter = 0;
	std::sort(freq, freq+seednum);

	tmp_freq_count[0].tmp_count = 1;
	tmp_freq_count[0].tmp_req = freq[0];


/* 	unsigned int zero_end = 0;
	while (freq[zero_end] == 0)
	{
		zero_end++;
	}
	tmp_freq_count[0].tmp_count = 1;
	tmp_freq_count[0].tmp_req = freq[zero_end]; */
	

	for(int i = 0; i < seednum-1; ++i)
    {
    	if(freq[i]!= freq[i+1])
		{
			++counter;
			tmp_freq_count[counter].tmp_count = 1;
			tmp_freq_count[counter].tmp_req = freq[i+1];

		}
    	else
			tmp_freq_count[counter].tmp_count ++;
    }

	size_t total = seednum;
	size_t now = 0;
	size_t index;
	if (tmp_freq_count[0].tmp_req == 0)
	{
		total = seednum - tmp_freq_count[0].tmp_count;
		
	}
	


	counter = counter + 1;
	for (index = 0; index < counter; index++)
	{
		if (tmp_freq_count[index].tmp_req > 0)
		{
			// fprintf(stdout, "%u %u\n",tmp_freq_count[index].tmp_req, tmp_freq_count[index].tmp_count);
			now += tmp_freq_count[index].tmp_count;
			if (now > total * 0.7) break;
			
		}		
	}
	
				

	uint32_t threshold;
	threshold = tmp_freq_count[index].tmp_req;
	free(tmp_freq_count);
	return threshold;
}



/*--------------------test minimizer-----------------------*/
// _seed *search_from_bwt_more_than_3(char *read)
// {
// 	uint32_t length_read;
// 	// char* read;
// 	unsigned int top, bot, t;
// 	unsigned int pre_top, pre_bot;


// 	length_read = strlen(read);//需要传入一个read

	
// 	// fprintf(stdout, "read length %u\n", length_read);
// 	uint32_t* locates;


// 	long long number_of_locations = 0;
// 	unsigned int tmp_SA_length = 0;
// 	long long number_of_hits = 0;


// 	setkvalue_sqrt();

	
// 	_seed *common_seed_list;    //存放了所有公共kmer信息，用于下一个阶段的相兼容的path的构造
// 	unsigned int tmp_comseednum;    //记录共同种子的个数
// 	tmp_comseednum = 0;
// 	common_seed_list = (_seed *)malloc(sizeof(_seed) * 1);
	
// 	char pattern[K + 1];
// 	uint32_t *hit_seq;
// 	// hit_seq = (uint32_t *)malloc((length_read - K + 1)* sizeof(uint32_t) );
// 	hit_seq = (uint32_t *)calloc((length_read - K + 1), sizeof(uint32_t) );


// 	unsigned int *mm;
// 	unsigned int mm_num;
// 	mm = (unsigned int*)malloc(sizeof(unsigned int) * (length_read - K + 1));
// 	mm_num = minimize(read, K, mm);//得到了所有的minimizer



	
	
// 	for (uint32_t i = 0; i < mm_num; ++i)
// 	{
// 		if (hit_seq[i] != 0)//说明已经统计过频率
// 		{
// 			continue;
// 		}
		
// 		convert_unint_to_chars(mm[i], pattern, K);
// 		// fprintf(stdout, "%s\n", pattern);

// 		number_of_hits 
// 			= count(pattern, K, &top, &bot, &pre_top, &pre_bot);


// 		if (number_of_hits <njob)
// 		{
// 			continue;
// 		}

// 		locates = (uint32_t *)malloc(number_of_hits * sizeof(uint32_t));
		

// 		tmp_SA_length = 0;
// 		locate(pattern, top, bot,
// 			pre_top, pre_bot, locates, K, &tmp_SA_length);


// 		bool ff = true;
// 		_seqseed *ss;               //存放的是一个公共kmer对应的信息
// 		ss = (_seqseed *)malloc(sizeof(_seqseed) * njob);

// 		//这里改成了不过滤高频
// 		ff = getcommonseed(number_of_hits, locates, ss);  
// 		uint32_t *tmp_loc;
// 		uint32_t tmp_loc_size;

// 		if (ff)//这个种子出现在所有序列上
// 		{
// 			tmp_loc = ss[minseqID].locLIST;
// 			tmp_loc_size = ss[minseqID].freqonSEQ;
// 			for (size_t k = 0; k < tmp_loc_size; k++)//统计频率
// 			{
// 				hit_seq[tmp_loc[k]] = tmp_loc_size;
// 				// fprintf(stdout, "hit_seq %u\n", tmp_loc_size);
// 			}
			
// 			common_seed_list[tmp_comseednum].loconSEQ = ss;
	
// 			tmp_comseednum ++;	
// 			common_seed_list = (_seed *)realloc(common_seed_list, sizeof(_seed) * (tmp_comseednum + 1));
// 		}
// 		else //如果这个不是公共种子，但是肯定出现在取种子的序列上
// 		{
// 			tmp_loc = ss[minseqID].locLIST;
// 			tmp_loc_size = ss[minseqID].freqonSEQ;
// 			for (size_t k = 0; k < tmp_loc_size; k++)//统计频率
// 			{
// 				hit_seq[tmp_loc[k]] = tmp_loc_size;
// 			}

			
// 			free(ss);
// 		}

// 		free(locates);
// 	}

// 	fprintf(stdout, "common seeds %u without filter\n",  tmp_comseednum);
// 	//过滤掉高频的种子，通过这个阈值再将高频过滤一遍
// 	bool flag;
// 	threshold = fileter_threshold(hit_seq, length_read - K + 1);//hit_seq每个种子的频率
// 	fprintf(stdout, "threshold %u\n", threshold);

// 	_seed *common_seed;    //存放了所有公共kmer信息，用于下一个阶段的相兼容的path的构造
// 	unsigned int comseednum;    //记录共同种子的个数
// 	comseednum = 0;
// 	// vector<_seed> common_seed(tmp_comseednum);
// 	common_seed = (_seed *)malloc(sizeof(_seed) * tmp_comseednum);

// 	free(hit_seq);
// 	for (uint32_t i = 0; i < tmp_comseednum; i++)
// 	{
// 		flag = true;
// 		for (int j = 0; j < njob; ++j)
// 		{
// 			if (common_seed_list[i].loconSEQ[j].freqonSEQ > threshold)//如果种子频率太高，就过滤
// 			{
// 				// fprintf(stdout, "%u freq %u\n",i,  common_seed_list[i].loconSEQ[j].freqonSEQ);
// 				flag = false;
// 				break;
// 			}	
// 		}

// 		if (flag)
// 		{
		
// 			memcpy(&common_seed[comseednum], &common_seed_list[i],sizeof(_seed));
			
// 			comseednum ++;
// 		}	
// 	}
	
	

// 	comseedsize = comseednum;
// 	fprintf(stdout, "common_seed=%d \n", comseednum);
// 	free(common_seed_list);
// 	// vector<uint32_t>().swap(hit_seq);
// 	return common_seed;
// }



/* ---------------------test longer seeds--------------------- */
_seed *search_from_bwt_more_than_3(char *read)
{
	uint32_t length_read;
	// char* read;
	unsigned int top, bot, t;
	unsigned int pre_top, pre_bot;


	length_read = strlen(read);//需要传入一个read

	
	// fprintf(stdout, "read length %u\n", length_read);
	uint32_t* locates;


	long long number_of_locations = 0;
	unsigned int tmp_SA_length = 0;
	long long number_of_hits = 0;


	setkvalue_sqrt();
	// K = 39;

	
	_seed *common_seed_list;    //存放了所有公共kmer信息，用于下一个阶段的相兼容的path的构造
	unsigned int tmp_comseednum;    //记录共同种子的个数
	tmp_comseednum = 0;
	common_seed_list = (_seed *)malloc(sizeof(_seed) * 1);
	
	char pattern[K + 1];
	uint32_t *hit_seq;
	// hit_seq = (uint32_t *)malloc((length_read - K + 1)* sizeof(uint32_t) );
	hit_seq = (uint32_t *)calloc(length_read - K + 1, sizeof(uint32_t) );

	for (uint32_t i = 0; i < length_read - K + 1; ++i)
	{
		if (hit_seq[i] != 0)//说明已经统计过频率
		{
			continue;
		}
		
		memcpy(pattern, read + i, K);

		number_of_hits 
			= count(pattern, K, &top, &bot, &pre_top, &pre_bot);


		if (number_of_hits <njob)
		{
			continue;
		}

		locates = (uint32_t *)malloc(number_of_hits * sizeof(uint32_t));
		

		tmp_SA_length = 0;
		locate(pattern, top, bot,
			pre_top, pre_bot, locates, K, &tmp_SA_length);


		bool ff = true;
		_seqseed *ss;               //存放的是一个公共kmer对应的信息
		ss = (_seqseed *)malloc(sizeof(_seqseed) * njob);

		//这里改成了不过滤高频
		ff = getcommonseed(number_of_hits, locates, ss);  
		uint32_t *tmp_loc;
		uint32_t tmp_loc_size;





			

		if (ff)//这个种子出现在所有序列上
		{
			tmp_loc = ss[minseqID].locLIST;
			tmp_loc_size = ss[minseqID].freqonSEQ;
			for (size_t k = 0; k < tmp_loc_size; k++)//统计频率
			{
				hit_seq[tmp_loc[k]] = tmp_loc_size;
			}


			common_seed_list[tmp_comseednum].loconSEQ = ss;
	
			tmp_comseednum ++;	
			common_seed_list = (_seed *)realloc(common_seed_list, sizeof(_seed) * (tmp_comseednum + 1));
		}
		else //如果这个不是公共种子，但是肯定出现在取种子的序列上
		{
/* 			tmp_loc = ss[minseqID].locLIST;
			tmp_loc_size = ss[minseqID].freqonSEQ;
			for (size_t k = 0; k < tmp_loc_size; k++)//统计频率
			{
				hit_seq[tmp_loc[k]] = tmp_loc_size;
				// fprintf(stdout, " %u\n", tmp_loc_size);
			}
 */
			
			free(ss);
		}

		free(locates);
	}


	bool flag;
	_seed *common_seed;    //存放了所有公共kmer信息，用于下一个阶段的相兼容的path的构造
	unsigned int comseednum;    //记录共同种子的个数
	comseednum = 0;




	// fprintf(stdout, "common seeds %u without filter\n",  tmp_comseednum);



		// for (size_t i = 0; i < length_read - K + 1; i++)
		// {
			
		// 	fprintf(stdout, " %u\n", hit_seq[i]);
		// } 

	if (tmp_comseednum == 0)
	{
		common_seed == NULL;
		return common_seed;
	}



	else
	{
		common_seed = (_seed *)malloc(sizeof(_seed) * tmp_comseednum);
		//过滤掉高频的种子，通过这个阈值再将高频过滤一遍
		threshold = fileter_threshold(hit_seq, length_read - K + 1);//hit_seq每个种子的频率


		

		free(hit_seq);
		for (uint32_t i = 0; i < tmp_comseednum; i++)
		{
			flag = true;
			for (int j = 0; j < njob; ++j)
			{
				if (common_seed_list[i].loconSEQ[j].freqonSEQ > threshold)//如果种子频率太高，就过滤
				{
					// fprintf(stdout, "%u freq %u\n",i,  common_seed_list[i].loconSEQ[j].freqonSEQ);
					flag = false;
					break;
				}	
			}

			if (flag)
			{
			
				memcpy(&common_seed[comseednum], &common_seed_list[i],sizeof(_seed));
				// for (size_t j = 0; j < njob; j++)
				// {
				// 	fprintf(stdout, "%u freq %u\n",j,  common_seed[comseednum].loconSEQ[j].freqonSEQ);
				// }
				
				comseednum ++;
			}	
		}
	
	

		comseedsize = comseednum;
		fprintf(stdout, "filtered common seed %d \n", comseednum);
		free(common_seed_list);
		// vector<uint32_t>().swap(hit_seq);
		return common_seed;
	}
	


}


void search_from_bwt_less_than_4(char *read)
{
	long long i, j;
	int length_read;
	// char* read;
	unsigned int top, bot, t;
	unsigned int pre_top, pre_bot;
	char *pattern;


	length_read = strlen(read);

	unsigned int* locates;

	// double start = clock();


	long long number_of_locations = 0;



	unsigned int tmp_SA_length = 0;





	long long number_of_hits = 0;




	// struct  timeval  start_timeval;
	// struct  timeval  end_timeval;
	// unsigned long timer;
	// gettimeofday(&start_timeval, NULL);

	// int k;
	// setkvalue();
	setkvalue_sqrt();
	// K = 3;
	for (int i = 0; i < length_read - K + 1; ++i)
	{
		pattern = (char *)malloc(K*sizeof(char));
		strncpy(pattern, read, K);
		number_of_hits = count_less_than_4(pattern, K, &top, &bot);
		if (number_of_hits < njob)
		{
			//read = read + 1;
			continue;
		}

		locates = (unsigned int *)malloc(number_of_hits*sizeof(unsigned int));

		tmp_SA_length = 0;
		locate_less_than_4(pattern, top, bot, pre_top, pre_bot, locates, K, &tmp_SA_length);
		free(locates);
		free(pattern);
		number_of_locations = number_of_locations + tmp_SA_length;
	}
	free(read);


}


//传出seed的所有信息
_seed *searchBWT(char *read)
{

	_seed *commonseed;
	load_index(inputfile);
	commonseed = search_from_bwt_more_than_3(read);//这暂时不管SA值

	// if (commonseed != NULL) fprintf(stderr, "not null\n");
	return commonseed;
	
	//return commonseed;
	// if (bitmapper_index_params.compress_sa >= 4)
	// {
	// 	search_from_bwt_more_than_3(read);
	// }
	// else
	// {
	// 	search_from_bwt_less_than_4(read);
	// }	

}
