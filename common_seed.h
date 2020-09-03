#ifndef _COMMON_SEED_H
#define _COMMON_SEED_H

#include <vector>
#include "Auxiliary.h"

using std::vector;



struct _seqseed//一个种子在每条序列上对应的位置和频次
{
	uint32_t freqonSEQ;
	uint32_t *locLIST;

};

struct _seed//一个公共种子的id及在所有序列上的信息
{
	// char *seedSTR;
	_seqseed *loconSEQ;

};



void setkvalue();
void setkvalue_sqrt();
_seed * search_from_bwt_more_than_3(char *read);
void search_from_bwt_less_than_4(char *read);
_seed * searchBWT(char *read);
int cmpfunc (const void * a, const void * b);
uint32_t fileter_threshold(uint32_t *freq, uint32_t seednum);
bool getcommonseed(unsigned int hits, unsigned int *locates, _seqseed *comseed);

extern int K;

inline uint32_t getseqid(unsigned int location)
{

	int low, high, mid;
	unsigned int tmpsp , tmpep ,tmploc;
	tmploc = location + 1;
	int flag = 0;
	low = 1;
	high = njob;


	while(low <= high)
	{
		mid = (low + high) / 2;
		// printf("low---high---mid :: %d---%d---%d\n", low, high, mid);
		tmpep = seqsrank[mid];
		tmpsp = seqsrank[mid - 1];
		if( tmploc <= tmpep - K + 1 && tmploc > tmpsp  )
		{
			flag = 1;
			break;
		}
		if (tmploc <= tmpsp) 
		{
			high = mid - 1;
		}
		if (tmploc > tmpep)
		{
			low = mid + 1;
		}
		if (tmploc > tmpep - K + 1 && tmploc <= tmpep) 
		{
			// fprintf(stderr, "location %u is join interval !\n", location);
			break;
		}

	}

	if (flag == 0)
	{
		// fprintf(stderr, "location %u is not exist !\n", location);
		return (0);
	}
	else
	{
		return mid;
	}
}





//下面两个方法用于_seqseed的初始化和重置
// inline void setss(_seqseed ss)
// {
// 	ss.freqonSEQ = 0;
// 	ss.locLIST = (unsigned int *)malloc(SEQ_MAX_LENGTH * sizeof(unsigned int));
// }
// inline void pushss()
// inline void refreshss(_seqseed ss)
// {
// 	ss.freqonSEQ = 0;
// 	// ss.locLIST.clear();
// 	free(ss.locLIST);
// }


//_seed的添加数据
// inline void set_seed(_seed sl)
// {
// 	sl.seqsID = 0;
// }
// inline void push_seed(_seed sl, int seqid, _seqseed ss)
// {
// 	sl.seqsID.push_back(seqid);
// 	sl.loconSEQ.push_back(&ss);
// }
// inline void clear_seed(_seed sl)
// {
// 	sl.seqsID.clear();
// 	sl.loconSEQ.clear();
// }
// inline void print_seed(_seed sl)
// {
// 	int count = sl.seqsID.size();
// 	for (int i = 0; i < count; ++i)
// 	{
// 		fprintf(stdout, "%d: ", sl.seqsID[i]);
// 		int c = sl.loconSEQ[i]->freqonSEQ;
// 		for (int j = 0; j < c; ++j)
// 		{
// 			fprintf(stdout, "%u, ", sl.loconSEQ[i]->locLIST[j]);
// 		}
// 		fprintf(stdout, "\n");
// 	}
// }

// int k;
// int average = seqsrank[njob-1]/njob;
// inline int getseqid(unsigned int location)
// {
// 	// int start = location / average;
// 	int low, high, mid;
// 	unsigned int tmpsp , tmpep;
// 	int flag = 0;
// 	low = 0;
// 	high = njob -1;


// 	while(low < high)
// 	{
// 		mid = (low + high) / 2;
// 		tmpep = seqsrank[mid];
// 		tmpsp = seqsrank[mid - 1];
// 		if( location <= tmpep - k + 1 && location > tmpsp  )
// 		{
// 			flag = 1;
// 			break;
// 		}
// 		else if (location <= tmpsp - k + 1) 
// 		{
// 			high = mid - 1;
// 		}
// 		else if (location > tmpep)
// 		{
// 			low = mid + 1;
// 		}
// 		else
// 		{
// 			fprintf(stderr, "this location is join interval !\n");
// 			break;
// 		}

// 	}

// 	if (flag == 1)
// 	{
// 		return -1;
// 	}
// 	else
// 	{
// 		return mid;
// 	}

// }

#endif 
