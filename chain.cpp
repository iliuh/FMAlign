#include<stdio.h>
#include<stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <time.h>
#include<stdint.h>
#include<ctype.h>
#include "load_seqs.h"
#include "Auxiliary.h"
#include "chain.h"

using std::vector;



bool compare_wide(const _chain &chain1, const _chain &chain2) 
{	

	 return chain1.wide > chain2.wide;
}

bool compare_pos(const _chain &chain1, const _chain &chain2) 
{	

	 return chain1.pos[minseqID] < chain2.pos[minseqID];
}


void circulate(_seed tmpseed, vector<vector<uint32_t>> &res)
{
	uint32_t total = 1;
	uint32_t index = 1;
	vector<uint32_t> s;
	vector<uint32_t> tmp;

	for (uint32_t i = 0; i < njob; ++i)
	{
		tmp.push_back(tmpseed.loconSEQ[i].locLIST[0]);
	}
	res.push_back(tmp);

/*	for (int i = 0; i < njob; ++i){
		total *= tmpseed.loconSEQ[i].freqonSEQ;
	}


	for (uint32_t j = 0; j < tmpseed.loconSEQ[0].freqonSEQ; ++j){
		tmp.push_back(tmpseed.loconSEQ[0].locLIST[j]);
	}


	while(index < njob - 1)
	{
		for (uint32_t i = 0; i < tmpseed.loconSEQ[index].freqonSEQ; ++i)
		{
			s.clear();
			for (unsigned int j = 0; j < tmp.size(); j++){
				s.push_back(tmp[j]);
			}
			s.push_back(tmpseed.loconSEQ[index].locLIST[i]);
		}

		tmp.assign(s.begin(), s.end());
		index ++;
	}
	for (int j = 0; j < tmpseed.loconSEQ[index].freqonSEQ;j++){
        tmp.push_back(tmpseed.loconSEQ[index].locLIST[j]);
		res.push_back(tmp);
		tmp.pop_back();
	}
*/
}


void productImplement(_seed tmpseed,vector<vector<unsigned int>> &res,int index,vector<unsigned int> tmp)
{
	if (index < njob - 1){
		for (int i = 0; i < tmpseed.loconSEQ[index].freqonSEQ; i++){
			vector<unsigned int> s;
			s.clear();
			
			for (unsigned int i = 0; i < tmp.size(); i++){
					s.push_back(tmp[i]);
			}
			s.push_back(tmpseed.loconSEQ[index].locLIST[i]);
			productImplement(tmpseed, res, index+1, s);
		}
	}
	else if (index == njob-1){
		for (int j = 0; j < tmpseed.loconSEQ[index].freqonSEQ;j++){
            tmp.push_back(tmpseed.loconSEQ[index].locLIST[j]);
			res.push_back(tmp);
			tmp.pop_back();
		}
	}
}



// inline bool offset(uint32_t *list1, uint32_t *list2)
inline bool offset(_chain chain_first, _chain chain_second)
{
	int tmp_offset;
	bool flag = true;

	tmp_offset = (int)chain_second.pos[0] - (int)chain_first.pos[0];
	if (tmp_offset 	> chain_first.wide)
	{
		return false;
	}
	
	for (int i = 1; i < njob; ++i)
	{
		if (tmp_offset != ((int)chain_second.pos[i] - (int)chain_first.pos[i]))
		{
			flag = false;
			break;
		}
	}
	return flag;
}

//主要判断前一条链是否跨过了第二条链
inline bool no_cross(_chain chain_first, _chain chain_second)
{
	int tmp_offset;
	bool flag = true;

	for (int i = 0; i < njob; ++i)
	{
		tmp_offset = (int)chain_second.pos[i] - (int)(chain_first.pos[i] + chain_first.wide);
		if (tmp_offset < 1001)
		{
			flag = false;
			break;
		}
	}
	return flag;
}


inline bool my_equal(uint32_t *list1, uint32_t *list2)
{
	bool flag = true;
	for (int i = 0; i < njob; ++i)
	{
		if (list1[i] != list2[i])
		{
			flag = flag;
			break;
		}

	}
	return flag;
}


std::vector<_chain> get_from_seed(_seed seed)
{

	// if (threshold == 1)
	// {
	// 	std::vector<_chain> chain_list(1);
	// 	for (size_t i = 0; i < njob; i++)
	// 	{
	// 		chain_list[0].pos.push_back(seed.loconSEQ[i].locLIST[0]);
	// 	}
	// 	return chain_list;	
	// }
	// else
	// {
		uint16_t size = 1;
		vector<std::vector<uint32_t > > list;
		vector<unsigned int> tmp;

		for (uint32_t i = 0; i < njob; ++i) 
			size  *= seed.loconSEQ[i].freqonSEQ;

		std::vector<_chain> chain_list(size);

		if (threshold == 1)
		{
			circulate(seed, list);
		}
		else
		{
			productImplement(seed, list, 0, tmp);
		}



		for (uint16_t i = 0; i < size; ++i)
		{

			chain_list[i].pos = list[i];
			// fprintf(stdout, "%d\n", list[i].size());
			chain_list[i].wide = K;
		}

		return chain_list;
	// }
}





std::vector<_chain> filter_noise_chain(_seed *commonseedlist, uint32_t seedsize)
{
	//种子目前是按出现在read上的顺序来排列的
	uint32_t tmp_chainsize = 0;
	uint32_t tmpsize;
	// _chain *noise_chain;
	vector<_chain> noise_chain;
	vector<_chain> tmp_noise_chain;
	vector<_chain> all_chain;


	_chain left;
	_chain right;
	_chain merge_chain;
	bool flag;



	for (uint32_t i = 0; i < seedsize; ++i)//计算所有链
	{
		tmp_noise_chain = get_from_seed(commonseedlist[i]);
		noise_chain.insert(noise_chain.end(), tmp_noise_chain.begin(), tmp_noise_chain.end());
		vector<_chain>().swap(tmp_noise_chain);
	}


	//按之前的read的位置排序
	sort(noise_chain.begin(), noise_chain.end(), compare_pos);


	all_chain.push_back(noise_chain[0]) ;
	for (uint32_t n = 1; n < noise_chain.size(); n++)
	{	
		flag = false;
		right = noise_chain[n];
		for (uint32_t i = all_chain.size(); i > 0; i--)
		{
			left = all_chain[i-1];
			// fprintf(stdout, "%u %u %u\n", n, all_chain.size(), i);
			if ( offset(left, right))
			{
				all_chain[i - 1].wide = right.pos[0] + K - left.pos[0];
				flag = true;
				break;
			}
		}
		if(!flag) all_chain.push_back(right);
	}

	vector<_chain>().swap(noise_chain);
	vector<_chain> best_chain;

	for (uint32_t i = 0; i < all_chain.size(); i++)
	{
		right = all_chain[i];
		if (right.wide <= K)
		{
			continue;
		}
		
		loop:
		if(best_chain.empty()) 
		{
			best_chain.push_back(right);
		}
		else//不为空就需要判断是否交叉
		{
			if (no_cross(best_chain.back(), right)) 
			{
				best_chain.push_back(right);
			}
			if ( !no_cross(best_chain.back(), right) && best_chain.back().wide < right.wide)
			{
				best_chain.pop_back();
				goto loop;
			}
		}
	}
	
	vector<_chain>().swap(all_chain);

	// sort(all_chain.begin(), all_chain.end(), compare_wide);
	// for (int i = 0; i < best_chain.size(); ++i)
	// {
	// 	for (int j = 0; j < njob; ++j)
	// 	{
	// 		fprintf(stdout, "%u ", best_chain[i].pos[j]);
	// 	}
	// 	fprintf(stdout, "\n" );
	// 	fprintf(stdout, "%d %u\n", i, best_chain[i].wide);
	// }

	// int tmp_size = (int)floor(log10(seqset[minseqID].seq_size));
	
	// if(all_chain.size() >= pow(2, tmp_size)) 
	// {
	// 	chain_size_ = pow(2, tmp_size);
	// 	best_chain.assign(all_chain.begin(), all_chain.begin() + chain_size_);
	// }
	// else
	// {
	// 	chain_size_ = all_chain.size();
	// 	best_chain.assign(all_chain.begin(), all_chain.end());
	// }


	// sort(best_chain.begin(), best_chain.end(), compare_pos);
	chain_size_ = best_chain.size();
	fprintf(stdout, "best_chain=%u\n", chain_size_);



	return best_chain;
	
}


