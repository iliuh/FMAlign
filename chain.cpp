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
inline bool no_cross(_chain chain_first, _chain chain_second, unsigned int mean)
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
			if (no_cross(best_chain.back(), right,77)) 
			{
				best_chain.push_back(right);
			}
			if ( !no_cross(best_chain.back(), right,77) && best_chain.back().wide < right.wide)
			{
				best_chain.pop_back();
				goto loop;
			}
		}
	}
	
	vector<_chain>().swap(all_chain);


	chain_size_ = best_chain.size();
	fprintf(stdout, "best_chain=%u\n", chain_size_);



	return best_chain;
	
}





inline float variance(_chain chain_first, _chain chain_second)
{
	float mean = 0.0;
	float sum = 0.0;
	float var = 0.0;
	int i;
	float* offset;
	offset = (float*)malloc(sizeof(float) * njob);

	for (i = 0; i < njob; ++i)
	{
		offset[i] = (float)chain_second.pos[i] - (float)(chain_first.pos[i] + chain_first.wide);
		sum += offset[i];
	}

	mean = sum / njob;

	i = 0;
	for ( i = 0; i < njob; i++)
	{
		var += (mean-offset[i])*(mean-offset[i]);
	}
	// fprintf(stdout, "   %f\n", var);
	free(offset);
	var = var/njob;
	
	return sqrt(var);
}


std::vector<_chain> creat_optimal_chain(_seed *commonseedlist, uint32_t seedsize)
{
	//种子目前是按出现在read上的顺序来排列的
	uint32_t tmp_chainsize = 0;
	uint32_t tmpsize;

	vector<_chain> noise_chain;
	vector<_chain> tmp_noise_chain;
	vector<_chain> all_chain;


	_chain left;
	_chain right;
	_chain merge_chain;
	bool flag;
	unsigned int i;
	unsigned int j;
	int min;
	unsigned int all_chain_num;
	unsigned int longest_wide;





	//计算所有链
	for (i = 0; i < seedsize; ++i)
	{
		tmp_noise_chain = get_from_seed(commonseedlist[i]);
		noise_chain.insert(noise_chain.end(), tmp_noise_chain.begin(), tmp_noise_chain.end());
		vector<_chain>().swap(tmp_noise_chain);
	}


	//按之前的read的位置排序
	sort(noise_chain.begin(), noise_chain.end(), compare_pos);


	//先拼接链all_chain存放所有的拼接链
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
	




	all_chain_num = all_chain.size();
	longest_wide = all_chain[0].wide;
	unsigned int long_flag = 0;
	for ( i = 1; i < all_chain_num; i++)
	{
		if (all_chain[i].wide > longest_wide )
		{
			longest_wide = all_chain[i].wide;
			long_flag = i;
		}
		
	}
	// mean_wide = mean_wide / all_chain_num;
	fprintf(stdout, "the total number of merged chains is %u \n", all_chain_num);
	fprintf(stdout, "chain %u have long wide %u\n",i,  longest_wide);



	//DP
	uint8_t iter_step_num = 5;
	uint8_t step = 0;

	
	min = INT32_MAX;
	int tmp_min;
	float ecah_anchor_score[all_chain_num];
	int p[all_chain_num];//用于回溯

	float max_pre_socre = __FLT_MIN__;//前驱节点的得分

	unsigned int max_pre_socre_id = 0;
	float tmp_score = 0.0;
	float score;
	float gap_cost = __FLT_MAX__;
	unsigned int k = 0;
	unsigned int start_pos;//开始计算得分的位置





	
	for ( i = 0; i < all_chain_num; i++)
	{
		 ecah_anchor_score[i] =  all_chain[i].wide;
		 p[i] = 0;
	}
	

	unsigned int tmp_step;
	for ( i = 1; i < all_chain_num; i++)
	{

		//首先找到不相交的位置，即开始计算得分的位置
		start_pos = i;
		while (start_pos > 0 && !no_cross( all_chain[start_pos-1], all_chain[i], longest_wide))
		{
			start_pos --;
			// fprintf(stdout, "start_pos %u\n", start_pos);
			//continue;
		}


		tmp_step = (iter_step_num > start_pos ? start_pos : iter_step_num);
		for (step = 0; step < tmp_step; step++)
		{
			//计算罚分
			//检查时候相交，相交罚分为负无穷，否则就是距离的方差
			/*
			if (!no_cross( all_chain[start_pos-step-1], all_chain[i], longest_wide))
			{
				gap_cost = __FLT_MAX__;
				// fprintf(stdout, "%u:%u continue!\n",i, step);
				continue;
			}
			else
			{
				gap_cost = variance(all_chain[start_pos-step-1], all_chain[i]);
				// fprintf(stdout, "%u:%f gap cost\n",i, gap_cost);
			}
			*/



			gap_cost = variance(all_chain[start_pos-step-1], all_chain[i]);

			//计算得分
			for (j = 0; j < njob; j++)
			{
				tmp_min = all_chain[i].pos[j] - all_chain[start_pos-step-1].pos[j] - all_chain[start_pos-step-1].wide;
				if (tmp_min < min)
				{
					min = tmp_min;
					
				}	
			}
			
			
			/* 
			这个最小匹配的得分就是两条链之间的最小距离；
			但是不能让这个距离太大，也不能太小，因此组要再处理一下

			 */
			tmp_score = min;

			
			
			tmp_score =  all_chain[i].wide + ecah_anchor_score[start_pos-step-1] - gap_cost;
			min = INT32_MAX;

			
			if (tmp_score > ecah_anchor_score[i])
			{
				p[i] = start_pos-step-1 + 1;//需要跳转到0号的肯定不能再存0,
				ecah_anchor_score[i] = tmp_score;
				// fprintf(stdout, "%u:%u pi\n",i, p[i]);
			}
			tmp_score = 0.0;
			
		}

	}





	vector<_chain> best_chain;
	vector<vector<int>> best_chain_index;
	vector<int> tmp_chain_index;
	int p_index;
	int tmp_p_index;

	int longest_chain_index = 0; 
	int chain_num = 0;
	float higest_score = 0.0;


	for (i = all_chain_num - 1; i > 0; i--)
	{
		// if (p[i] == -1 && p[i] == 0) continue;
		if (p[i] <= 0) continue;

		
		

		//push完一整个链
		p_index = i;
		// p[p_index] = -1;
		tmp_chain_index.push_back(p_index);
		
		if (ecah_anchor_score[i] > higest_score)
		{
			higest_score = ecah_anchor_score[i];
			longest_chain_index = chain_num;
		}
		
		
		while (p[p_index] > 0)
		{
			tmp_p_index = p[p_index];
			// p[p_index] = -1;
			tmp_chain_index.push_back(tmp_p_index - 1);
			
			p_index = tmp_p_index - 1;	
			
		}
		
		


		best_chain_index.push_back(tmp_chain_index);
		chain_num ++;
		tmp_chain_index.clear();		
		
	}

	
	vector<int>().swap(tmp_chain_index);
	/*
	long_flag = 0;
	for ( i = 0; i < best_chain_index.size(); i++)
	{
		// if (best_chain_index[i].size() > long_flag)//最长的
		if (ecah_anchor_score[best_chain_index[i][0]] > long_flag)//得分最大的
		{
			// long_flag = best_chain_index[i].size();
			long_flag = ecah_anchor_score[best_chain_index[i][0]];
			longest_chain_index = i;
		}
	}
	*/




	unsigned int lll = 0;

	if(all_chain_num == 1) best_chain.push_back(all_chain[0]);
	else{
		for ( i = best_chain_index[longest_chain_index].size() - 1; i >= 0; i--)
		{
			
			best_chain.push_back(all_chain[best_chain_index[longest_chain_index][i]]);
			lll ++;
		}
		vector<vector<int>>().swap(best_chain_index);
	}

	vector<_chain>().swap(all_chain);
	
	
	// fprintf(stdout, "best chain %u\n", best_chain.size());
	// for (size_t i = 0; i < best_chain.size(); i++)
	// {
	// 	for (size_t j = 0; j < njob; j++)
	// 	{
	// 		fprintf(stdout, "%u ", best_chain[i].pos[j]);
	// 	}
	// 	fprintf(stdout, "\n");
	// }
	
	chain_size_ = best_chain.size();
	fprintf(stdout, "best_chain=%u\n", chain_size_);
	return best_chain;

}


void creat_optimal_chain_minimizer(_seed *commonseedlist, uint32_t seedsize)
{
	//种子目前是按出现在read上的顺序来排列的
	uint32_t tmp_chainsize = 0;
	uint32_t tmpsize;

	vector<_chain> noise_chain;
	vector<_chain> tmp_noise_chain;
	vector<_chain> all_chain;


	_chain left;
	_chain right;
	_chain merge_chain;
	bool flag;
	unsigned int i;
	unsigned int j;
	int min;
	unsigned int chain_num;
	unsigned int mean_wide;





	//计算所有链
	for (i = 0; i < seedsize; ++i)
	{
		tmp_noise_chain = get_from_seed(commonseedlist[i]);
		noise_chain.insert(noise_chain.end(), tmp_noise_chain.begin(), tmp_noise_chain.end());
		vector<_chain>().swap(tmp_noise_chain);
	}


	//按之前的read的位置排序
	sort(noise_chain.begin(), noise_chain.end(), compare_pos);



	chain_num = noise_chain.size();
	fprintf(stdout, "chain number %u\n",chain_num);
	mean_wide = noise_chain[0].wide;
	unsigned int long_flag = 0;
	for ( i = 1; i < chain_num; i++)
	{

		if (noise_chain[i].wide > mean_wide)//以前想找个最大的干嘛？
		{
			mean_wide = noise_chain[i].wide;
			long_flag = i;
		}
		
	}
	fprintf(stdout, "chain %u have long wide %u\n",i,  mean_wide);
	



	//DP
	uint8_t iter_step_num = 5;
	uint8_t step = 0;

	
	min = INT32_MAX;
	int tmp_min;
	float ecah_chain_score[chain_num+1];
	unsigned int p[chain_num+1];//用于回溯

	float max_pre_socre = __FLT_MIN__;//前驱节点的得分

	unsigned int max_pre_socre_id = 0;
	float tmp_score = 0.0;
	float score;
	float gap_cost = __FLT_MAX__;
	unsigned int k = 0;
	unsigned int start_pos = 0;//开始计算得分的位置







	// p[chain_num] = {0};
	for ( i = 0; i < chain_num; i++)//每个链的初始得分
	{
		ecah_chain_score[i] =  0.01 * noise_chain[i].wide;
		p[i] = 0;
	}
	

	unsigned int tmp_step;
	for ( i = 1; i < chain_num; i++)
	{

		//首先找到不相交的位置，即开始计算得分的位置
		start_pos = i - 1;
		if (!no_cross( noise_chain[start_pos], noise_chain[i], mean_wide))
		{
			start_pos --;
			continue;
		}

		tmp_step = (iter_step_num > i ? i : iter_step_num);
		for (step = 0; step < tmp_step; step++)
		{
			//计算罚分
			//检查时候相交，相交罚分为负无穷，否则就是距离的方差
			if (!no_cross( noise_chain[i-step-1], noise_chain[i], mean_wide))
			{
				gap_cost = __FLT_MAX__;
				// fprintf(stdout, "%u:%u continue!\n",i, step);
				continue;
			}
			else
			{
				gap_cost = variance(noise_chain[i-step-1], noise_chain[i]);
				// fprintf(stdout, "%u:%f gap cost\n",i, gap_cost);
			}




			//计算得分
			for (j = 0; j < njob; j++)
			{
				tmp_min = noise_chain[i].pos[j] - noise_chain[i-step-1].pos[j] - noise_chain[i-step-1].wide;
				if (tmp_min < min)
				{
					min = tmp_min;
					
				}	
			}
			// fprintf(stdout, "min %u:%u min\n",i, step);
			
			
			/* 
			这个最小匹配的得分就是两条链之间的最小距离；
			但是不能让这个距离太大，也不能太小，因此组要再处理一下

			 */
			// tmp_score = 100.0*exp(-(abs((int)(min - 1000)) / 3));
			// fprintf(stdout, "%u:%.3f distance score\n",i, tmp_score);

			
			tmp_score =  0.01 * min + 0.01 * noise_chain[i].wide + ecah_chain_score[i-step-1] - gap_cost;
			// fprintf(stdout, "%.3f:%f:%f score\n",0.1 * min, ecah_chain_score[i-step-1], gap_cost);
			min = INT32_MAX;

			
			if (tmp_score > ecah_chain_score[i])
			{
				p[i] = i - step-1;
				ecah_chain_score[i] = tmp_score;
				// fprintf(stdout, "%u:%u pi\n",i, p[i]);
			}
			tmp_score = 0.0;
			
		}

	}





	unsigned used_flag_num;
	used_flag_num = chain_num / 32 + 1;
	unsigned int used_vector[used_flag_num] = {0};
	unsigned int shift_length = 0;




	vector<_chain> best_chain;

	

}
