/*
 * @Author: your name
 * @Date: 2021-03-18 20:54:26
 * @LastEditTime: 2021-04-06 20:30:33
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /FMAlign-0.2.6/chain.h
 */
#ifndef _CHAIN_H
#define _CHAIN_H

#include <vector>
#include "common_seed.h"
#include "Auxiliary.h"

using std::vector;
struct _chain//一个链就是一个公共种子在每条序列上的位置组成
{
	// unsigned int chainID;
	uint16_t wide;
	// unsigned int *pos;
	std::vector<uint32_t> pos;
};



double score(unsigned int *chain);
inline bool offset(uint32_t *list1, uint32_t *list2);
inline bool no_cross(uint32_t *list1, uint32_t *list2, unsigned int mean);
inline bool my_equal(uint32_t *list1, uint32_t *list2);
inline float variance(_chain chain_first, _chain chain_second);
std::vector<_chain> merge(std::vector<_chain> left, std::vector<_chain> right);
std::vector<_chain> get_from_seed(_seed seed);
std::vector<_chain> getcompatiblechain(_seed *commonseedlist, unsigned int seedsize);
bool compare_wide(const _chain &chain1, const _chain &chain2) ;
bool compare_pos(const _chain &chain1, const _chain &chain2) ;
std::vector<_chain> filter_noise_chain(_seed *commonseedlist, uint32_t seedsize);
void best_chain(std::vector<_chain> chain);
void getonechain(_seed tmpseed, _chain chain);
std::vector<_chain> getchain(_seed * commonseedlist, unsigned int seedsize);
std::vector<_chain> creat_optimal_chain(_seed *commonseedlist, uint32_t seedsize);
void creat_optimal_chain_minimizer(_seed *commonseedlist, uint32_t seedsize);



void circulate(_seed tmpseed, vector<vector<uint32_t>> &res);
void productImplement(_seed tmpseed,vector<vector<unsigned int>> &res,int index,vector<unsigned int> tmp);
unsigned int **recursion(_seed tmpseed, int index, unsigned int **list, unsigned int listsize);

#endif