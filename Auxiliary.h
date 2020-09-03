#ifndef _AUXILIARY_H
#define _AUXILIARY_H


#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>



#define SEQ_MAX_LENGTH	5000000
#define SEQNAME_MAX_SIZE 2000
#define KMER_MAX_SIZE 39
#define N 5000000 
#define PATH_MAX_LEN 1000
#define PER_READ_LEN 128

enum PACKAGE {  MAFFT = 0, HALIGN = 1 };



typedef unsigned     short int uint16_t;
typedef unsigned           int uint32_t;




extern int threads;
extern char *inputfile;
extern char *outputfile;
extern char *inputfile;
extern char *outputfile;
extern int njob;
extern char *joinseqs;
extern uint32_t *seqsrank;//这个就是记录上面拼接串中的每个串的结束位置，相当于索引的一部分
// extern int nlenmax;
// extern int nlenmim;
extern unsigned int comseedsize;
extern uint32_t chain_size_;
extern PACKAGE pkg;
extern uint32_t threshold;






#endif