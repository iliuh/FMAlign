#ifndef _LOAD_SEQS_H
#define _LOAD_SEQS_H

#include "Auxiliary.h"





typedef struct
{
	char* name;
	char* seq;
	uint16_t name_size;
	uint32_t seq_size;
} Seq;

extern Seq *seqset;
extern uint32_t maxseqID;
extern uint32_t minseqID;


uint32_t get_file_length(FILE *fp);
char *AllocateCharVec( int l1 );
int countKUorWA( FILE *fp );
void searchKUorWA( FILE *fp );
void getnumlen_nogap_countn( FILE *fp, int *nlenminpt, double *nfreq );
int myfgets(char *s, int l, FILE *fp) ;\
int charfilter( char *str );
char *load1SeqWithoutName_realloc( FILE *fpp );
void allocatewithoutsapce(char *seq);
int countnogaplen( char *seq );
int countATGCandN( char *s, int *countN, int *total );
void replace_N(char* seqs, unsigned int length);

void getSeqs_obo(FILE *fp);
void getRef(FILE *fp);//得到序列作为建索引的基因


void freeSeqSet(Seq *s);


#endif