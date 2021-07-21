#include <stdio.h>
#include <string.h>
#include <cstring>
#include <stdlib.h>
#include <cstdlib>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <algorithm>
#include <limits.h>
#include "load_seqs.h"



int njob;
char *joinseqs;
uint32_t *seqsrank;

uint32_t maxseqID;
uint32_t minseqID;
Seq *seqset;

uint32_t get_file_length(FILE *fp)
{
	// fp = fopen (filename, "r");

	// if (fp == NULL)
	// {
	// 	fprintf(stdout, "Failed to open %s! will exit ...\n", filename);
	// 	return 0;
	// }

	uint32_t file_length = 0;//文件的总长度
	fseek(fp, 0, SEEK_END);
	file_length = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	fprintf(stderr, "file length is %u\n", file_length);
	return file_length ;
}


char *AllocateCharVec( int l1 )
{
	char *cvec;
	
	cvec = (char *)calloc( l1, sizeof( char ) );
	if( cvec == NULL )
	{
		fprintf( stderr, "Cannot allocate %d character vector.\n", l1 );
		exit( 1 );
	}
	// else{
	// 	fprintf(stderr, "allocate %d character vector\n", l1);
	// }
	return( cvec );
}


int countKUorWA( FILE *fp )
{
	int value;
	char c, b;

	value= 0;
	b = '\n';
	//fprintf(stderr, "%d\n", value);
	while(!feof( fp ))
	{
		c = getc( fp );
		//fprintf(stderr, "%c\n", c);
		if( b == '\n' && ( c == '>' ) )
			value++;
			//fprintf(stderr, "%d\n", value);
		b = c;
	}
	rewind( fp );

	return( value );
}

void searchKUorWA( FILE *fp )
{
	int c, b;
	b = '\n';
	while( !( ( ( c = getc( fp ) ) == '>' || c == EOF ) && b == '\n') )
		b = c;
	ungetc( c, fp );
}


int myfgets(char *s, int l, FILE *fp)  /* l°Ê¾å¤Ï¡¢¹ÔËö¤Þ¤ÇÆÉ¤ßÈô¤Ð¤¹ */
{
        int     c = 0, i = 0 ;

		if( feof( fp ) ) return( 1 );

		for( i=0; i<l && ( c=getc( fp ) ) != '\n'; i++ ) 
        	*s++ = c;
        *s = '\0' ;
		if( c != '\n' ) 
			while( getc(fp) != '\n' )
				;
		return( 0 );
}

int charfilter( char *str )
{
	char tmp;
	char *res = str;
	char *bk = str;

	while( (tmp=*str++) )
	{
//		if( tmp == '=' || tmp == '*' || tmp == '<' || tmp == '>' || tmp == '(' || tmp == ')' )
		if( tmp == '=' || tmp == '<' || tmp == '>' || tmp == '\n')
		{
			fprintf( stderr, "\n" );
			fprintf( stderr, "Characters '= < >' can be used only in the title lines in the --anysymbol or --text mode.\n" );
			fprintf( stderr, "\n" );
			exit( 1 );
		}
//		if( 0x20 < tmp && tmp < 0x7f )
//		if( 0x0 <=tmp && tmp < 0x100 && 
		if( tmp != 0x0a && tmp != 0x20 && tmp != 0x0d )
//		if( tmp != '\n' && tmp != ' ' && tmp != '\t' ) // unprintable characters mo ok.
		{
			*res++ = tmp;
//			reporterr( "tmp=%d (%c)\n", tmp, tmp );
		}
	}
	*res = 0;
	return( res - bk );
}

char *load1SeqWithoutName_realloc( FILE *fpp )
{
	int c, b;
	char *cbuf;
	int size = N;
	char *val;

	val = (char *)malloc( (size+1) * sizeof( char ) );
	cbuf = val;

	b = '\n';
	while( ( c = getc( fpp ) ) != EOF &&        
          !( ( c == '>' || c == 'EOF') && b == '\n' ) )
	{

		if (!isspace((char)c) && (char)c != '\n')*cbuf++ = (char)c;  
		if( cbuf - val == size )
		{
			size += N;
			fprintf( stderr, "reallocating...\n" );
			val = (char *)realloc( val, (size+1) * sizeof( char ) );
			if( !val )
			{
				fprintf( stderr, "Allocation error in load1SeqWithoutName_realloc \n" );
				exit( 1 );
			}
			fprintf( stderr, "done.\n" );
			cbuf = val + size - N;
		}
		b = c;
	}
	ungetc( c, fpp );
	*cbuf = 0;

	return( val );
}

int countnogaplen( char *seq )
{
	int val = 0;
	while( *seq )
		if( *seq++ != '-' ) val++;
	return( val );
}


void _replace(char *seq, unsigned int length)
{
	char *p = seq;
	int i = 0;

	while(*p)
	{
		if(*p != '\n')
			seq[i++]=*p;
		p++;
	}
}

int countATGCandN( char *s, int *countN, int *total )
{
	int nATGC;
	int nChar;
	int nN;
	char c;
	nN = nATGC = nChar = 0;

	if( *s == 0 ) 
	{
		*total = 0;
		return( 0 );
	}

	do
	{
		c = tolower( *s );
		if( isalpha( c ) )
		{
			nChar++;
			if( c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'u' || c == 'n' )
				nATGC++;
			if( c == 'n' )
				nN++;
		}
	}
	while( *++s );

//	reporterr( "nN = %d", nN );

	*total = nChar;
	*countN = nN;
	return( nATGC );
}

void replace_N(char* seqs, unsigned int length)
{
	long long i = 0;
	char buffer;

	srand((unsigned int)time(0));
	char randome_char[4];
	randome_char[0] = 'A';
	randome_char[1] = 'C';
	randome_char[2] = 'G';
	randome_char[3] = 'T';

	for (i = 0; i < length; ++i)
	{
		if (seqs[i] != 'A'&&
			seqs[i] != 'C'&&
			seqs[i] != 'G'&&
			seqs[i] != 'T'&&
			seqs[i] != 'a'&&
			seqs[i] != 'c'&&
			seqs[i] != 'g'&&
			seqs[i] != 't')
		{
			buffer = randome_char[rand() % 4];
			seqs[i] = buffer;
		}
	}
}

void getnumlen_nogap_countn( FILE *fp, int *nlenminpt, double *nfreq )
{
	int total;
	int nsite = 0;
	int atgcnum, nnum, nN;
	int i, tmp, njob, nlenmax;
	char *tmpseq, *tmpname;
	double atgcfreq;
	tmpname = AllocateCharVec( SEQNAME_MAX_SIZE );
	njob = countKUorWA( fp );//得到序列总数
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	nnum = 0;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, SEQNAME_MAX_SIZE-1, fp );//得到一个序列
		seqset[i].name = tmpname;
		tmpseq = load1SeqWithoutName_realloc( fp );
		// allocatewithoutsapce(tmpseq);
		tmp = countnogaplen( tmpseq );
		fprintf(stderr, "%s\n%s\n", tmpname, tmpseq);
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		atgcnum += countATGCandN( tmpseq, &nN, &nsite );//nN是序列中N的个数，nsite是序列的总长度
		total += nsite;
		nnum += nN;
		free( tmpseq );
	}
	free( tmpname );
	atgcfreq = (double)atgcnum / total;
	*nfreq = (double)nnum / atgcnum;

}

void getSeqs_obo(FILE *fp)
{
	uint32_t tmp;
	char *tmpseq, *tmpname;
	
	njob = countKUorWA( fp );

	searchKUorWA( fp );

	unsigned int nlenmax = 0;
	unsigned int nlenmin = INT_MAX;
	
	seqset = (Seq *)malloc(sizeof(Seq) * njob);

	for (int i = 0; i < njob; ++i)
	{
		tmpname = AllocateCharVec( SEQNAME_MAX_SIZE );
		myfgets( tmpname, SEQNAME_MAX_SIZE-1, fp );//得到一个序列
		
		tmpseq = load1SeqWithoutName_realloc( fp );
		tmp = charfilter(tmpseq);


		// nlenmim = std::min(nlenmim, tmp);
		// nlenmax = std::max(nlenmax, tmp);
		if (nlenmin > tmp) 
			minseqID = i;
			nlenmin = tmp;
		if (nlenmax < tmp) 
			maxseqID = i;
			nlenmax = tmp;

		seqset[i].name = tmpname;
		seqset[i].name_size = strlen(tmpname);
		seqset[i].seq = tmpseq;
		seqset[i].seq_size = tmp;

	}
	//处理最后一个序列
	// myfgets( tmpname, SEQ_MAX_LENGTH-1, fp );//得到一个序列

	// tmpseq = load1SeqWithoutName_realloc( fp );
	// tmp = charfilter(tmpseq);

	// nlenmim = std::min(nlenmim, tmp);
	// nlenmax = std::max(nlenmax, tmp);

	// seqset[i].name = tmpname;
	// seqset[i].seq = tmpseq;
	// seqset[i].seqsize = tmp;

}


void getRef(FILE *fp)//得到序列作为建索引的基因
{
	uint32_t tmplen, filelen;
	tmplen = 0;
	filelen = get_file_length( fp );
	joinseqs = AllocateCharVec( filelen );
	seqsrank = (uint32_t *)malloc(( njob + 1) * sizeof(uint32_t));
	seqsrank[0] = 0;

	for (int i = 0; i < njob; ++i)
	{
		memcpy(joinseqs + tmplen, seqset[i].seq, seqset[i].seq_size);
		tmplen += seqset[i].seq_size;
		seqsrank[i + 1] = tmplen;
		fprintf(stdout, "%u\n", seqset[i].seq_size);
	}
	joinseqs[tmplen] = '\0';
}


void freeSeqSet(Seq *s)
{
	for (int i = 0; i < njob; ++i)
	{
		free(s[i].seq);
		free(s[i].name);
	}
	free(s);
}
