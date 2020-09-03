#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <stdint.h>
#include <ctype.h>
#include <vector>
#include "FMtree/bwt.h"
#include "Auxiliary.h"
#include "load_seqs.h"
#include "common_seed.h"
#include "subseqs.h"
#include "chain.h"
#include "command_line.h"


int main( int argc, char *argv[] )
{
	FILE *infp;
	int nlenmin;
	double nfreq;

	// inputfile = "data/Variola.fasta";
	CommandLine_parse ( argc, argv);

	if( inputfile )
	{
		infp = fopen( inputfile, "r" );
		if( !infp )
		{
			fprintf( stderr, "Cannot open %s\n", inputfile );
			exit( 1 );
		}
	}
	else
		infp = stdin;



	struct timeval start, fm_time, chain_time, end;
	double timeconsumed;
	gettimeofday(&start, NULL);


	//处理文件得到序列信息
	getSeqs_obo(infp);
	getRef(infp);
	fclose( infp );

	//这里直接调用bwt中的indenpendent_creadte_index，需要传入的参数是序列长，序列，SA值，输入文件名
	//这里的inputfile是拷贝了.index这个串的
	char tmpfile[100];
	strcpy(tmpfile, inputfile);
	puts(tmpfile);
	
	indenpendent_creadte_index(seqsrank[njob], 
	                          &joinseqs, 5, tmpfile);//这个方法以后，joinseqs已经没有字符了

	free( joinseqs );


	char *tmpread;	
	uint32_t readlen = seqset[minseqID].seq_size;
	tmpread = (char *)malloc((readlen + 1) * sizeof(char));
	memcpy(tmpread, seqset[minseqID].seq, readlen);
	tmpread[readlen] = '\0';


	_seed *commonseed;
	commonseed = searchBWT(tmpread);

	free(tmpread);
	free( seqsrank );



	gettimeofday(&fm_time, NULL);
	timeconsumed = fm_time.tv_sec-start.tv_sec +(fm_time.tv_usec-start.tv_usec)/1000000.0;
	fprintf(stdout, "FM-index and search common seed time: %5.3f seconds\n",  timeconsumed);


	std::vector<_chain> chain;
	chain = filter_noise_chain(commonseed, comseedsize);


	gettimeofday(&chain_time, NULL);
	timeconsumed = chain_time.tv_sec-fm_time.tv_sec +(chain_time.tv_usec-fm_time.tv_usec)/1000000.0;
	fprintf(stdout, "Get best chains time: %5.3f seconds\n",  timeconsumed);
	
	
	write_sub_file(chain);




	// finish = clock();
	// duration = (double)(finish - start) / CLOCKS_PER_SEC;
	// fprintf(stdout, "Time: %5.3f seconds\n",  duration);
	// double cost_t = time(NULL) - t;
	// fprintf(stdout, "Time: %5.3f seconds\n",  cost_t);
	
	gettimeofday(&end, NULL);
	timeconsumed = end.tv_sec-start.tv_sec +(end.tv_usec-start.tv_usec)/1000000.0;
	fprintf(stdout, "Time: %5.3f seconds\n",  timeconsumed);


	free(seqset);
	return( 0 );
}





