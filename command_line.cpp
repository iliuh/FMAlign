#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <string.h>
#include "Auxiliary.h"
#include "command_line.h"


char *inputfile;
char *outputfile;
int threads;
uint32_t chain_size_;
PACKAGE pkg;
int CommandLine_parse (int argc, char *argv[])
{
	chain_size_ = 0;
	threads = 1;
	pkg = MAFFT;
	int o;
	int index = 0;


	static struct option longOpts[] = {
		// {"version",	no_argument,	    0,			'v'},
		{"help",	no_argument,	    0,			'h'},
		{"in",      required_argument,  0,          'i'},
		{"out",     required_argument,  0,          'o'},
		{"thread",	required_argument,  0,			't'},
		{"package",	required_argument,  0,			'p'},
		{0, 0, 0, 0}
	
	};
	
	if (argc == 1)
  	{
    	Print_H();
    	return 0;
  	}

  	while( (o = getopt_long ( argc, argv, "h:i:o:t:p:", longOpts, &index)) != -1 )
  	{
		printf(" %c\n",o);
		switch (o)
		{
		case 'i':
			// inputfile = optarg;
			inputfile = (char*)malloc(SEQNAME_MAX_SIZE);
	  		strcpy(inputfile, optarg);
			break;
		case 'o':
			// outputfile = optarg;
			outputfile = (char*)malloc(SEQNAME_MAX_SIZE);
	  		strcpy(outputfile, optarg);
			break;
		case 't':
			threads = atoi(optarg);
			break;
		case 'p':
			if (!strcmp( optarg, "mafft")) pkg = MAFFT;
			else pkg = HALIGN;
			break;
		// case 's':
		// 	chain_size_ = atoi(optarg) - 1;
		// 	break;
		case 'h':
			Print_H();
			exit (0);
		case '?':
	  		fprintf(stderr, "Unavailable parameter: %s\n", longOpts[index].name);
			exit (1);
		default: 
			Print_H();
			exit (1);
		}
  	}
	
	return 0;

}


void Print_H()
{


	fprintf(stdout,"FMAlign: a fast and accurate multiple sequence aligner.\n\n");
	fprintf(stdout,"Usage: FMAlign [options]\n\n");


	fprintf(stdout,"General Options:\n");
	fprintf(stdout," -h\t\t\tShow the help file.\n");
	fprintf(stdout,"\n");


	fprintf(stdout,"Options of read mapping step:\n");
	fprintf(stdout," --in [inputfile]\tInput the sequences file in FASTA format\n");
	fprintf(stdout," --out [file]\t\tOutput of the aligned sequences in FASTA format (default: inputfile.fmalign).\n");
	fprintf(stdout," --thread [int]\t\tSet the number of CPU threads (default: 1).\n");
	fprintf(stdout," --package \t\tAligning the sequences with MSA tools. (default: mafft).\n");
	

	  
	fprintf(stdout,"\n\n");
  
}
