#include <stdio.h>
#include <string.h>
#include <cstring>
#include <stdlib.h>
#include <limits.h>
#include "minimizer.h"
#include "Auxiliary.h"






void convert_unint_to_chars(unsigned int input, char* output, unsigned int length_of_chars)
{
  unsigned int i = 0, j = 0;
  char ACGT_C[4];

  ACGT_C[0] = 'A';
  ACGT_C[1] = 'C';
  ACGT_C[2] = 'G';
  ACGT_C[3] = 'T';


  output[length_of_chars] = '\0';

  j = length_of_chars - 1;

  for (i = 0; i <length_of_chars; i++)
  {
    output[j] = ACGT_C[input % 4];
    input = input / 4;
    j--;
  }
}





unsigned int *contain(char *str, unsigned int len, unsigned int k)
{
	// fprintf(stderr, "%s\n", str);
  	unsigned int ctoi[256];

  	ctoi['A'] = 0;
  	ctoi['C'] = 1;
  	ctoi['G'] = 2;
  	ctoi['T'] = 3;
    ctoi['a'] = 0;
    ctoi['c'] = 1;
    ctoi['g'] = 2;
    ctoi['t'] = 3;

    unsigned int *output;
    // unsigned int len = strlen(str);
    output = (unsigned int*)malloc(sizeof(unsigned int) * (len-k+1));


    if(len < k)return 0;
   
  
    unsigned int strh = 0;
    unsigned int B = 4;
    unsigned int t = 1;


    for (int i = 0; i < k; ++i) t*=B;

    for (int i = 0; i < k; ++i)
    {
        strh = strh * B + ctoi[ str[i] ];
        output[0] = strh;
    }
    // fprintf(stderr, "%u\n", strh);

    for (int i = 0; i + k < len; ++i)
    {
        strh = strh * B + ctoi[ str[i + k] ] - ctoi[ str[i] ] * t;
        output[i + 1] = strh;
        // fprintf(stderr, "%u\n", strh);
    }

    return output;
    // free(output);
}

unsigned int minimize(char *read, unsigned int k, unsigned int *mm)
{
  
  	// unsigned int *wnd;
  	// puts(read);

  	unsigned int min_seed = INT_MAX;
  	unsigned int flag_min;

  	unsigned int readlen = strlen(read);


  	unsigned int *seed;
  	unsigned int seedlen;
  	seedlen = readlen - k + 1;

    seed = (unsigned int*)malloc(sizeof(unsigned int) * seedlen);
  	seed = contain(read, readlen, k);          //这里将所有kmer处理成了无符号整数


  	unsigned int wnd_size;
  	wnd_size  = 7;


  	// mm = (unsigned int *)malloc(sizeof(unsigned int) * seedlen);//记录了所有无符号整形的minimizer

  
  	unsigned int id = 0;


  	for (int i = 0; i < wnd_size; ++i)
  	{
    	if (seed[i] < min_seed) 
    	{
      		min_seed = seed[i];
      		flag_min = i;
    	}
  	}

  	mm[id] = min_seed;
  	++ id;
  	// fprintf(stderr, "--%u\n", min_seed);


  	for (int i = 0; i + wnd_size < seedlen  ; ++i)
  	{
    	if (flag_min == i)
    	{
      		min_seed = INT_MAX;
      		for (int j = 1; j <= wnd_size; ++j)
      		{
        		if (seed[i] < min_seed) 
        		{
          			min_seed = seed[i + j];
          			flag_min = i + j;
        		}
      		}

      		mm[id] = min_seed;
      		++ id;
   		}

    	else
    	{
      		if (min_seed > seed[i + wnd_size])
      		{
        		min_seed = seed[i + wnd_size];
        		flag_min = i + wnd_size;
        		mm[id] = min_seed;
        		++ id;
        		// fprintf(stderr, "--%u\n", min_seed);
      		}
      		else
      		{
        		continue;
      		}
    	}
  	}

  	fprintf(stdout, "minimizer ----------------- %u\n", id);
  	free (seed);
    return id;
}


