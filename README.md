**FMAlign**  
a fast multiple nucleotide sequence alignment method based on FM-index  

**How to use FMAlign**  
**Download:**  
>git clone https://github.com/iliuh/FMAlign.git  
>or download ZIP package  

**Build:**  
>cd FMAlign   
>make  

**Usage:**  
>./FMAlign [options]

* General Options:  
**--help**                  Show the help file.  

* Options of alignment step:  
**--in [filename]**         Input the sequences file in FASTA format  
**--out [file]**           Output of the aligned sequences in FASTA format. The default "output" is filename.fmalign.  
**--thread [int]**          Set the number of CPU threads (default: 1).  
**--package**               Aligning the sub-sequences by MSA tool, mafft or halign. The default "package" is mafft.  

**notice：** if package mafft can not work, please type
>cd packages/MAFFT
>make

**Related tools:**   
* SP is a jar file that was used in our experiment, it is used to compute the average value of sum-of-pair score for the alignment result. http://lab.malab.cn/soft/halign/  
* The FMtree library that was used in our work, FMtree is an improved FM-INDEX based on tree. https://github.com/chhylp123/FMtree


