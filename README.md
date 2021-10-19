# NOVOPlasty - The organelle assembler and heteroplasmy caller

NOVOPlasty is a de novo assembler and heteroplasmy/variance caller for short circular genomes.  

**Last updates: 04/02/21 version 4.3**  
- UPDATED CONFIG FILE! 
- Added an output directory option to the config file
- Added an store hash locally option to the config file to speed up the assembly when using the same dataset for multiple runs               
**03/06/20**                                                                                                                                     
- Improved assembly method                                                                       
- Improved heteroplasmy detection.                                                               
- Heteroplasmy output now includes AF and DP for indels.                                                        
**07/01/20**                                                                     
- UPDATED CONFIG FILE!                                                          
- Improved heteroplasmy calling
- Improved organelle assembly
- NEW FUNCTION: It is now possible to extend a seed directly, this will not correct the given seed, so it may not have mismatches                                                                                                    
**24/04/19** 
- UPDATED CONFIG FILE!                                                                     
- New heteroplasmy function (can link heteroplasmic mutations by local assembly)
- New batch function   
- Small improvements in assembly                                                                                         
                                                                                                                                                                                                                                                         
## Cite

<a href="http://nar.oxfordjournals.org/content/early/2016/10/24/nar.gkw955.full">Dierckxsens N., Mardulyn P. and Smits G. (2016) NOVOPlasty: De novo assembly of organelle genomes from whole genome data. <i>Nucleic Acids Research</i>, doi: 10.1093/nar/gkw955<a>
 
 If you use the heteroplasmy or variance calling option, please cite both articles:
 
  <a href="https://academic.oup.com/nargab/article/doi/10.1093/nargab/lqz011/5606586?guestAccessKey=46546992-7b07-4a28-82f4-89baea9527f5">Dierckxsens N., Mardulyn P. and Smits G. (2019) Unraveling heteroplasmy patterns with NOVOPlasty. <i>NAR Genomics and Bioinformatics</i>, https://doi.org/10.1093/nargab/lqz011<a>


## Getting help

Any issues/requests/problems/comments that are not yet addressed on this page can be posted on [Github issues](https://github.com/ndierckx/NOVOPlasty/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens at hotmail dot com 

If your assembly was unsuccessful, you could already add the log file and the configuration file to the mail, this could help me to identify the problem!

<a href="http://www.huderf.be/" target="_blank"><img border="0" src="http://cdn.prezly.com/08/500080108e11e68f98cbaed0afad97/HUDERF-LOGO_cmyb_4_www.jpg" width=auto height="95" ></a>
<a href="http://ebe.ulb.ac.be/ebe/ebe-Welcome.html" target="_blank"><img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/9/99/Universit%C3%A9_libre_de_Bruxelles_logo.svg/1024px-Universit%C3%A9_libre_de_Bruxelles_logo.svg.png" width=auto height="95" ></a>
<a href="http://www.belgiankidsfund.be" target="_blank"><img border="0" src="http://www.belgiankidsfund.be/images/logo.png" width=auto height="95" ></a>  


## Instructions

**!A more complete manual can found under the [Wiki](https://github.com/ndierckx/NOVOPlasty/wiki) section!**

### 1. Find a suitable seed

There are different types of seed possible:
- A single read from the dataset that originates from the organelle genome.
- A organelle sequence derived from the same or a related species.
- A complete organelle sequence of a more distant species (recommended when there is no close related sequence available)

The format should be like a standard fasta file (first line: >Id_sequence)

Be cautious for seed sequences that are similar in both mitochondrial and chloroplast genomes.  
We observed good results with RUBP sequences as seeds for chloroplast assembly.

### 2. Create configuration file

You can download the example file (config.txt) and adjust the settings to your liking.  
Every parameter of the configuration file is explained below. 


### 3. Run NOVOPlasty

No further installation is necessary:

<code>perl NOVOPlasty4.3.pl -c config.txt</code>
    
To run with conda:
    
<code>NOVOPlasty4.2.pl -c config.txt</code>

The input reads have to be uncompressed Illumina reads (fastq/fasta files) or gz/bz2 zipped files.  
There is also an Ion Torrent option, but it does not produce the best results. 
Either two separate files(forward and reverse) or a merged fastq/fasta file.  
Multiple libraries as input is not yet supported.

DO NOT filter or quality trim the reads!!! Use the raw whole genome dataset (Only adapters should be removed)!

You can subsample to speed up the process and to reduce the memory requirements. This also possible by using the max memory option in the config file. But it is recommended to use as much reads as possible, especially when the organelle genome contains AT-rich stretches.

You can always try different K-mer's. In the case of low coverage problems or seed errors, it's recommended to lower the K-mer (set between 21-39)!!!.


----------------------------------------------------------------------------------------------------------

## Configuration file

This is an example of a configuration file for the assembly of a chloroplast.
To make the assembler work, your configuration file has to have the exact same structure.
(Make sure there is always a space after the equals sign and every parameter is captured in one single line)

**1. Example of configuration file:**
<pre>

Project:
-----------------------
Project name          = Test
Type                  = mito
Genome Range          = 12000-22000
K-mer                 = 33
Max memory            = 
Extended log          = 0
Save assembled reads  = no
Seed Input            = /path/to/seed_file/Seed.fasta
Extend seed directly  = no
Reference sequence    = /path/to/reference_file/reference.fasta (optional)
Variance detection    = 
Chloroplast sequence  = /path/to/chloroplast_file/chloroplast.fasta (only for "mito_plant" option)

Dataset 1:
-----------------------
Read Length           = 151
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        = 
Forward reads         = /path/to/reads/reads_1.fastq
Reverse reads         = /path/to/reads/reads_2.fastq
Store Hash            =

Heteroplasmy:
-----------------------
MAF                   = 
HP exclude list       = 
PCR-free              = 

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
Output path           = 

</pre>

**2. Explanation parameters:**
<pre>

Project:
-----------------------
Project name         = Choose a name for your project, it will be used for the output files.
Type                 = (chloro/mito/mito_plant) "chloro" for chloroplast assembly, "mito" for mitochondrial assembly and 
                       "mito_plant" for mitochondrial assembly in plants.
Genome Range         = (minimum genome size-maximum genome size) The expected genome size range of the genome.
                       Default value for mito: 12000-20000 / Default value for chloro: 120000-200000
                       If the expected size is know, you can lower the range, this can be useful when there is a repetitive
                       region, what could lead to a premature circularization of the genome.
K-mer                = (integer) This is the length of the overlap between matching reads (Default: 33). 
                       If reads are shorter then 90 bp or you have low coverage data, this value should be decreased down to 23. 
                       For reads longer then 101 bp, this value can be increased, but this is not necessary.
Max memory           = You can choose a max memory usage, suitable to automatically subsample the data or when you have limited                      
                       memory capacity. If you have sufficient memory, leave it blank, else write your available memory in GB
                       (if you have for example a 8 GB RAM laptop, put down 7 or 7.5 (don't add the unit in the config file))
Extended log         = Prints out a very extensive log, could be useful to send me when there is a problem  (0/1).
Save assembled reads = All the reads used for the assembly will be stored in seperate files (yes/no)
Seed Input           = The path to the file that contains the seed sequence.
Extend seed directly = This gives the option to extend the seed directly, in stead of finding matching reads. Only use this when your seed 
                       originates from the same sample and there are no possible mismatches (yes/no)
Reference (optional) = If a reference is available, you can give here the path to the fasta file.
                       The assembly will still be de novo, but references of the same genus can be used as a guide to resolve 
                       duplicated regions in the plant mitochondria or the inverted repeat in the chloroplast. 
                       References from different genus haven't beeen tested yet.
Variance detection   = If you select yes, you should also have a reference sequence (previous line). It will create a vcf file                
                       with all the variances compared to the give reference (yes/no)
Chloroplast sequence = The path to the file that contains the chloroplast sequence (Only for mito_plant mode).
                       You have to assemble the chloroplast before you assemble the mitochondria of plants!

Dataset 1:
-----------------------
Read Length          = The read length of your reads.
Insert size          = Total insert size of your paired end reads, it doesn't have to be accurate but should be close enough.
Platform             = illumina/ion - The performance on Ion Torrent data is significantly lower
Single/Paired        = For the moment only paired end reads are supported.
Combined reads       = The path to the file that contains the combined reads (forward and reverse in 1 file)
Forward reads        = The path to the file that contains the forward reads (not necessary when there is a merged file)
Reverse reads        = The path to the file that contains the reverse reads (not necessary when there is a merged file)
Store Hash           = If you want several runs on one dataset, you can store the hash locally to speed up the process (put "yes" to store the hashes locally)
                       To run local saved files, goto te wiki section of the github page

Heteroplasmy:
-----------------------
MAF                  = (0.007-0.49) Minor Allele Frequency: If you want to detect heteroplasmy, first assemble the genome without this option. Then give the resulting                         
                       sequence as a reference and as a seed input. And give the minimum minor allele frequency for this option 
                       (0.01 will detect heteroplasmy of >1%)
HP exclude list      = Option not yet available  
PCR-free             = (yes/no) If you have a PCR-free library write yes

Optional:
-----------------------
Insert size auto     = (yes/no) This will finetune your insert size automatically (Default: yes)                               
Use Quality Scores   = It will take in account the quality scores, only use this when reads have low quality, like with the    
                       300 bp reads of Illumina (yes/no)
Output path          = You can change the directory where all the output files wil be stored.

</pre>
</html>
