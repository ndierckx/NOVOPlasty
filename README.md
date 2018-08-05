# NOVOPlasty - The organelle assembler

NOVOPlasty is a de novo assembler and variance caller for short circular genomes.  
For the moment NOVOPlasty only supports whole genome Illumina paired-end reads as input.

**Last updates: 04/08/18 version 2.7.2**  
- A bug in the variance caller is resolved
- Added new read ids                                                                       
**07/06/18**                                                                                      
- There was a bug since version 2.2.3 that made seed files, where the bases weren't in CAPs (actg in stead of ACTG), unreadable    and outputted INVALID SEED ERROR.
- A new option was added to the config file (Use Quality Scores), Use this option for the 300 bp Illumina reads (as long as their high error rates aren't fixed) UPDATE CONFIG FILE!                                                                      
**24/05/18**                                                                                                  
- Version 2.6.8 and 2.6.9 had a bug that read the discriptions of the config file in stead of your input                  
**23/05/18**                                                                                                         
- Improved assembly in SNR regions                                                              
**19/05/18**                                                                                            
- Forgot to delete a line of code that prevented to use more than 15 GB of RAM, use version 2.6.8. or higher! 
- Problems with incorrect reading of the config file (especially with Windows) should be resolved
- UPDATED CONFIG FILE!   
**16/04/18**                                                                                            
- Read ids that end with '1:N:0:1' lead to errenous results with versions 2.6.5 and 2.6.6                                  
**11/03/18**                                                                                       
- Automatic insert size was inaccurate for version 2.6.4 (fixed)
- Resolved some bugs and small improvements.  
- Improved Heteroplamsy calling                                                          
**20/02/18**                                                                                     
- UPDATED CONFIG FILE!                                                                                  
- There are two new options, a basic variance and heteroplasmy caller
- Improved reference guidance
- Resolved some bugs and small improvements.                                                                                                                                                                                                                                                                                                                                                                          

## Cite

<a href="http://nar.oxfordjournals.org/content/early/2016/10/24/nar.gkw955.full">Dierckxsens N., Mardulyn P. and Smits G. (2016) NOVOPlasty: De novo assembly of organelle genomes from whole genome data. <i>Nucleic Acids Research</i>, doi: 10.1093/nar/gkw955<a>


## Getting help

Any issues/requests/problems/comments that are not yet addressed on this page can be posted on [Github issues](https://github.com/ndierckx/NOVOPlasty/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens at hotmail dot com 

If your assembly was unsuccessful, you could already add the log file and the configuration file to the mail, this could help me to identify the problem!

<a href="http://ibsquare.be/" target="_blank"><img border="0" src="http://ibsquare.be/sites/default/files/logo_1.png" width="100" height="100" ></a>
<a href="http://www.huderf.be/" target="_blank"><img border="0" src="http://cdn.prezly.com/08/500080108e11e68f98cbaed0afad97/HUDERF-LOGO_cmyb_4_www.jpg" width=auto height="95" ></a>
<a href="http://ebe.ulb.ac.be/ebe/ebe-Welcome.html" target="_blank"><img border="0" src="https://www.ulb.ac.be/dre/com/docs/logo3lg.jpg" width=auto height="95" ></a>


## Prerequisites

Perl


## Instructions

### 1. Find a suitable seed

There are different types of seed possible:
- A single read from the dataset that originates from the organelle plastid.
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

<code>perl NOVOPlasty.pl -c config.txt</code>

The input reads have to be uncompressed Illumina reads (fastq/fasta files) or gz zipped files.  
Either two separate files(forward and reverse) or a merged fastq/fasta file.  
Multiple libraries as input is not yet supported.

DO NOT filter or quality trim the reads!!! Use the raw whole genome dataset (Only adapters should be removed)!

You can subsample to speed up the process and to reduce the memory requirements. But it is recommended to use as much reads as possible, especially when the organelle genome contains AT-rich stretches.

You can always try different K-mer's. In the case of low coverage problems or seed errors, it's recommended to lower the K-mer (set between 25-39)!!!.

### 4. Output files

NOVOPlasty outputs four types of files:

#### 1. Contigs_projectname.txt

This file contains the contigs of the assemblies.

#### 2. Merged_contigs_projectname.txt

When there are multiple contigs, NOVOPlasty will try to combine all contigs in to a complete circular genome, all the different possibilities can be found in this file.

#### 3. Option_nr_projectname.txt

All possible contig combinations will have a seperate fasta file.

#### 4. contigs_tmp_projectname.txt

If non of the above files are outputted or are empty, you can retrieve some contigs from this file.

----------------------------------------------------------------------------------------------------------
### 5. Heteroplasmy detection

Only use heteropalsmy detection when you have enough coverage (>100X)!

#### 1. Assemble the organelle genome

First assemble the organelle genome with NOVOPlasty but without the heteroplasmy function.

#### 2. Prepare the reference sequence

For heteroplasmy detection, you have to use the same sequence as reference AND as seed!
This sequence should always be from the assembled genome in the previous step.
If there are no repetitive regions, you can use the complete sequence, otherwise you should remove these regions.

#### 3. Minimum minor allele frequency

Only potential mutations above this frequency will be detected (Value for the "Heteroplasmy" option in the config file)
Heteroplasmy detection will be activated by giving a value for this option in the config file.

A value of 0.01 will detect heteroplasmy above 1%, a value of 0.2 will detect heteropalsmy above 20%.
If you want to detect low frequencies around 1%, you should not go lower than 0.007. And only go this low when you have very high coverage (>1000X). The higher the coverage the more accurate!

&nbsp;
## Interpretation and post-processing

### 1. General

#### *  
A '*' in the fasta output files indicates that the nucleotide before is a possible deletion/insertion. This can occur when the exact length of single nucleotide repeat can't be determined exactly due to systemic Illumina sequencing errors or within repetitive regions. Since this sign can interfere with post processing algorithms it is best resolve them manually or to delete them. 

#### Gaps  
Most gaps are caused by Single Nucleotide Repeats (SNR). Illumina seqeuncers have a high rate of systemic errors after SNR's and are therefore hard to assemble. NOVOPlasty is cabale of assembling these regions as correct as possible by approaching these regions from both sides (sequencing errors commence once in the SNR). If this region is not too long, NOVOPlasty can automaticallly merge both sides, otherwise it will output a gap. Although this gap can often be closed automatically (if both sides overlap).  
You can find an example (form the Avicennia officinalis chloroplast assembly) below how to do this manually:

These regions are indicated by 15 N's:
```
TTCTTGTCATTTCTCCCCCCCCCCCCCBTTTTTTTTTTHAAAAAAAAAAAANNNNNNNNNNNNNNNTTTTTTTCCTTTCCCCCCCCCCCCCCCTTTTTTTTTTCAAAAAAAAAAGAGACGAGAAACTC
```
Remove the N's and align both sides (Remember that the end of the first sequence and the start of the second sequence are not reliable):
```
TTCTTGTCATTTCTCCCCCCCCCCCCCBTTTTTTTTTTHAAAAAAAAAAAA
 TTTTTTTCCTTTCCCCCCCCCCCCCCCTTTTTTTTTTCAAAAAAAAAAGAGACGAGAAACTCTGAA
```
The most likable consensus sequence would be this, so you should correct the assembly like this:
```
TTCTTGTCATTTCTCCCCCCCCCCCCCCTTTTTTTTTTCAAAAAAAAAAGAGACGAGAAACTCTGAA
```
It is adviced to make these corrections before you verify the assembly by realigning the reads

### 2. Chloroplast assembly

Ideally you will have one or two outputted assemblies. When you have two assemblies from the same length, the only difference will be the orientation of the inverted repeat. This can be resolved manually by mapping the assemblies to the closest reference. (On NCBI's BLAST you can further examine your mapping by clicking on 'Graphics', this will show you which orientation is correct.) Otherwise you can first annotate the two assemblies and compare the gene order.


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
K-mer                 = 39
Max memory            = 
Extended log          = 0
Save assembled reads  = no
Seed Input            = Seed.fasta
Reference sequence    = /path/to/reference_file/reference.fasta (optional)
Variance detection    = no
Heteroplasmy          = 
HP exclude list       =
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

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.8
Insert Range strict   = 1.3
Use Quality Scores    = no
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
K-mer                = (integer) This is the length of the overlap between matching reads (Default: 39). 
                       If reads are shorter then 90 bp or you have low coverage data, this value should be decreased down to 23. 
                       For reads longer then 101 bp, this value can be increased, but this is not necessary.
Max memory           = You can choose a max memory usage, suitable to automatically subsample the data or when you have limited                      
                       memory capacity. If you have sufficient memory, leave it blank, else write your available memory in GB
                       (if you have for example a 8 GB RAM laptop, put down 7 or 7.5 (don't add the unit in the config file))
Extended log         = Prints out a very extensive log, could be useful to send me when there is a problem  (0/1).
Save assembled reads = All the reads used for the assembly will be stored in seperate files (yes/no)
Seed Input           = The path to the file that contains the seed sequence.
Reference (optional) = If a reference is available, you can give here the path to the fasta file.
                       The assembly will still be de novo, but references of the same genus can be used as a guide to resolve 
                       duplicated regions in the plant mitochondria or the inverted repeat in the chloroplast. 
                       References from different genus haven't beeen tested yet.
Variance detection   = If you select yes, you should also have a reference sequence (previous line). It will create a vcf file                
                       with all the variances compared to the give reference (yes/no)
Heteroplasmy         = If you want to detect heteroplasmy,first assemble the genome without this option. Then give the resulting                         
                       sequence as a reference and as a seed input. And give the minimum minor allele frequency for this option 
                       (0.01 will detect heteroplasmy of >1%)
HP exclude list      = Option not yet available                      
Chloroplast sequence = The path to the file that contains the chloroplast sequence (Only for mito_plant mode).
                       You have to assemble the chloroplast before you assemble the mitochondria of plants!

Dataset 1:
-----------------------
Read Length          = The read length of your reads.
Insert size          = Total insert size of your paired end reads, it doesn't have to be accurate but should be close enough.
Platform             = illumina is for now the only option
Single/Paired        = For the moment only paired end reads are supported.
Combined reads       = The path to the file that contains the combined reads (forward and reverse in 1 file)
Forward reads        = The path to the file that contains the forward reads (not necessary when there is a merged file)
Reverse reads        = The path to the file that contains the reverse reads (not necessary when there is a merged file)

Optional:
-----------------------
Insert size auto     = (yes/no) This will finetune your insert size automatically (Default: yes)
Insert Range         = This variation on the insert size, could lower it when the coverage is very high or raise it when the
                       coverage is too low (Default: 1.6). 
Insert Range strict  = Strict variation to resolve repetitive regions (Default: 1.2).          
Use Quality Scores   = It will take in account the quality scores, only use this when reads have low quality, like with the    
                       300 bp reads of Illumina (yes/no)
</pre>
</html>
