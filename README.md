#NOVOPlasty - The organelle assembler                            

NOVOPlasty is a de novo assembler for short circular genomes.  
For the moment NOVOPlasty only supports whole genome Illumina paired-end reads as input.

**Last updates: 25/01/17 version 2.3**
- SEED RETRIEVAL bug is fixed, if the assembly got stuck at seed retrieval use version after 2.3!!!
- More accurate assembly of erroneous single nucleotide repeat regions
**24/01/17**  
- Fixed several bugs
- Reduced cross assembly of choroplast and plant mitochondrial sequences
**04/01/17**  
- Improved assembly of duplicated and repetitive regions.
- 'chloro' setting is replaced by the 'chloro2' setting (now the only options are 'chloro' and 'mito')  
**11/11/16**  
- Improved seed retrieval with low coverage
- Data files without quality scores can be used as input file  
**09/11/16**  
- IUPAC codes are correctly called, in stead of to many 'N's
- Bugs fixed
- Additional metrics are automatically calculated (Total reads, aligned reads, assembled reads, percentage organelle reads, average organelle coverage)
- Bugs in contig output fixed


## Cite

<a href="http://nar.oxfordjournals.org/content/early/2016/10/24/nar.gkw955.full">Dierckxsens N., Mardulyn P. and Smits G. (2016) NOVOPlasty: De novo assembly of organelle genomes from whole genome data. <i>Nucleic Acids Research</i>, doi: 10.1093/nar/gkw955<a>


## Getting help

Any issues/requests/problems/comments that are not yet addressed on this page can be posted on [Github issues](https://github.com/ndierckx/NOVOPlasty/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens@hotmail.com 

If your assembly was unsuccessful, you could already add the log file and the configuration file to the mail, this could help me to identify the problem!

<a href="http://ibsquare.be/" target="_blank"><img border="0" src="http://ibsquare.be/sites/default/files/logo_1.png" width="100" height="100" ></a>           


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

The input reads have to be uncompressed Illumina reads (fastq/fasta files).  
Either two separate files(forward and reverse) or a merged fastq/fasta file.  
Multiple libraries as input is not yet supported.

DO NOT filter or quality trim the reads!!! Use the raw whole genome dataset (Only adapters should be removed)!

You can subsample to speed up the process and to reduce the memory requirements. But it is recommended to use as much reads as possible, especially when the organelle genome contains AT-rich stretches.

You can always try different K-mer's. In the case of low coverage problems or seed errors, it's recommended to lower the K-mer (set to 39)!!!.

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

## Interpretation and post-processing

### 1. General

#### *  
A '*' in the fasta output files indicates that the nucleotide before is a possible deletion/insertion. This can occur when the exact length of single nucleotide repeat can't be determined exactly due to systemic Illumina sequencing errors. Since this sign can interfere with post processing algorithms it is best resolve them manually or to delete them. 

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
Project name         = AOB_chloro
Insert size          = 300
Insert size aut      = yes
Read Length          = 101
Type                 = chloro
Genome Range         = 120000-200000
K-mer                = 39
Insert Range         = 1.5
Insert Range strict  = 1.2
Single/Paired        = PE
Coverage Cut off     = 1000
Extended log         = 0
Combined reads       = /path/to/reads/AOB_reads.fastq
Forward reads        = 
Reverse reads        = 
Seed Input           = Seed_AOB.fasta
</pre>

**2. Explanation parameters:**
<pre>
Project name         = Choose a name for your project, it will be used for the output files.
Insert size          = Total insert size of your paired end reads, it doesn't have to be accurate but should be close enough.
Insert size auto     = (yes/no) This will finetune your insert size automatically (Default: yes)
Read Length          = The read length of your reads.
Type                 = (chloro/mito) "chloro" for chloroplast assembly and "mito for mitochondrial assembly.
Genome Range         = (minimum genome size-maximum genome size) The expected genome size range of the genome.
                       Default value for mito: 12000-20000 / Default value for chloro: 120000-200000
                       If the expected size is know, you can lower the range, this can be useful when there is a repetitive
                       region, what could lead to a premature circularization of the genome.
K-mer                = (integer) This is the length of the overlap between matching reads (Default: 39). 
                       If reads are shorter then 90 bp, this value should be decreased. 
                       For reads longer then 101 bp, this value can be increased, but this is not necessary.
Insert Range         = This variation on the insert size, could lower it when the coverage is very high or raise it when the
                       coverage is too low (Default: 1.5). 
Insert Range strict  = Strict variation to resolve repetitive regions (Default: 1.2). 
Single/Paired        = For the moment only paired end reads are supported.
Coverage Cut off     = You can speed up the assembly by lowering the coverage cut off, standard it will use up to 1000 coverage
Extended log         = Prints out a very extensive log, could be useful to send me when there is a problem  (0/1).
Combined reads       = The path to the file that contains the combined reads (forward and reverse in 1 file)
Forward reads        = The path to the file that contains the forward reads (not necessary when there is a merged file)
Reverse reads        = The path to the file that contains the reverse reads (not necessary when there is a merged file)
Seed Input           = The path to the file that contains the seed sequence
</pre>
</html>
