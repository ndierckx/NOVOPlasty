<html>
#NOVOPlasty - The organelle assembler                            

NOVOPlasty is a de novo assembler for short circular genomes.</br>
For the moment NOVOPlasty only supports whole genome Illumina paired-end reads as input.

<strong>Last updates: 04/01/16 version 2.1</strong></br>
- Improved assembly of duplicated and repetitive regions.
- 'chloro' setting is replaced by the 'chloro2' setting (now the only options are 'chloro' and 'mito')</br>
<strong>11/11/16</strong>
- Improved seed retrieval with low coverage
- Data files without quality scores can be used as input file</br>
<strong>09/11/16</strong>
- IUPAC codes are correctly called, in stead of to many 'N's
- Bugs fixed
- Additional metrics are automatically calculated (Total reads, aligned reads, assembled reads, percentage organelle reads, average organelle coverage)
- Bugs in contig output fixed</br>


# Cite

<a href="http://nar.oxfordjournals.org/content/early/2016/10/24/nar.gkw955.full">Dierckxsens N., Mardulyn P. and Smits G. (2016) NOVOPlasty: De novo assembly of organelle genomes from whole genome data. <i>Nucleic Acids Research</i>, doi: 10.1093/nar/gkw955<a>


# Getting help

Any issues/requests/problems/comments can be posted on [Github issues](https://github.com/ndierckx/NOVOPlasty/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens@hotmail.com 

If your assembly was unsuccessful, you could already add the log file and the configuration file to the mail, this could help me to identify the problem!

<a href="http://ibsquare.be/" target="_blank"><img border="0" src="http://ibsquare.be/sites/default/files/logo_1.png" width="100" height="100" ></a>           


# Prerequisites

Perl


# Instructions

<strong>1\. Find a suitable seed</strong>

&nbsp;&nbsp;&nbsp;There are different types of seed possible:</br>
&nbsp;&nbsp;&nbsp;- A single read from the dataset that originates from the organelle plastid.</br>
&nbsp;&nbsp;&nbsp;- A organelle sequence derived from the same or a related species.</br>
&nbsp;&nbsp;&nbsp;- A complete organelle sequence of a more distant species (recommended when there is no close related</br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sequence available)

&nbsp;&nbsp;&nbsp;The format should be like a standard fasta file (first line: >Id_sequence)

&nbsp;&nbsp;&nbsp;Be cautious for seed sequences that are similar in both mitochondrial and chloroplast genomes.</br>
&nbsp;&nbsp;&nbsp;We observed good results with RUBP sequences as seeds for chloroplast assembly

<strong>2\. Create configuration file</strong>

&nbsp;&nbsp;&nbsp;You can download the example file (config.txt) and adjust the settings to your liking.</br>
&nbsp;&nbsp;&nbsp;Every parameter of the configuration file is explained below. 


<strong>3\. Run NOVOPlasty</strong>

&nbsp;&nbsp;&nbsp;No further installation is necessary:

&nbsp;&nbsp;&nbsp;<code>perl NOVOPlasty.pl -c config.txt</code>

&nbsp;&nbsp;&nbsp;The input reads have to be uncompressed Illumina reads (fastq/fasta files).</br>
&nbsp;&nbsp;&nbsp;Either two separate files(forward and reverse) or a merged fastq/fasta file.</br>
&nbsp;&nbsp;&nbsp;Multiple libraries as input is not yet supported.

&nbsp;&nbsp;&nbsp;DO NOT filter or quality trim the reads!!! Use the raw whole genome dataset!</br>

&nbsp;&nbsp;&nbsp;You can subsample to speed up the process and to reduce the memory requirements. But it is recommended </br> &nbsp;&nbsp;&nbsp;to use as much reads as possible, especially when the organelle genome contains AT-rich stretches.

&nbsp;&nbsp;&nbsp;Recommended maximum K-mer lengths:</br>

&nbsp;&nbsp;&nbsp;100 bp reads: +/- 39</br>
&nbsp;&nbsp;&nbsp;150 bp reads: +/- 49</br>
&nbsp;&nbsp;&nbsp;250 bp reads: +/- 73</br>

&nbsp;&nbsp;&nbsp;You can always try different K-mer's. In the case of low coverage problems or seed errors, </br> 
&nbsp;&nbsp;&nbsp;it's recommended to lower the K-mer (set to 39)!!!.</br>

&nbsp;&nbsp;&nbsp;<strong>4\. Output files</strong>

&nbsp;&nbsp;&nbsp;NOVOPlasty outputs four types of files:

&nbsp;&nbsp;&nbsp;1\. Contigs_projectname.txt

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; This file contains the contigs of the assemblies.</br> 

&nbsp;&nbsp;&nbsp;2\. Merged_contigs_projectname.txt

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;When there are multiple contigs, NOVOPlasty will try to combine all contigs </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;in to a complete circular genome, all the different possibilities can be found in this file.</br> 

&nbsp;&nbsp;&nbsp;3\. Option_nr_projectname.txt

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;All possible contig combinations will have a seperate fasta file. </br>

&nbsp;&nbsp;&nbsp;4\. contigs_tmp_projectname.txt

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If non of the above files are outputted or are empty, you can retrieve some contigs from this file.

&nbsp;&nbsp;&nbsp;<strong>5\. Interpretation</strong>

&nbsp;&nbsp;&nbsp;1\. Chloroplast assembly

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Ideally you will have one or two outputted assemblies. When you have two assemblies from the same length,</br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the only difference will be the orientation of the inverted repeat. This can be resolved manually by mapping</br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the assemblies to the closest reference. (On NCBI's BLAST you can further examine your mapping by clicking on</br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'Graphics', this will show you which orientation is correct.) </br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Otherwise you can first annotate the two assemblies and compare the gene order.

# Configuration file

This is an example of a configuration file for the assembly of a chloroplast.
To make the assembler work, your configuration file has to have the exact same structure.
(Make sure there is always a space after the equals sign and every parameter is captured in one single line)

<strong>1\. Example of configuration file:</strong>
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

<strong>2\. Explanation parameters:</strong>
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
