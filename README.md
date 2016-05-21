<html>
#NOVOPlasty - The plastid assembler                            

Is currently being beta tested, you can always contact if you want already a copy

NOVOPlasty is a de novo assembler for short circular genomes.</br>
For the moment NOVOPlasty only supports the assembly of mitochondrial and plastid genomes with Illumina paired-end reads as input.

# Contact

Any issues/requests/problems/comments can be posted on [Github issues](https://github.com/ndierckx/NOVOPlasty/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens@hotmail.com 

<a href="http://ibsquare.be/" target="_blank"><img border="0" alt="IB2" src="http://ibsquare.be/sites/default/files/logo_1.png" width="100" height="100"></a>                


# Prerequisites

Perl


# Instructions

<strong>1\. Find a suitable seed</strong>

&nbsp;&nbsp;&nbsp;There are different types of seed possible:</br>
&nbsp;&nbsp;&nbsp;- A single read from the dataset that originates from the targeted plastid.</br>
&nbsp;&nbsp;&nbsp;- A plastid sequence derived from the same or a close related species.</br>
&nbsp;&nbsp;&nbsp;- A complete plastid sequence of a more distant species (recommended when there is no close related sequence &nbsp;&nbsp;&nbsp;&nbsp;available)

&nbsp;&nbsp;&nbsp;The format should be like a standard fasta file (first line: >Id_sequence)

<strong>2\. Create configuration file</strong>

&nbsp;&nbsp;&nbsp;You can download the example file (config.txt) and adjust the settings to your liking.</br>
&nbsp;&nbsp;&nbsp;Every parameter of the configuration file is explained below. 


<strong>3\. Run NOVOPlasty</strong>

&nbsp;&nbsp;&nbsp;No futher installation is necessary:

&nbsp;&nbsp;&nbsp;<code>perl NOVOPlasty.pl -c config.txt</code>

&nbsp;&nbsp;&nbsp;The input reads have to be uncompressed Illumina reads (fastq files).</br>
&nbsp;&nbsp;&nbsp;Either two seperate files(forward and reverse) or a merged fastq file.</br>
&nbsp;&nbsp;&nbsp;Multiple libraries as input is not yet supported.


<strong>4\. Output files</strong>

&nbsp;&nbsp;&nbsp;NOVOPlasty outputs four types of files:

&nbsp;&nbsp;&nbsp;1\. Circularized_assembly_projectname.fasta

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This file contains the assemblies that were succesfully circularized.

&nbsp;&nbsp;&nbsp;2\. Uncircularized_assembly_projectname.fasta

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This file contains the assemblies that were not circularized.


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
Genome Range         = 120000-160000
K-mer                = 38
Insert Range         = 1.45
Insert Range strict  = 1.2
Single/Paired        = PE
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
Type                 = (chloro/chloro2/mito) "chloro" for chloroplast assembly and "mito for mitochondrial assembly. It is     
                       recommended to first try "chloro", and if you don't get a circularized genome or it splits in to many
                       versions, to try "chloro2". This will create different contigs without circularizing
Genome Range         = (minimum genome size-maximum genome size) The expected genome size range of the genome.
                       Default value for mito: 12000-20000 / Default value for chloro: 120000-170000
                       If the expected size is know, you can lower the range, this can be useful when there is a repetitive
                       region, what could lead to a premature circularization of the genome.
K-mer                = (integer) This is the length of the overlap between matching reads (Default: 38). 
                       If reads are shorter then 90 bp, this value should be decreased. 
                       For reads longer then 101 bp, this value can be increased, but this is not necessary.
Insert Range         = This variation on the insert size, could lower it when the coverage is very high or raise it when the
                       coverage is too low (Default: 1.45). 
Insert Range strict  = Strict variation to resolve repetitive regions (Default: 1.2). 
Single/Paired        = For the moment only paired end reads are supported.
Extended log         = Prints out a very extensive log, could be useful to send me when there is a problem  (0/1).
Combined reads       = The path to the file that contains the combined reads (forward and reverse in 1 file)
Forward reads        = The path to the file that contains the forward reads
Reverse reads        = The path to the file that contains the reverse reads
Seed Input           = The path to the file that contains the seed sequence
</pre>
</html>
