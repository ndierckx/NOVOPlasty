<html>
#NOVOPlasty - The plastid assembler 

Will be available by then end of 2015

NOVOPlasty is a de novo assembler for short circular genomes.</br>
For the moment NOVOPlasty only supports the assembly of plastid genomes with Illumina paired-end reads as input.

# Contact

For any questions, you can contact me directly through the following email address:

nicolasdierckxsens@hotmail.com 


# Prerequisites

Perl


# Instructions

<strong>1\. Find a suitable seed</strong>

&nbsp;&nbsp;&nbsp;There are different types of seed possible:</br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- A single read from the dataset that originates from the targeted plastid.</br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- A plastid sequence derived from the same or a close related species.</br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- A complete plastid sequence (recommended when there is no close related sequence available)

&nbsp;&nbsp;&nbsp;The format should be like a standard fasta file (first line: >Id_sequence)

<strong>2\. Create configuration file</strong>

&nbsp;&nbsp;&nbsp;You can download the example file (config.txt) and adjust the settings to your liking.</br>
&nbsp;&nbsp;&nbsp;Every parameter of the configuration file is explained below. 


<strong>3\. Run NOVOPlasty</strong>

&nbsp;&nbsp;&nbsp;Make sure all the files are in the same folder.</br>
&nbsp;&nbsp;&nbsp;No futher installation is necessary:

&nbsp;&nbsp;&nbsp;<code>perl NOVOPlasty.pl -c config.txt</code>


<strong>4\. Output files</strong>

&nbsp;&nbsp;&nbsp;NOVOPlasty outputs two different files:

&nbsp;&nbsp;&nbsp;1\. Circularized_assemblies_projectname.fasta

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This file contains the assemblies that were succesfully circularized.

&nbsp;&nbsp;&nbsp;2\. Uncircularized_assemblies_projectname.fasta

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This file contains the assemblies that were not circularized.


# Configuration file

This is an example of a configuration file for the assembly of a chloroplast.
To make the assembler work, your configuration file has to have the exact same structure.
(Make sure there is always a space after the equals sign and every parameter is captured in one single line)

<strong>1\. Example of config file:</strong>
<pre>
Project name     = AOB_chloro
Insert size      = 300
Insert size auto = yes
Read Length      = 101
Type             = chloro
Genome Range     = 120000-160000
K-mer            = 38
Combined reads   = /path/to/reads/AOB_reads.fastq
Forward reads    = 
Reverse reads    = 
Seed Input       = Seed_AOB.fasta
</pre>

<strong>2\. Explanation parameters:</strong>
<pre>
Project name     = Choose a name for your project, it will be used for the output files.
Insert size      = Total insert size of your paired end reads, it doesn't have to be accurate but should be close enough.
Insert size auto = (yes/no) This will finetune your insert size automatically (Default: yes)
Read Length      = The read length of your reads.
Type             = (chloro/mito) "chloro" for chloroplast assembly and "mito for mitochondrial assembly
Genome Range     = (minimum genome size-maximum genome size) The expected genome size range of the genome.
K-mer            = 38
Combined reads   = The path to the file that contains the combined reads (forward and reverse in 1 file)
Forward reads    = The path to the file that contains the forward reads
Reverse reads    = The path to the file that contains the reverse reads
Seed Input       = The path to the file that contains the seed sequence
</pre>
</html>
