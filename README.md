<html>
#NOVOPlasty - The plastid assembler 

Will be available by then end of 2015

NOVOPlasty is a de novo assembler for short circular genomes.</br>
For the moment NOVOPlasty supports only the assembly of plastid genomes and with Illumina paired-end reads as input.

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

<pre>
Project name     = AOB_chloro
Insert size      = 300
Insert size auto = yes
Read Length      = 101
Type             = chloro
Encrypt          = no
Reference        = no
Mito Range       = 12000-20000
Chloro Range     = 120000-160000
K-mer            = 38
Combined reads   = /path/to/reads/AOB_reads.fastq
Forward reads    = 
Reverse reads    = 
Seed Input       = Seed_AOB.fasta
</pre>
</html>


