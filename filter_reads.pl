#!/usr/bin/env perl
######################################################
#         SOFTWARE COPYRIGHT NOTICE AGREEMENT        #
#  Copyright (C) {2019}  {Nicolas Dierckxsens}  #
#              All Rights Reserved                   #
#         See file LICENSE for details.              #
######################################################
#           nicolasdierckxsens@hotmail.com
use strict;
use Getopt::Long;

my $project = "";
my $reads12 = "";
my $reads1 = "";
my $reads2 = "";
my $reference = "";
my $kmer = '16';

GetOptions (
            "1=s" => \$reads1,
            "2=s" => \$reads2,
            "ref=s" => \$reference,
            "out=s" => \$project,
            ) or die "Incorrect usage!\n";

my $USAGE = "\nUsage: perl filter_reads.pl -1 forward_reads.fastq -2 reverse_reads.fastq -ref reference_sequence.fasta -out project_name";
if ($reads1 eq "" || $reads2 eq "" || $reference eq "")
{
    die "\n\nYou need to give a reference sequence and forward/reverse read files!\n".$USAGE."\n";
}

my %hash_pairs_ids;
my %hash_pairs_ids2;

my %hashref;
my %hashref2;
undef %hashref;
undef %hashref2;

my %quality1;
my %id1;
my %read1;
my %quality2;
my %id2;
my %read2;

my @reads_tmp = undef;

if ($reads12 eq "")
{
    @reads_tmp = ($reads1, $reads2);
    if ($reads1 eq $reads2)
    {
        die "\nThe two input files are identical, please check the configuration file!\n";
    }
}
else
{
    @reads_tmp = ($reads12);
}

my $check_zip = substr $reads_tmp[0], -2;
my $check_zip2 = substr $reads_tmp[0], -3;
my $firstLine = "";
my $secondLine = "";
my $thirdLine = "";
my $fourthLine = "";
my $fifthLine = "";

if ($check_zip eq "gz")
{
    open (my $FILE, '-|', 'gzip', '-dc', $reads_tmp[0]) or die "Can't open file $reads_tmp[0], $!\n";
    $firstLine = <$FILE>;
    chomp $firstLine;
    $secondLine = <$FILE>;
    chomp $secondLine;
    $thirdLine = <$FILE>;
    chomp $thirdLine;
    $fourthLine = <$FILE>;
    chomp $fourthLine;
    $fifthLine = <$FILE>;
    chomp $fifthLine;
    close $FILE;
}
elsif ($check_zip2 eq "bz2")
{
    open (my $FILE, '-|', 'bzip2', '-dc', $reads_tmp[0]) or die "Can't open file $reads_tmp[0], $!\n";
    $firstLine = <$FILE>;
    chomp $firstLine;
    $secondLine = <$FILE>;
    chomp $secondLine;
    $thirdLine = <$FILE>;
    chomp $thirdLine;
    $fourthLine = <$FILE>;
    chomp $fourthLine;
    $fifthLine = <$FILE>;
    chomp $fifthLine;
    close $FILE;
}
else
{
    open(INPUT, $reads_tmp[0]) or die "No input file found, make sure it are fastq files $!\n";
    $firstLine = <INPUT>;
    chomp $firstLine;
    $secondLine = <INPUT>;
    chomp $secondLine;
    $thirdLine = <INPUT>;
    chomp $thirdLine;
    $fourthLine = <INPUT>;
    chomp $fourthLine;
    $fifthLine = <INPUT>;
    chomp $fifthLine;
    close INPUT;
}

my $firstLine_reverse = "";
if ($reads12 eq "")
{
    if ($check_zip eq "gz")
    {
        open (my $FILE, '-|', 'gzip', '-dc', $reads_tmp[1]) or die "Can't open file $reads_tmp[1], $!\n";
        $firstLine_reverse = <$FILE>;
        chomp $firstLine_reverse;
        close $FILE;
    }
    elsif ($check_zip2 eq "bz2")
    {
        open (my $FILE, '-|', 'bzip2', '-dc', $reads_tmp[1]) or die "Can't open file $reads_tmp[1], $!\n";
        $firstLine_reverse = <$FILE>;
        chomp $firstLine_reverse;
        close $FILE;
    }
    else
    {
        open(INPUT2, $reads_tmp[1]) or die "\n\nNo input file found, make sure it are fastq files $!\n";
        $firstLine_reverse = <INPUT2>;
        chomp $firstLine_reverse;
        close INPUT2;
    }
}
my $no_quality_score = substr $thirdLine, 0, 1;
my $type_of_file = "";
my $code_before_end = substr $firstLine, -2,1;
my $code_before_end0 = substr $firstLine, -1,1;
my $code_before_end_reverse = substr $firstLine_reverse, -2,1;
my $code_before_end0_reverse = substr $firstLine_reverse, -1,1;

my $SRA = "";


if (($code_before_end eq "/" || $code_before_end eq "R" || $code_before_end eq "#") && $firstLine_reverse ne $firstLine && $code_before_end0 eq "1")
{
    $type_of_file = '-1';
}
elsif ($code_before_end eq ":" && $code_before_end0 eq "1" && $firstLine_reverse ne $firstLine && $code_before_end_reverse eq ":" && $code_before_end0_reverse eq "2")
{
    $type_of_file = '-1';
}
elsif($firstLine =~ m/.*(_|\s)(1)(:\w.*\d+:*(\s.*)*\s*\t*)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = "yes";
}
elsif ($code_before_end eq ":" && $code_before_end0 eq "1" && $firstLine_reverse ne $firstLine && $firstLine_reverse eq "")
{
    $type_of_file = '-1';
}
elsif($firstLine =~ m/.*\s(1)(\S*)$/ && $firstLine_reverse ne $firstLine)
{
    my $firstLine_tmp = $firstLine;
    my $test_space = $firstLine_tmp =~ tr/ //;
    $type_of_file = -length($2)-1;
    if ($test_space eq '1')
    {
        $type_of_file = "split";
    }
}
elsif($firstLine =~ m/.*_(1)(:N.*)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = -length($2)-1;
}
elsif($firstLine =~ m/\S*\.(1)(\s(\d+)\s.*)$/ && $firstLine_reverse ne $firstLine)
{
    my $test1 = $3;
    if($fifthLine =~ m/\S*\.(1)(\s(\d+)(\s.*))$/ && $firstLine_reverse ne $firstLine)
    {
         my $test2 = $3;
         if ($test2 eq $test1)
         {
            $type_of_file = -length($2)-1;
         }
         else
         {
            $type_of_file = -length($4);
            $SRA = "yes";
         }
    }
}
elsif($fifthLine =~ m/.*\.(1)(\s(.+)\s.+)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = -length($2)-1;
}
elsif($firstLine_reverse eq $firstLine)
{
    $type_of_file = "identical";
}
elsif($reads12 ne "")
{
    print "\n\nCOMBINED FILE NOT SUPPORTED, PLEASE TRY SEPERATE FILES FOR THE FORWARD AND REVERSE READS!\n\n";
    print OUTPUT4 "\n\nCOMBINED FILE NOT SUPPORTED, PLEASE TRY SEPERATE FILES FOR THE FORWARD AND REVERSE READS!\n\n";
    exit;
}
else
{
    print "\n\nTHE INPUT READS HAVE AN INCORRECT FILE FORMAT!\nPLEASE SEND ME THE ID STRUCTURE!\n\n";
    print OUTPUT4 "\n\nTHE INPUT READS HAVE AN INCORRECT FILE FORMAT!\nPLEASE SEND ME THE ID STRUCTURE!\n\n";
    exit;
}
my $last_character = substr $secondLine, -1;
my $space_at_end = "";
if ($last_character =~ m/\s|\t/g)
{
    #print OUTPUT5 $last_character." LAST2\n";
    $space_at_end = "yes";
}
select(STDERR);
$| = 1;
select(STDOUT); # default
$| = 1;
print "\nScan reference sequence...";

open(INPUT5, $reference) or die "\n\nCan't open reference file $reference, $!\n";
my $ff2 = '0';
my $value_ref2 = "";

while (my $line = <INPUT5>)
{
    if ($ff2 < 1)
    {
        $ff2++;
        next;
    }
    chomp $line;    
    $line =~ tr/actgn/ACTGN/;
    my $first = substr $line, 0, 1;
    my $line3;
    if ($first eq '>' || $first eq '@')
    {
        $line3 = $value_ref2."NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
    }
    else
    {
        $line3 = $value_ref2.$line;
    }
    
    while (length($line3) > (($kmer*3)-1))
    {
        my $value_ref2b = substr $line3, 0, $kmer;
        my $line2 = $line3;
        $line3 = substr $line2, 1;
        my $A = $value_ref2b =~ tr/A/A/;
        my $C = $value_ref2b =~ tr/C/C/;
        my $T = $value_ref2b =~ tr/T/T/;
        my $G = $value_ref2b =~ tr/G/G/;
        my $N = $value_ref2b =~ tr/N/N/;
        if ($A < $kmer-1 && $C < $kmer-1 && $T < $kmer-1 && $G < $kmer-1 && $N < $kmer-1)
        {
            my $value_ref2b_reverse = reverse($value_ref2b);
            $value_ref2b_reverse =~ tr/ACTG/TGAC/;
            
            $hashref{$value_ref2b} .= exists $hashref{$value_ref2b} ? ",$ff2" : $ff2;
            $hashref{$value_ref2b_reverse} .= exists $hashref{$value_ref2b_reverse} ? ",$ff2" : $ff2;
            $hashref2{$ff2} .= exists $hashref2{$ff2} ? "$value_ref2b" : $value_ref2b;
        }
        
        $ff2++;
    }
    $value_ref2 = $line3;
}
while (length($value_ref2) > 1)
{
    my $value_ref2b = substr $value_ref2, 0, $kmer;
    my $value_ref2bc = $value_ref2;
    $value_ref2 = substr $value_ref2bc, 1;
    my $value_ref2b_reverse = reverse($value_ref2b);
    $value_ref2b_reverse =~ tr/ACTG/TGAC/;
    
    $hashref{$value_ref2b} .= exists $hashref{$value_ref2b} ? ",$ff2" : $ff2;
    $hashref{$value_ref2b_reverse} .= exists $hashref{$value_ref2b_reverse} ? ",$ff2" : $ff2;
    $hashref2{$ff2} .= exists $hashref2{$ff2} ? "$value_ref2b" : $value_ref2b;
    $ff2++;
}
close INPUT5;
print "...OK\n";

my $output_file10 = "Filtered_reads_".$project."_R1.fastq";
my $output_file11 = "Filtered_reads_".$project."_R2.fastq";
open(OUTPUT10, ">" .$output_file10) or die "Can't open saved reads1 file $output_file10, $!\n";
open(OUTPUT11, ">" .$output_file11) or die "Can't open saved reads2 file $output_file11, $!\n";

my $type_of_file2 = "";
my $file = "0";
my $count_all_reads = '0';
my $count_filtered_reads = '0';
foreach my $reads_tmp (@reads_tmp)
{
    my $FILE;
    $file++;
    if ($file eq '1')
    {
        select(STDERR);
        $| = 1;
        select(STDOUT); # default
        $| = 1;
        print "\nFilter forward reads...";
    }
    else
    {
        select(STDERR);
        $| = 1;
        select(STDOUT); # default
        $| = 1;
        print "\nFilter reverse reads..."; 
    }
    
    if ($check_zip eq "gz")
    {
        open ($FILE, '-|', 'gzip', '-dc', $reads_tmp) or die "Can't open file $reads_tmp, $!\n";
    }
    elsif ($check_zip2 eq "bz2")
    {
        open ($FILE, '-|', 'bzip2', '-dc', $reads_tmp) or die "Can't open file $reads_tmp, $!\n";
    }
    else
    {
        open($FILE, $reads_tmp) or die "\n\nCan't open file $reads_tmp, $!\n";
    }
    my $line_order = '1';
    my $id = "";
    my $id2 = "";
    my $quality_score = "";
    my $sequence = "";
    my $save_read = "";

    while (my $line = <$FILE>)
    {
        chomp $line;
        if ($line_order eq '1')
        {
            $id = $line;
            #$id2 = substr $line, 0, $type_of_file;
            
            my $code0_SRA = "";
            if ($SRA eq "yes")
            {
                my $code_SRA = substr $id, 0, $type_of_file;
                my @code_SRA = split / /, $code_SRA;
                $code0_SRA = $code_SRA[0];
            }                  
            if ($type_of_file eq "yes")
            {
                if($id =~ m/.*(_|\s)(1|2)(:\w.*\d+:*(\s.*)*\s*\t*)$/)
                {
                    $type_of_file2 = -length($3)-1;
                }   
            }
            if ($type_of_file eq "split")
            {
                my @split = split / /, $id;
                $type_of_file2 = $split[0];
            }
            
            $id2 = substr $id, 0, $type_of_file;
            if ($SRA eq "yes")
            {
                $id2 = $code0_SRA;
            }
            
            if ($type_of_file eq "identical")
            {
               $id2 = $id;
            }
            if ($type_of_file eq "yes")
            {
               $id2 = substr $id, 0, $type_of_file2;
            }
            if ($type_of_file eq "split")
            {
                $id2 = $type_of_file2;
            }
            if ($file eq '2')
            {
                if (exists($hash_pairs_ids{$id2}))
                {
                    $save_read = "yes2";
                }
            }
        }
        elsif ($line_order eq '2')
        {
            my $t = '0';
            $count_all_reads++;
            $sequence = $line;
            
            if ($save_read ne "yes2")
            {
    READ:       while ($t < length($line)-$kmer && $save_read eq "")
                {
                    my $line_part = substr $line, $t, $kmer;
    
                    if (exists($hashref{$line_part}))
                    {                   
                        if ($file eq '2')
                        {
                            $hash_pairs_ids2{$id2} = undef;
                        }
                        else
                        {
                            $hash_pairs_ids{$id2} = undef;
                            
                        }
                        if ($save_read eq "")
                        {
                            $save_read = "yes";
                        }
                        last READ;
                    }
                    $t += 4;
                }
            }    
        }
        elsif ($line_order eq '4')
        {
            if ($save_read ne "")
            {
                $quality_score = $line;
                if ($file eq '1')
                {
                    #$id1{$id2} = $id;
                    #$read1{$id2} = $sequence;
                    #$quality1{$id2} = $quality_score;
                    print OUTPUT10 $id."\n";
                    print OUTPUT10 $sequence."\n";
                    print OUTPUT10 "+\n";
                    print OUTPUT10 $quality_score."\n";
                }
                elsif ($file eq '2' && $save_read eq "yes")
                {
                    $id2{$id2} = $id;
                    $read2{$id2} = $sequence;
                    $quality2{$id2} = $quality_score;
                }
                elsif ($file eq '2' && $save_read eq "yes2")
                {
                    print OUTPUT11 $id."\n";
                    print OUTPUT11 $sequence."\n";
                    print OUTPUT11 "+\n";
                    print OUTPUT11 $quality_score."\n";
                }
                $count_filtered_reads++;
                $save_read = "";
            }
            $line_order = '0';
        }
        $line_order++;
    }
    close $FILE;
    print "...OK\n";
}
%hash_pairs_ids = ();
my $FILE2;
if ($check_zip eq "gz")
{
    open ($FILE2, '-|', 'gzip', '-dc', $reads_tmp[0]) or die "Can't open file $reads_tmp[0], $!\n";
}
elsif ($check_zip2 eq "bz2")
{
    open ($FILE2, '-|', 'bzip2', '-dc', $reads_tmp[0]) or die "Can't open file $reads_tmp[0], $!\n";
}
else
{
    open($FILE2, $reads_tmp[0]) or die "\n\nCan't open file $reads_tmp[0], $!\n";
}
my $line_order = '1';
my $id = "";
my $id2 = "";
my $quality_score = "";
my $sequence = "";
my $save_read = "";
while (my $line = <$FILE2>)
{
    chomp $line;
    if ($line_order eq '1')
    {
        $id = $line;
        #$id2 = substr $line, 0, $type_of_file;
        
        my $code0_SRA = "";
        if ($SRA eq "yes")
        {
            my $code_SRA = substr $id, 0, $type_of_file;
            my @code_SRA = split / /, $code_SRA;
            $code0_SRA = $code_SRA[0];
        }                  
        if ($type_of_file eq "yes")
        {
            if($id =~ m/.*(_|\s)(1|2)(:\w.*\d+:*(\s.*)*\s*\t*)$/)
            {
                $type_of_file2 = -length($3)-1;
            }   
        }
        if ($type_of_file eq "split")
        {
            my @split = split / /, $id;
            $type_of_file2 = $split[0];
        }
        
        my $id2 = substr $id, 0, $type_of_file;
        if ($SRA eq "yes")
        {
            $id2 = $code0_SRA;
        }
        
        if ($type_of_file eq "identical")
        {
           $id2 = $id;
        }
        if ($type_of_file eq "yes")
        {
           $id2 = substr $id, 0, $type_of_file2;
        }
        if ($type_of_file eq "split")
        {
            $id2 = $type_of_file2;
        }

        if (exists($hash_pairs_ids2{$id2}))
        {
            $save_read = "yes";
        }
    }
    elsif ($line_order eq '2')
    {
        my $t = '0';
        $sequence = $line;
    }
    elsif ($line_order eq '4')
    {
        if ($save_read eq "yes")
        {
            $quality_score = $line;
            #print OUTPUT10 $id."\n";
            #print OUTPUT10 $sequence."\n";
            #print OUTPUT10 "+\n";
            #print OUTPUT10 $quality_score."\n";
            $id1{$id2} = $id;
            $read1{$id2} = $sequence;
            $quality1{$id2} = $quality_score;
            $count_filtered_reads++;
            $save_read = "";
        }
        $line_order = '0';
    }
    $line_order++;
}
print "\n------------------------------------\n";
my $filterd_percentage = ($count_filtered_reads/$count_all_reads)*100;
my $filterd_percentage2 = sprintf("%.2f",$filterd_percentage);
print "\nFiltered reads       : ".$count_filtered_reads."\n";
print "Filtered percentage  : ".$filterd_percentage2." %\n";

foreach my $id2 (keys %id1)
{
    if (exists($id2{$id2}))
    {
        print OUTPUT10 $id1{$id2}."\n";
        print OUTPUT10 $read1{$id2}."\n";
        print OUTPUT10 "+\n";
        print OUTPUT10 $quality1{$id2}."\n";
        
        print OUTPUT11 $id2{$id2}."\n";
        print OUTPUT11 $read2{$id2}."\n";
        print OUTPUT11 "+\n";
        print OUTPUT11 $quality2{$id2}."\n";
    }
}
close OUTPUT10;
close OUTPUT11;
