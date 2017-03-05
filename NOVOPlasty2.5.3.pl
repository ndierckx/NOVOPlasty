#!/usr/bin/env perl
######################################################
#         SOFTWARE COPYRIGHT NOTICE AGREEMENT        #
#  Copyright (C) {2015-2017}  {Nicolas Dierckxsens}  #
#              All Rights Reserved                   #
#         See file LICENSE for details.              #
######################################################
#           NOVOPlasty - The Organelle Assembler
#           nicolasdierckxsens@hotmail.com
use Getopt::Long;
use strict;

my $left = '9';
my $right = '15';
my $insert_range = '1.8';
my $insert_range_back = '1.8';
my $noseed = "";  
my $iterations = "1000000";
my $startprint = "1000000";
my $startprint2 = "1000000";
my $high_coverage = "yes";
my $coverage_cut_off = '1000';
my $option = '1';

my $simple_read = "yes";
my $insert_size_correct;
my $overlap;
my $insert_range_b;
my $insert_range_c;
my $genome_range_low;
my $genome_range_high;
my $genome_range ;
my $seed_input0;
my $seed_input;  
my $insert_size;
my $project;
my $read_length;
my $type;
my $encrypt = "no";
my $paired;
my $print_log;
my $bad_read;
my %contigs_id;
my %contigs_end;
my $contig_num = '1';
my $id = "";
my $id_pair = "";
my $id_match = "";
my $id_pair_match = "";
my $contig_end_check;

my %seed;
my %seed_old;
my $counttest1 = '0';
my $indel_split = '0';
my $indel_split_skip;
my %indel_split_skip;
my %indel_split;
my $indel_split_back = '0';
my $indel_split_skip_back;
my %indel_split_skip_back;
my %indel_split_back;
my $insert_range2 = $insert_range;
my $SNR_read;
my $SNR_read_back;
my $SNR_read2;
my $SNR_read_back2;
my %best_extension_prev;
my %best_extension_back_prev;
my %contig_gap_min;
my %contig_gap_max;
my %contig_count;
my $contig_count;
my $progress_before;
my $contig_id2;
my $contig_id1;
my $contig_id3;
my $contig_id4;
my %seed_split;
my $lastbit_contig_prev;
my $split_forward;
my %old_id;
my %old_id2;
my %old_rep;
my %old_rep_old;
my %SNR_length;
my @insert_size;
my %regex;
my %regex_back;
my %regex1;
my %regexb;
my %regex2b;
my %noback;
my %noforward;
my %SNR;
my %SNR_back;
my $SNR_nucleo;
my $SNR_nucleo_back;
my %id_bad;
my $y ='1';
my %position;
my %position_back;
my %bad_read;
my %position_adjust;
my %read_end_tmp;
my %read_end_b_tmp;
my %read_short_end_tmp;
my %read_short_zone_tmp;
my %read_start_tmp;
my %read_start_b_tmp;
my %read_short_start_tmp;
my %read_short_zone_start_tmp;
my %tree;
my %contig_id2;
my %contig_id1;
my %contig_id3;
my %contig_id4;
my %hash;
my %hash2b;
my %hash2c;
my @row;
my %row;
my %match_rep;
my %count_rep;
my %before;
my %before_back;
my %first_before;
my %SNP_active;
my %first_before_back;
my %SNP_active_back;
my %filter_before1;
my %filter_before2;
my %filter_before3;
my %filter_before4;
my %filter_before1_pair;
my %filter_before2_pair;
my %filter_before3_pair;
my %filter_before4_pair;
my %nosecond;
my %repetitive_pair;
my %count_reads;
my %count_reads_all;
my $assembly_length = '1';
my $assembly_success;
my %seeds_check;

my $reads12;
my $reads1;
my $reads2;
my $config;
my $read;
my $deletion;
my $deletion_back;
my $reference_guided;
my $reference_guided_back;
my $contig_read2;
my $contig_read1;
my $first_contig_start;
my $first_contig_start_reverse;
my $finish;
my $repetitive_detect;
my $repetitive_detect_back;
my $repetitive_detect2;
my $repetitive_detect_back2;
my $contig_end;
my $repetitive;
my $before_repetitive;
my $before_repetitive_short;
my $CP_check;
my $before_extension1;
my $before_extension2;
my $before_extension_back1;
my $before_extension_back2;
my $id_split1;
my $id_split2;
my $id_split3;
my $before;
my $before_back;
my $end_repetitive_tmp4;
my $first_before;
my $SNP_active;
my $first_before_back;
my $SNP_active_back;
my $nosecond;
my $first_contig_id;
my $no_contig_id2;
my $no_contig_id1;
my $no_contig_id3;
my $no_contig_id4;
my $rep_detect2;
my $hasL;
my $no_next_seed;
my $count_split;
my @firstSNP_max;
my %yuyu_option;
my %yuyu_option_back;
my $yuyu_option_A;
my $yuyu_option_C;
my $yuyu_option_T;
my $yuyu_option_G;
my $yuyu_option_A_back;
my $yuyu_option_C_back;
my $yuyu_option_T_back;
my $yuyu_option_G_back;
my $extensions_before;
my $before_shorter_skip;
my %before_shorter_skip;
my $before_shorter_skip_back;
my %before_shorter_skip_back;
my $SNR_next_seed;
my $jump_rep;
my %jump_rep;
my $jump_rep_back;
my %jump_rep_back;
my $AT_rich_before;
my $AT_rich_before_back;
my $insert_range_shorter;
my $cp_input;
my $merge_now;

GetOptions (
            "c=s" => \$config,
            ) or die "Incorrect usage!\n";

open(CONFIG, $config) or die "Error: $!\nCan't open the configuration file, please check the manual!\n";

my $ln ='0';

while (my $line = <CONFIG>)
{
    if ($ln eq '0')
    {
        $project = substr $line, 23;
        chomp $project;
    }
    if ($ln eq '1')
    {
        $insert_size = substr $line, 23;
        chomp $insert_size;
    }
    if ($ln eq '2')
    {
        $insert_size_correct = substr $line, 23;
        chomp $insert_size_correct;
    }
    if ($ln eq '3')
    {
        $read_length = substr $line, 23;
        chomp $read_length;
    }
    if ($ln eq '4')
    {
        $type = substr $line, 23;
        chomp $type;
    }
    if ($ln eq '5')
    {
        $genome_range = substr $line, 23;
        chomp $genome_range;
        my @words = split /-/, $genome_range;
        $genome_range_low = $words[0];
        $genome_range_high = $words[1];
    }
    if ($ln eq '6')
    {
        $overlap = substr $line, 23;
        chomp $overlap;
    }
    if ($ln eq '7')
    {
        $insert_range_b = substr $line, 23;
        chomp $insert_range_b;
    }
    if ($ln eq '8')
    {
        $insert_range_c = substr $line, 23;
        chomp $insert_range_c;
    }
    if ($ln eq '9')
    {
        $paired = substr $line, 23;
        chomp $paired;
    }
    if ($ln eq '10')
    {
        $coverage_cut_off = substr $line, 23;
        chomp $coverage_cut_off;
    }
    if ($ln eq '11')
    {
        $print_log = substr $line, 23;
        chomp $print_log;
    }
    if ($ln eq '12')
    {
        $reads12 = substr $line, 23;
        chomp $reads12;
    }
    if ($ln eq '13')
    {
        $reads1 = substr $line, 23;
        chomp $reads1;
    }
    if ($ln eq '14')
    {
        $reads2 = substr $line, 23;
        chomp $reads2;
    }
    if ($ln eq '15')
    {
        $seed_input0 = substr $line, 23;
        chomp $seed_input0;
    }
    if ($ln eq '16')
    {
        $cp_input = substr $line, 23;
        chomp $seed_input0;
    }
    $ln++;
}

close CONFIG;
if ($print_log eq '1' || $print_log eq '2')
{
    $startprint2 = '0';
    $startprint = '1000000';
}

my $USAGE = "\nUsage: perl NOVOPlasty.pl -c config_example.txt";

print "\n\n-----------------------------------------------";
print "\nNOVOPlasty: The Organelle Assembler\n";
print "Version 2.5\n";
print "Author: Nicolas Dierckxsens, (c) 2015-2016\n";
print "-----------------------------------------------\n\n";
print OUTPUT4 "\n\n-----------------------------------------------";
print OUTPUT4 "\nNOVOPlasty: The Organelle Assembler\n";
print OUTPUT4 "Version 2.5\n";
print OUTPUT4 "Author: Nicolas Dierckxsens, (c) 2015-2016\n";
print OUTPUT4 "-----------------------------------------------\n\n";


sub build_partial {
 
my $A = "";
my $G = "";
my $T = "";
my $C = "";

my ($str) = (@_);
    my @re;
    undef @re;
        my $v = length($str);
        my $m = '2';
        my $star = substr $str, 0,1;
        if ($star eq "*")
        {
            substr $str, 0,1, ".";
        }
        
        while ($m < $v) 
        {
            my $str9 = substr $str, $m+1;
            my $str6 = substr $str, 0, $m; 
            my $x = '0';
            my $y = length($str6);
            
            while ($x < $y) 
            {     
                my $str8 = substr $str6, $x+1;
                my $str7 = substr $str6, 0, $x;
                $A = $str7.".".$str8.".".$str9;
                push @re, $A;
                $x++;
            }
            $m++;
        }
    my $re = join ('|' , @re);
    qr/$re/;
}
sub build_partialb {
 
my $A = "";
my $G = "";
my $T = "";
my $C = "";

my ($str) = (@_);
    my @re;
    undef @re;
        my $v = length($str);
        my $m = '2';
    
        while ($m < $v) 
        {
            my $str9 = substr $str, $m+1;
            my $str6 = substr $str, 0, $m; 
            my $x = '0';
            my $y = length($str6);
            
            while ($x < $y) 
            {     
                my $str8 = substr $str6, $x+1;
                my $str7 = substr $str6, 0, $x;
                $A = $str7.".".$str8.".".$str9;
                push @re, $A;
                $x++;
            }
            $m++;
        }
    @re;
}
sub build_partial2b
{
    my $A = "";
    my (%str) = (@_);
    my %re;
    undef %re;

    foreach my $str (keys %str) 
    {
        my $v = length($str)-1;
        my $m = '2';
    
        while ($m < $v) 
        {
            my $str9 = substr $str, $m+1;
            my $str6 = substr $str, 0, $m; 
            my $x = '0';
            my $y = length($str6);
            
            while ($x < $y) 
            {     
                my $str8 = substr $str6, $x+1;
                my $str7 = substr $str6, 0, $x;
                $re{$str7."A".$str8."A".$str9} = undef;
                $re{$str7."A".$str8."C".$str9} = undef;
                $re{$str7."A".$str8."T".$str9} = undef;
                $re{$str7."A".$str8."G".$str9} = undef;
                $re{$str7."C".$str8."A".$str9} = undef;
                $re{$str7."C".$str8."C".$str9} = undef;
                $re{$str7."C".$str8."T".$str9} = undef;
                $re{$str7."C".$str8."G".$str9} = undef;
                $re{$str7."T".$str8."A".$str9} = undef;
                $re{$str7."T".$str8."C".$str9} = undef;
                $re{$str7."T".$str8."T".$str9} = undef;               
                $re{$str7."T".$str8."G".$str9} = undef;
                $re{$str7."G".$str8."A".$str9} = undef;
                $re{$str7."G".$str8."C".$str9} = undef;
                $re{$str7."G".$str8."T".$str9} = undef;
                $re{$str7."G".$str8."G".$str9} = undef;
                $x++;
            }
            $m++;
        }
    }
    %re;
}
sub build_partial3b {
 
my $A = "";
my @str = @_;
my $str = $str[0];
my $reverse = $str[1];
my %re;
my %re2;
my %re3;
undef %re;
undef %re2;
undef %re3;
my $hasdot = $str =~ tr/\.//;
my $hasstar = $str =~ tr/\*//;

    
    $re2{$str} = undef;

    if ($hasstar > 0)
    {
        foreach my $str (keys %re2)
        {
            if ($reverse eq "reverse")
            {
                my $str2 = $str;
                $str2 =~ s/\*.//g;
                $re3{$str2} = undef;
            }
            else
            {
                my $str2 = $str;
                $str2 =~ s/.\*//g;
                $re3{$str2} = undef;
            }
            if ($reverse eq "reverse")
            {
                $str =~ tr/\*//d;
                my $temp_sre = $str;
                $str = substr $temp_sre, 0, -$hasstar;
                $re3{$str} = undef;
            }
            else
            {
                $str =~ tr/\*//d;
                my $temp_sre = $str;
                $str = substr $temp_sre, $hasstar;
                $re3{$str} = undef;
            }
        }  
    }
    else
    {
        %re3 = %re2;
    }

if ($hasdot > 0)
{
    foreach my $str (keys %re3) 
    {
        $hasdot = $str =~ tr/\.//;
        if ($hasdot eq '1')
        {
            if ($str =~ m/^(\w*)\.(\w*)$/)
            {
                my $str7 = $1;
                my $str8 = $2;
                
                $re{$str7."A".$str8} = undef;
                $re{$str7."C".$str8} = undef;
                $re{$str7."T".$str8} = undef;
                $re{$str7."G".$str8} = undef;
            }
        }
        elsif ($hasdot eq '2')
        {
            if ($str =~ m/^(\w*)\.(\w*)\.(\w*)$/)
            {
                my $str7 = $1;
                my $str8 = $2;
                my $str9 = $3;
                
                    $A = $str7."A".$str8."A".$str9; #two substitutions
                    $re{$A} = undef;
                    $A = $str7."A".$str8."C".$str9;  
                    $re{$A} = undef;
                    $A = $str7."A".$str8."T".$str9;  
                    $re{$A} = undef;
                    $A = $str7."A".$str8."G".$str9;  
                    $re{$A} = undef;
                    $A = $str7."C".$str8."A".$str9; #two substitutions
                    $re{$A} = undef;
                    $A = $str7."C".$str8."C".$str9;  
                    $re{$A} = undef;
                    $A = $str7."C".$str8."T".$str9;  
                    $re{$A} = undef;
                    $A = $str7."C".$str8."G".$str9;  
                    $re{$A} = undef;
                    $A = $str7."T".$str8."A".$str9; #two substitutions
                    $re{$A} = undef;
                    $A = $str7."T".$str8."C".$str9;  
                    $re{$A} = undef;
                    $A = $str7."T".$str8."T".$str9;  
                    $re{$A} = undef;
                    $A = $str7."T".$str8."G".$str9;  
                    $re{$A} = undef;
                    $A = $str7."G".$str8."A".$str9; #two substitutions
                    $re{$A} = undef;
                    $A = $str7."G".$str8."C".$str9;  
                    $re{$A} = undef;
                    $A = $str7."G".$str8."T".$str9;  
                    $re{$A} = undef;
                    $A = $str7."G".$str8."G".$str9;  
                    $re{$A} = undef;
            }
        }
        elsif ($hasdot eq '3')
        {
            my $str1 = $str;
            my $x = '0';

                my @combi = ('A','C','T','G');
                while ($x < 4)
                {
                    $str1 = $str;
                    $str1 =~ s/\./$combi[$x]/;
                    my $x2 = '0';
                    while ($x2 < 4)
                    {
                        my $str2 = $str1;
                        $str2 =~ s/\./$combi[$x2]/;
                        my $x3 = '0';
                        while ($x3 < 4)
                        {
                            my $str3 = $str2;
                            $str3 =~ s/\./$combi[$x3]/;
                            $re{$str3} = undef;
                            $x3++; 
                        }
                        $x2++; 
                    }
                    $x++; 
                }
        }
        elsif ($hasdot eq '4')
        {
            my $str1 = $str;
            my $x = '0';
                my @combi = ('A','C','T','G');
                while ($x < 4)
                {
                    $str1 = $str;
                    $str1 =~ s/\./$combi[$x]/;
                    my $x2 = '0';
                    while ($x2 < 4)
                    {
                        my $str2 = $str1;
                        $str2 =~ s/\./$combi[$x2]/;
                        my $x3 = '0';
                        while ($x3 < 4)
                        {
                            my $str3 = $str2;
                            $str3 =~ s/\./$combi[$x3]/;
                            my $x4 = '0';
                            while ($x4 < 4)
                            {
                               my $str4 = $str3;
                               $str4 =~ s/\./$combi[$x4]/;
                               $re{$str4} = undef;
                               $x4++; 
                            }
                            $x3++; 
                        }
                        $x2++; 
                    }
                    $x++; 
                }
        }
        elsif ($hasdot eq '5')
        {
            my $str1 = $str;
            my $x = '0';
                my @combi = ('A','C','T','G');
                while ($x < 4)
                {
                    $str1 = $str;
                    $str1 =~ s/\./$combi[$x]/;
                    my $x2 = '0';
                    while ($x2 < 4)
                    {
                        my $str2 = $str1;
                        $str2 =~ s/\./$combi[$x2]/;
                        my $x3 = '0';
                        while ($x3 < 4)
                        {
                            my $str3 = $str2;
                            $str3 =~ s/\./$combi[$x3]/;
                            my $x4 = '0';
                            while ($x4 < 4)
                            {
                                my $str4 = $str3;
                                $str4 =~ s/\./$combi[$x4]/;
                                my $x5 = '0';
                                while ($x5 < 4)
                                {
                                    my $str5 = $str4;
                                    $str5 =~ s/\./$combi[$x5]/;
                                    $re{$str5} = undef;
                                    $x5++; 
                                }
                               $x4++; 
                            }
                            $x3++; 
                        }
                        $x2++; 
                    }
                    $x++; 
                }
        }
        elsif ($hasdot eq '6')
        {
            my $str1 = $str;
            my $x = '0';
                my @combi = ('A','C','T','G');
                while ($x < 4)
                {
                    $str1 = $str;
                    $str1 =~ s/\./$combi[$x]/;
                    my $x2 = '0';
                    while ($x2 < 4)
                    {
                        my $str2 = $str1;
                        $str2 =~ s/\./$combi[$x2]/;
                        my $x3 = '0';
                        while ($x3 < 4)
                        {
                            my $str3 = $str2;
                            $str3 =~ s/\./$combi[$x3]/;
                            my $x4 = '0';
                            while ($x4 < 4)
                            {
                                my $str4 = $str3;
                                $str4 =~ s/\./$combi[$x4]/;
                                my $x5 = '0';
                                while ($x5 < 4)
                                {
                                    my $str5 = $str4;
                                    $str5 =~ s/\./$combi[$x5]/;
                                    my $x6 = '0';
                                    while ($x6 < 4)
                                    {
                                        my $str6 = $str5;
                                        $str6 =~ s/\./$combi[$x6]/;
                                        $re{$str6} = undef;
                                        $x6++;
                                    }
                                    $x5++;
                                }
                               $x4++; 
                            }
                            $x3++; 
                        }
                        $x2++; 
                    }
                    $x++; 
                }
                
        }
        elsif ($hasdot eq '7')
        {
            my $str1 = $str;
            my $x = '0';
                my @combi = ('A','C','T','G');
                while ($x < 4)
                {
                    $str1 = $str;
                    $str1 =~ s/\./$combi[$x]/;
                    my $x2 = '0';
                    while ($x2 < 4)
                    {
                        my $str2 = $str1;
                        $str2 =~ s/\./$combi[$x2]/;
                        my $x3 = '0';
                        while ($x3 < 4)
                        {
                            my $str3 = $str2;
                            $str3 =~ s/\./$combi[$x3]/;
                            my $x4 = '0';
                            while ($x4 < 4)
                            {
                                my $str4 = $str3;
                                $str4 =~ s/\./$combi[$x4]/;
                                my $x5 = '0';
                                while ($x5 < 4)
                                {
                                    my $str5 = $str4;
                                    $str5 =~ s/\./$combi[$x5]/;
                                    my $x6 = '0';
                                    while ($x6 < 4)
                                    {
                                        my $str6 = $str5;
                                        $str6 =~ s/\./$combi[$x6]/;
                                        my $x7 = '0';
                                        while ($x7 < 4)
                                        {
                                            my $str7 = $str6;
                                            $str7 =~ s/\./$combi[$x7]/;
                                            $re{$str7} = undef;
                                            $x7++;
                                        }
                                        $x6++;
                                    }
                                    $x5++;
                                }
                               $x4++; 
                            }
                            $x3++; 
                        }
                        $x2++; 
                    }
                    $x++; 
                }
        }
        elsif ($hasdot eq '8')
        {
            my $str1 = $str;
            my $x = '0';
                my @combi = ('A','C','T','G');
                while ($x < 4)
                {
                    $str1 = $str;
                    $str1 =~ s/\./$combi[$x]/;
                    my $x2 = '0';
                    while ($x2 < 4)
                    {
                        my $str2 = $str1;
                        $str2 =~ s/\./$combi[$x2]/;
                        my $x3 = '0';
                        while ($x3 < 4)
                        {
                            my $str3 = $str2;
                            $str3 =~ s/\./$combi[$x3]/;
                            my $x4 = '0';
                            while ($x4 < 4)
                            {
                                my $str4 = $str3;
                                $str4 =~ s/\./$combi[$x4]/;
                                my $x5 = '0';
                                while ($x5 < 4)
                                {
                                    my $str5 = $str4;
                                    $str5 =~ s/\./$combi[$x5]/;
                                    my $x6 = '0';
                                    while ($x6 < 4)
                                    {
                                        my $str6 = $str5;
                                        $str6 =~ s/\./$combi[$x6]/;
                                        my $x7 = '0';
                                        while ($x7 < 4)
                                        {
                                            my $str7 = $str6;
                                            $str7 =~ s/\./$combi[$x7]/;
                                            my $x8 = '0';
                                            while ($x8 < 4)
                                            {
                                                my $str8 = $str7;
                                                $str8 =~ s/\./$combi[$x8]/;
                                                $re{$str8} = undef;
                                                $x8++;
                                            }
                                            $x7++;
                                        }
                                        $x6++;
                                    }
                                    $x5++;
                                }
                               $x4++; 
                            }
                            $x3++; 
                        }
                        $x2++; 
                    }
                    $x++; 
                }
        }
        else
        {
            $re{$str} = undef;
        }
    }
}
else
{
    %re = %re3;
}
    %re;
}
sub build_partial3c {
 
my $A = "";
my @str = @_;
my $str = $str[0];
my $reverse = $str[1];
my %re2;
my %re;
undef %re2;
undef %re;
my $hasstar = $str =~ tr/\*//;

    $re2{$str} = undef;

    if ($hasstar > 0)
    {
        foreach my $str (keys %re2)
        {
            if ($reverse eq "reverse")
            {
                my $str2 = $str;
                $str2 =~ s/\*.//g;
                $re{$str2} = undef;
            }
            else
            {
                my $str2 = $str;
                $str2 =~ s/.\*//g;
                $re{$str2} = undef;
            }
            if ($reverse eq "reverse")
            {
                $str =~ tr/\*//d;
                my $temp_sre = $str;
                $str = substr $temp_sre, 0, -$hasstar;
                $re{$str} = undef;
            }
            else
            {
                $str =~ tr/\*//d;
                my $temp_sre = $str;
                $str = substr $temp_sre, $hasstar;
                $re{$str} = undef;
            }
        } 
    }
    else
    {
        %re = %re2;
    }
    %re;
}
sub uniq
{
    my %seen;
    grep !$seen{$_}++, @_;
}
sub encrypt
{
    my @value = @_;
    my $value = $value[0];
    
    $value =~ s/AA/0/g;
    $value =~ s/CC/1/g;
    $value =~ s/TT/2/g;
    $value =~ s/GG/3/g;
    $value =~ s/AC/4/g;
    $value =~ s/AG/5/g;
    $value =~ s/AT/6/g;
    $value =~ s/CT/7/g;
    $value =~ s/CA/8/g;
    $value =~ s/CG/9/g;
    $value =~ s/TA/Z/g;
    $value =~ s/TC/U/g;
    $value =~ s/TG/I/g;
    $value =~ s/GA/O/g;
    $value =~ s/GC/P/g;
    $value =~ s/GT/Q/g;
    return $value;
}
sub decrypt
{
    my @value = @_;
    my $value = $value[0];
    $value =~ s/Q/GT/g;
    $value =~ s/P/GC/g;
    $value =~ s/O/GA/g;
    $value =~ s/I/TG/g;
    $value =~ s/U/TC/g;
    $value =~ s/Z/TA/g;
    $value =~ s/9/CG/g;
    $value =~ s/8/CA/g;
    $value =~ s/7/CT/g;
    $value =~ s/6/AT/g;
    $value =~ s/5/AG/g;
    $value =~ s/4/AC/g;
    $value =~ s/3/GG/g;
    $value =~ s/2/TT/g;
    $value =~ s/1/CC/g;
    $value =~ s/0/AA/g;
    return $value;
}
sub IUPAC
{
    my @snps = @_;
    my $A = $snps[0];
    my $C = $snps[1];
    my $T = $snps[2];
    my $G = $snps[3];
    my $iupac;
    
    if ($A + $C > 0.93*($A+$C+$T+$G))
    {
        $iupac = "M";
    }
    elsif ($A + $G > 0.93*($A+$C+$T+$G))
    {
        $iupac = "R"; 
    }
    elsif ($A + $T > 0.93*($A+$C+$T+$G))
    {
        $iupac = "W"; 
    }
    elsif ($C + $G > 0.93*($A+$C+$T+$G))
    {
        $iupac = "S"; 
    }
    elsif ($C + $T > 0.93*($A+$C+$T+$G))
    {
        $iupac = "Y"; 
    }
    elsif ($G + $T > 0.93*($A+$C+$T+$G))
    {
        $iupac = "K"; 
    }
    elsif ($A + $C + $G eq ($A+$C+$T+$G))
    {
        $iupac = "V"; 
    }
    elsif ($A + $C + $T eq ($A+$C+$T+$G))
    {
        $iupac = "H"; 
    }
    elsif ($A + $T + $G eq ($A+$C+$T+$G))
    {
        $iupac = "D"; 
    }
    elsif ($T + $C + $G eq ($A+$C+$T+$G))
    {
        $iupac = "B"; 
    }
    else
    {
        $iupac = "."; 
    }
    return $iupac;
}
sub correct
{
        my @str = @_;
        my $read_correct = $str[0];
        my $cc = '0';
        
        $cc = '0';
        my %read_part;
        my %read_part_reverse;
        my %read_matches;
        my @read_matches;
        undef %read_part;
        undef %read_part_reverse;
        undef %read_matches;
        undef @read_matches;
        if ($y > $startprint2)
        {
            print OUTPUT5 $read_correct."\n";
        }
        while ($cc < (length($read_correct)/2)-$overlap+($overlap/2))
        {
            my $read_part = substr $read_correct, $cc, $overlap;
            $read_part =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
            %read_part = build_partial3b $read_part;
            foreach my $read_part2 (keys %read_part)
            {
                if (exists($hash2b{$read_part2}))
                {
                    my $id_part = $hash2b{$read_part2};
                    $id_part = substr $id_part, 1;
                    my @id_part = split /,/,$id_part;
                    
                    foreach my $id_part2 (@id_part)
                    {                    
                        my $id_part2_end = substr $id_part2, -1 , 1, "";
                        if (exists($hash{$id_part2}))
                        {
                            my @id_part2_tmp = split /,/,$hash{$id_part2};
                            my $id_part2_tmp_new;         
                            if (exists($hash{$id_part2}))
                            {
                                my @id_part2_tmp = split /,/,$hash{$id_part2};
                                my $id_part2_tmp_new;

                                if ($id_part2_end eq "1")
                                {
                                    $id_part2_tmp_new = $id_part2_tmp[0];
                                }
                                elsif ($id_part2_end eq "2")
                                {
                                    $id_part2_tmp_new = $id_part2_tmp[1];
                                }
                                if ($encrypt eq "yes")
                                {
                                    $id_part2_tmp_new = decrypt $id_part2_tmp_new;
                                }
                                my $id_part2_tmp_new3 = substr $id_part2_tmp_new,0, length($id_part2_tmp_new)- $cc+$left;                          
                                while (length($id_part2_tmp_new3) < length($id_part2_tmp_new))
                                {
                                    $id_part2_tmp_new3 = "N".$id_part2_tmp_new3;
                                }
                                if ($cc < $left)
                                {
                                    $id_part2_tmp_new3 = substr $id_part2_tmp_new,$left-$cc, length($id_part2_tmp_new)- $cc+$left;
                                }                        
                                $read_matches{$id_part2_tmp_new3} = $cc;
                                push @read_matches, $id_part2_tmp_new3;
                            }
                        }
                    }
                }
            }
            my $read_part_d = $read_part;
                                    
            $read_part_d =~ tr/ATCG/TAGC/;
            my $read_part_reverse = reverse($read_part_d);
            %read_part_reverse = build_partial3b $read_part_reverse;
            
            foreach my $read_part2 (keys %read_part_reverse)
            {
                if (exists($hash2c{$read_part2}))
                {
                    my $id_part = $hash2c{$read_part2};
                    $id_part = substr $id_part, 1;
                    my @id_part = split /,/,$id_part;
                    
                    foreach my $id_part2 (@id_part)
                    {
                        my $id_part2_end = substr $id_part2, -1 , 1, "";
                        if (exists($hash{$id_part2}))
                        {
                            my @id_part2_tmp = split /,/,$hash{$id_part2};
                            my $id_part2_tmp_new;

                            if ($id_part2_end eq "1")
                            {
                                $id_part2_tmp_new = $id_part2_tmp[0];
                            }
                            elsif ($id_part2_end eq "2")
                            {
                                $id_part2_tmp_new = $id_part2_tmp[1];
                            }
                            if ($encrypt eq "yes")
                            {
                                $id_part2_tmp_new = decrypt $id_part2_tmp_new;
                            } 
                            my $id_part2_tmp_new2b = $id_part2_tmp_new;
                            $id_part2_tmp_new2b =~ tr/ATCG/TAGC/;
                            $id_part2_tmp_new = reverse($id_part2_tmp_new2b);
                            my $id_part2_tmp_new3 = substr $id_part2_tmp_new,0, length($id_part2_tmp_new)- $cc+$right;
                            while (length($id_part2_tmp_new3) < length($id_part2_tmp_new))
                            {
                                $id_part2_tmp_new3 = "N".$id_part2_tmp_new3;
                            }
                            if ($cc < $right)
                            {
                                $id_part2_tmp_new3 = substr $id_part2_tmp_new,$right-$cc, length($id_part2_tmp_new)- $cc+$right;
                            }                        
                            $read_matches{$id_part2_tmp_new3} = $cc;
                            push @read_matches, $id_part2_tmp_new3;
                        }
                    }
                }
            }                                                                
            $cc++;
        }
        $cc -= ($overlap/2);
        while ($cc < length($read_correct)-$overlap)
        {
            my $read_part = substr $read_correct, $cc, $overlap;
            $read_part =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
            %read_part = build_partial3b $read_part;
            foreach my $read_part2 (keys %read_part)
            {
                if (exists($hash2c{$read_part2}))
                {
                    my $id_part = $hash2c{$read_part2};
                    $id_part = substr $id_part, 1;
                    my @id_part = split /,/,$id_part;
                    
                    foreach my $id_part2 (@id_part)
                    {
                        my $id_part2_end = substr $id_part2, -1 , 1, "";
                        if (exists($hash{$id_part2}))
                        {
                            my @id_part2_tmp = split /,/,$hash{$id_part2};
                            my $id_part2_tmp_new;

                            if ($id_part2_end eq "1")
                            {
                                $id_part2_tmp_new = $id_part2_tmp[0];
                            }
                            elsif ($id_part2_end eq "2")
                            {
                                $id_part2_tmp_new = $id_part2_tmp[1];
                            }
                            if ($encrypt eq "yes")
                            {
                                $id_part2_tmp_new = decrypt $id_part2_tmp_new;
                            } 
                            my $id_part2_tmp_new3 = substr $id_part2_tmp_new, -$cc -$overlap-$right;
                            if (length($id_part2_tmp_new)-$cc-$overlap < $right)
                            {
                                $id_part2_tmp_new3 = substr $id_part2_tmp_new, -$cc -$overlap-$right, -($right-length($id_part2_tmp_new)+$cc+$overlap);
                                while (length($id_part2_tmp_new3) < length($id_part2_tmp_new))
                                {
                                    $id_part2_tmp_new3 = "N".$id_part2_tmp_new3;
                                }
                            }
                            $read_matches{$id_part2_tmp_new3} = $cc;
                            push @read_matches, $id_part2_tmp_new3;
                        }
                    }
                }
            }
            my $read_part_d = $read_part;
                                    
            $read_part_d =~ tr/ATCG/TAGC/;
            my $read_part_reverse = reverse($read_part_d);
            %read_part_reverse = build_partial3b $read_part_reverse;
            
            foreach my $read_part2 (keys %read_part_reverse)
            {
                if (exists($hash2b{$read_part2}))
                {
                    my $id_part = $hash2b{$read_part2};
                    $id_part = substr $id_part, 1;
                    my @id_part = split /,/,$id_part;
                    foreach my $id_part2 (@id_part)
                    {
                        my $id_part2_end = substr $id_part2, -1 , 1, "";
                        if (exists($hash{$id_part2}))
                        {
                            my @id_part2_tmp = split /,/,$hash{$id_part2};
                            my $id_part2_tmp_new;
                            
                            if ($id_part2_end eq "1")
                            {
                                $id_part2_tmp_new = $id_part2_tmp[0];
                            }
                            elsif ($id_part2_end eq "2")
                            {
                                $id_part2_tmp_new = $id_part2_tmp[1];
                            }
                            if ($encrypt eq "yes")
                            {
                                $id_part2_tmp_new = decrypt $id_part2_tmp_new;
                            }
                            my $id_part2_tmp_new2b = $id_part2_tmp_new;
                            $id_part2_tmp_new2b =~ tr/ATCG/TAGC/;
                            $id_part2_tmp_new = reverse($id_part2_tmp_new2b);
                            
                            my $id_part2_tmp_new3 = substr $id_part2_tmp_new, -$cc -$overlap-$left;
                            
                            if (length($id_part2_tmp_new)-$cc-$overlap < $left)
                            {
                                $id_part2_tmp_new3 = substr $id_part2_tmp_new, -$cc -$overlap-$left, -($left-length($id_part2_tmp_new)+$cc+$overlap);
                                while (length($id_part2_tmp_new3) < length($id_part2_tmp_new))
                                {
                                    $id_part2_tmp_new3 = "N".$id_part2_tmp_new3;
                                }
                            }
                            $read_matches{$id_part2_tmp_new3} = $cc;
                            push @read_matches, $id_part2_tmp_new3;
                        }
                    }
                }
            }                  
            $cc++;
        }
        my $count_cor = '0';
        $count_cor++ for (keys %read_matches);
        if ($count_cor < 4)
        {
            if ($y > $startprint2)
            {
                print OUTPUT5 "\nBAD_READ\n";
            }
            if ($y eq '1')
            {
                $bad_read = "yes";
                goto FIRST_SEED;
            }
            else
            {
                goto CORRECT_END;
            }
        }
        my $l = '0';
        my $corrected_read = "";
        my @charso = split//, $read_correct;

        while ($l < length($read_correct))
        {
            my $A = '0';
            my $C = '0';
            my $T = '0';
            my $G = '0';
            
POINT:      foreach my $extensions (@read_matches)
            {                                
                my @chars = split //, $extensions;
                
                if ($chars[$l] eq "A")
                {
                    $A++;
                }
                elsif ($chars[$l] eq "C")
                {
                    $C++;
                }
                elsif ($chars[$l] eq "T")
                {
                    $T++;
                }
                elsif ($chars[$l] eq "G")
                {
                    $G++;
                }
                if (($A + $C + $T + $G) > 500)
                {
                    last POINT;
                }
            }
            if ($A >= ($C + $T + $G)*1.8)
            {
                $corrected_read = $corrected_read."A";
            }
            elsif ($C >= ($A + $T + $G)*1.8)
            {
                $corrected_read = $corrected_read."C";
            }
            elsif ($T >= ($A + $C + $G)*1.8)
            {
                $corrected_read = $corrected_read."T";
            }
            elsif ($G >= ($C + $T + $A)*1.8)
            {
                $corrected_read = $corrected_read."G";
            }
            else
            {
                $corrected_read = $corrected_read.$charso[$l];
            }
            $l++;
        }
        if ($y > $startprint2)
        {
            print OUTPUT5 $corrected_read." CORRECTED READ\n";
        }
        my $dot_bad = $corrected_read =~ tr/\./\./;
        if ($dot_bad > 3)
        {
            if ($y > $startprint2)
            {
                print OUTPUT5 "\nBAD_READ2\n";
                delete $seed{$id};
            }
            if ($y eq '1')
            {
                $bad_read = "yes";
                goto FIRST_SEED;
            }
        }
            
        $read_correct = $corrected_read;
        if ($y eq '1')
        {
            $read = $read_correct;
        }

        $seed{$id} = $read_correct;
        $seeds_check{$id} = undef;
        
        delete $seed_old{$id_pair};
        return $read_correct;
CORRECT_END:
}


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
my $output_file4  = "log_".$project.".txt";
my $output_file5  = "log_extended_".$project.".txt";
my $output_file6  = "contigs_tmp_".$project.".txt";
my $output_file7  = "Merged_contigs_".$project.".txt";

open(INPUT, $reads_tmp[0]) or die "\nNo input file found, make sure it are fastq files $!\n";
open(INPUT3, $seed_input0)  or die "\nCan't open the seed file, $!\n";
open(OUTPUT4, ">" .$output_file4) or die "\nCan't open file $output_file4, $!\n";
open(OUTPUT6, ">" .$output_file6) or die "\nCan't open file $output_file6, $!\n";

if ($print_log eq '1' || $print_log eq '2')
{
    open(OUTPUT5, ">" .$output_file5) or die "\nCan't open file $output_file5, $!\n";
}

open(OUTPUT7, ">" .$output_file7) or die "\nCan't open file $output_file7, $!\n";

select(STDERR);
$| = 1;
select(STDOUT); # default
$| = 1;
print "\nReading Input...";

my $firstLine = <INPUT>;
chomp $firstLine;
my $secondLine = <INPUT>;
chomp $secondLine;
my $thirdLine = <INPUT>;
close INPUT;
my $firstLine_reverse = "";
if ($reads12 eq "")
{
    open(INPUT2, $reads_tmp[1]) or die "\n\nNo input file found, make sure it are fastq files $!\n";
    $firstLine_reverse = <INPUT2>;
    chomp $firstLine_reverse;
}
close INPUT2;
my $no_quality_score = substr $thirdLine, 0, 1;
my $type_of_file;
my $code_before_end = substr $firstLine, -2,1;
if ($code_before_end eq "/" && $firstLine_reverse ne $firstLine)
{
    $type_of_file = '-1';
}
elsif($firstLine =~ m/.*(_|\s)(1)(:N.*\d+\s*\t*)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = "yes";
}
elsif($firstLine =~ m/.*\s(1)(.*)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = -length($2)-1;
}
elsif($firstLine =~ m/.*_(1)(:N.*)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = -length($2)-1;
}
elsif($firstLine =~ m/.*\.(1)(\s(\d+)\s.*)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = -length($2)-1;
}
elsif($firstLine =~ m/.*\.(1)(\s(.+)\s.+)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = -length($2)-1;
}
elsif($firstLine_reverse eq $firstLine)
{
    $type_of_file = "identical";
}
else
{
    print "\n\nTHE INPUT READS HAVE AN INCORRECT FILE FORMAT!\nPLEASE SEND ME THE ID STRUCTURE!\n\n";
    print OUTPUT4 "\n\nTHE INPUT READS HAVE AN INCORRECT FILE FORMAT!\nPLEASE SEND ME THE ID STRUCTURE!\n\n";
    exit;
}

my $last_character = substr $secondLine, -1;
my $space_at_end;
if ($last_character =~ m/\s|\t/g)
{
    print OUTPUT5 $last_character." LAST2\n";
    $space_at_end = "yes";
}
print "...OK\n";
print "\nInput parameters from the configuration file:   *** Verify if everything is correct ***\n\n";
print "Project name         = ".$project."\n";
print "Insert size          = ".$insert_size."\n";
print "Insert size auto     = ".$insert_size_correct."\n";
print "Read Length          = ".$read_length."\n";
print "Type                 = ".$type."\n";
print "Genome range         = ".$genome_range."\n";
print "K-mer                = ".$overlap."\n";
print "Insert range         = ".$insert_range_b."\n";
print "Insert range strict  = ".$insert_range_c."\n";
print "Paired/Single        = ".$paired."\n";
print "Coverage Cut off     = ".$coverage_cut_off."\n";
print "Extended log         = ".$print_log."\n";
print "Combined reads       = ".$reads12."\n";
print "Forward reads        = ".$reads1."\n";
print "Reverse reads        = ".$reads2."\n";
print "Seed Input           = ".$seed_input0."\n";
print "Chloroplast sequence = ".$cp_input."\n\n";

print OUTPUT4 "\n\n-----------------------------------------------";
print OUTPUT4 "\nNOVOPlasty: The Organelle Assembler\n";
print OUTPUT4 "Version 2.5\n";
print OUTPUT4 "Author: Nicolas Dierckxsens, (c) 2015-2016\n";
print OUTPUT4 "-----------------------------------------------\n\n";

print OUTPUT4 "\nInput parameters from the configuration file:   *** Verify if everything is correct ***\n\n";
print OUTPUT4 "Project name         = ".$project."\n";
print OUTPUT4 "Insert size          = ".$insert_size."\n";
print OUTPUT4 "Insert size auto     = ".$insert_size_correct."\n";
print OUTPUT4 "Read Length          = ".$read_length."\n";
print OUTPUT4 "Type                 = ".$type."\n";
print OUTPUT4 "Genome range         = ".$genome_range."\n";
print OUTPUT4 "K-mer                = ".$overlap."\n";
print OUTPUT4 "Insert range         = ".$insert_range_b."\n";
print OUTPUT4 "Insert range strict  = ".$insert_range_c."\n";
print OUTPUT4 "Paired/Single        = ".$paired."\n";
print OUTPUT4 "Coverage Cut off     = ".$coverage_cut_off."\n";
print OUTPUT4 "Extended log         = ".$print_log."\n";
print OUTPUT4 "Combined reads       = ".$reads12."\n";
print OUTPUT4 "Forward reads        = ".$reads1."\n";
print OUTPUT4 "Reverse reads        = ".$reads2."\n";
print OUTPUT4 "Seed Input           = ".$seed_input0."\n";
print OUTPUT4 "Chloroplast sequence = ".$cp_input."\n\n";


my %contigs;
my %hashref;
my %hashref2;
my $ff = '0';
my $value_ref = "";

if ($cp_input ne "")
{
    select(STDERR);
    $| = 1;
    select(STDOUT); # default
    $| = 1;
    print "\nScan chloroplast sequence...";
    open(INPUT2, $cp_input) or die "\n\nCan't open chloroplast file $cp_input, $!\n";
    
    while (my $line = <INPUT2>)
    {
        if ($ff < 1)
        {
            $ff++;
            next;
        }
        chomp $line;    
        $line =~ tr/actgn/ACTGN/;
        
        my $line3 = $value_ref.$line;
        
        while (length($line3) > (($overlap*3)-1))
        {
            my $value_ref2 = substr $line3, 0, $overlap;
            my $line2 = $line3;
            $line3 = substr $line2, 1;
            
            $hashref{$value_ref2} .= exists $hashref{$value_ref2} ? ",$ff" : $ff;
            $hashref2{$ff} .= exists $hashref2{$ff} ? "$value_ref2" : $value_ref2;
            $ff++;
        }
        $value_ref = $line3;
    }
    close INPUT2;
    print "...OK\n";
}

select(STDERR);
$| = 1;
select(STDOUT); # default
$| = 1;
print "\nBuilding Hash Table...";

my $file1;

if ($reads12 eq "")
{
    $file1 = $reads1;
}
else
{
    $file1 = $reads12;
}
    
my $file_count= '0';
my $type_of_file2 = "";

foreach my $reads_tmp (@reads_tmp)
{
    open(INPUT, $reads_tmp) or die "\n\nCan't open file $reads_tmp, $!\n";

    my $N = '0';
    my $f = "";
    my $code = "";
    my $code_new = '0';
    my $value = "";
    
    while (my $line = <INPUT>)
    {
        chomp $line;
       
        if ($f eq "yes" && $no_quality_score ne "@" && $no_quality_score ne ">")
        {
            $f = "yes2";
            next;
        }
        if ($f eq "yes" && ($no_quality_score eq "@" || $no_quality_score eq ">"))
        {
            $code = $line;
            $f = "no";
        }
        elsif ($f eq "yes2")
        {
            $f = "";
            next;
        }
        elsif ($f eq "no")
        {
            $value = $line;
            my $containN = $line =~ tr/N/N/;
            
            if ($containN < 5 && (length($line) > $read_length/1.5 || length($line) > $overlap+$right+$left+10))
            {
                if ($space_at_end eq "yes")
                {
                    chop($value);
                }
                my $code2 = $code;
                my $code_end = substr $code2, $type_of_file,1;
                            
                if ($type_of_file eq "identical")
                {
                    if ($file_count eq '0')
                    {
                        $code_end = "1";
                    }
                    elsif ($file_count eq '1')
                    {
                        $code_end = "2";
                    }    
                }
                if ($type_of_file eq "yes")
                {
                    if($code2 =~ m/.*(_|\s)(1|2)(:N.*\d+\s*\t*)$/)
                    {
                        $code_end = $2;
                        $type_of_file2 = -length($3)-1;
                    }   
                }             
                if ($code_end eq "1" )
                {
                    my $code0 = substr $code2, 0, $type_of_file;
                    
                    if ($type_of_file eq "identical")
                    {
                       $code0 = $code2;
                    }
                    if ($type_of_file eq "yes")
                    {
                       $code0 = substr $code2, 0, $type_of_file2;
                    }
                  
                    $hash{$code0."1"} = $value;         
                }
                if ($code_end eq "2" )
                {
                    my $code1 = substr $code2, 0, $type_of_file;
                    
                    if ($type_of_file eq "identical")
                    {
                       $code1 = $code2;
                    }
                    if ($type_of_file eq "yes")
                    {
                       $code1 = substr $code2, 0, $type_of_file2;
                    }
     
                    if (exists($hash{$code1."1"}))
                    {                    
                        my $value2 = $hash{$code1."1"};
                        delete $hash{$code1."1"};
                        my $code_new1 = ($code_new*10)+1;
                        my $code_new2 = ($code_new*10)+2;
                        my $value3;
                        my $value4;
                       
                        if ($encrypt eq "yes")
                        {
                            $value3 = encrypt $value2;
                            $value4 = encrypt $value;
                        }
                        else
                        {
                            $value3 = $value2;
                            $value4 = $value;
                        }
                        $hash{$code_new} = $value3;
                        $hash{$code_new} .= exists $hash{$code_new} ? ",$value4" : $value4;
 
                        my $first = substr $value, $left, $overlap;
                        my $second = substr $value, -($overlap+$right), $overlap;
                    
                        $hash2b{$first} .= exists $hash2b{$first} ? ",$code_new2" : $code_new2;
                        $hash2c{$second} .= exists $hash2c{$second} ? ",$code_new2" : $code_new2;

                        $first = substr $value2, $left, $overlap;
                        $second = substr $value2, -($overlap+$right), $overlap;
                    
                        $hash2b{$first} .= exists $hash2b{$first} ? ",$code_new1" : $code_new1;
                        $hash2c{$second} .= exists $hash2c{$second} ? ",$code_new1" : $code_new1;
                        $code_new++;
                    }             
                }
            }
            $f = "yes";
        }
        else
        {
            $code = $line;
            $f = "no";
        } 
    }
    $file_count++;
    close INPUT;
}
print "...OK\n";

select(STDERR);
$| = 1;
select(STDOUT); # default
$| = 1;
print "\nRetrieve Seed...";


my $si = '0';
my $space_at_end2;
while (my $line = <INPUT3>)
{
    if ($si eq 0)
    {
        my $last_character = substr $line, -1;       
        if ($last_character =~ m/\s|\t/g)
        {
            $space_at_end2 = "yes";
        }
    }
    if ($si > 0)
    {
        chomp($line);
        if ($space_at_end2 eq "yes")
        {
            chop($line);
        }
        my $seed_input_tmp = $seed_input;
        $seed_input = $seed_input_tmp.$line;
    }
    $si++;
}

my %read1;

FIRST_SEED:

my $seed_input_new2;
my $seed_input_tmp = $seed_input;

if ($bad_read eq "yes" && keys %contigs )
{
    $seed_input_tmp = $read;
}

if ($seed_input_tmp ne "")
{
    my $build = "";
REF0:
    my $n = '0';
    my $overlap_tmp = $overlap;
    if ($build eq "yes2" && $overlap > 42)
    {
        $overlap_tmp = '33';
    }
    while ($n < length($seed_input_tmp) - $overlap_tmp)
    {
        my $first_seed = substr $seed_input_tmp, $n, $overlap_tmp;
        my %first_seed;
        undef %first_seed;
        if ($build eq "yes")
        {
            %first_seed = build_partial2b $first_seed;
        }
        else
        {
            $first_seed{$first_seed} = undef;
        }
        if ($build eq "yes2" && $overlap > '42' && $n < 3000)
        {          
            my @match1 = grep {$_ =~ /.*($first_seed).+$/} keys %hash2b;
            my @match2 = grep {$_ =~ /.*($first_seed).+$/} keys %hash2c;
            my @match3 = (@match1,@match2);
            undef %first_seed;
            foreach $first_seed (@match3)
            {
                $first_seed{$first_seed} = undef;
            }
        }
FIRST_SEED2:foreach my $first_seed (keys %first_seed)
        {
            if (exists($hash2b{$first_seed}))
            {
                my $seed_input_id2 = substr $hash2b{$first_seed}, 1;
                my @seed_input_id = split /,/, $seed_input_id2;
                my $seed_input_id = $seed_input_id[0];
                my $seed_input_id_tmp = substr $seed_input_id, 0, -1;
                my $seed_input_id_end = substr $seed_input_id, -1;
                
                if (exists($hash{$seed_input_id_tmp}))
                {
                    my @seed_input_id_tmp = split /,/,$hash{$seed_input_id_tmp};
                    my $seed_input_new;
                    
                    if ($seed_input_id_end eq "1")
                    {
                        $seed_input_new = $seed_input_id_tmp[0];
                    }
                    elsif ($seed_input_id_end eq "2")
                    {
                        $seed_input_new = $seed_input_id_tmp[1];
                    }
                    if ($encrypt eq "yes")
                    {
                       $seed_input_new2 = decrypt $seed_input_new;
                    }
                    else
                    {
                        $seed_input_new2 = $seed_input_new;
                    }
                    if (exists($read1{$seed_input_new2}))
                    {
                        next FIRST_SEED2;
                    }
                    my $pp = '0';
                    my $pp2= '0';
                    my $part2 = substr $seed_input_new2, -$overlap-8;

                    while ($pp < length($part2)-$overlap)
                    {
                        my $seed_check = substr $part2, $pp, $overlap;
                        my $seed_check_reverse = reverse($seed_check);
                        $seed_check_reverse =~ tr/ACTG/TGAC/;

                        if (exists($hash2b{$seed_check}))
                        {
                            $pp2++;
                        }
                        if (exists($hash2c{$seed_check}))
                        {
                            $pp2++;
                        }
                        if (exists($hash2b{$seed_check_reverse}))
                        {
                            $pp2++;
                        }
                        if (exists($hash2c{$seed_check_reverse}))
                        {
                            $pp2++;
                        }
                        if ($pp2 > 4)
                        {
                            $seed{$seed_input_id} = $seed_input_new2;
                            $seeds_check{$seed_input_id} = undef;
                            $read1{$seed_input_new2} =undef;
                            $contig_count{$seed_input_id} = '0';
                            $position{$seed_input_id} = length($seed{$seed_input_id});
                            $bad_read{$seed_input_id} = "yes";
                            print "...OK\n";
                            print "\nInitial read retrieved successfully: ".$seed_input_new2."\n";
                            print OUTPUT4 "\nInitial read retrieved successfully: ".$seed_input_new2."\n";
                            print OUTPUT5 "\nInitial read retrieved successfully: ".$seed_input_new2."\n";
                            if ($bad_read eq "yes" && keys %contigs)
                            {
                                $noback{$seed_input_id} = "stop";
                                $tree{$id} = $seed_input_id;
                            }
                            else
                            {
                                $tree{"START"} = $seed_input_id;
                                $first_contig_id = $seed_input_id;
                            }
                            if (exists($old_id{$id}))
                            {
                                $old_id{$old_id{$id}} = $seed_input_id;
                            }
                            goto REF2;
                        }  
                        $pp++;
                    }   
                }
            }
            if (exists($hash2c{$first_seed}))
            {
                my $seed_input_id2 = substr $hash2c{$first_seed}, 1;
                my @seed_input_id = split /,/, $seed_input_id2;
                my $seed_input_id = $seed_input_id[0];
                my $seed_input_id_tmp = substr $seed_input_id, 0, -1;
                my $seed_input_id_end = substr $seed_input_id, -1;
                
                if (exists($hash{$seed_input_id_tmp}))
                {
                    my @seed_input_id_tmp = split /,/,$hash{$seed_input_id_tmp};
                    my $seed_input_new;
                    
                    if ($seed_input_id_end eq "1")
                    {
                        $seed_input_new = $seed_input_id_tmp[0];
                    }
                    elsif ($seed_input_id_end eq "2")
                    {
                        $seed_input_new = $seed_input_id_tmp[1];
                    }
                    if ($encrypt eq "yes")
                    {
                       $seed_input_new2 = decrypt $seed_input_new;
                    }
                    else
                    {
                        $seed_input_new2 = $seed_input_new;
                    }
                    if (exists($read1{$seed_input_new2}))
                    {
                        next FIRST_SEED2;
                    }
                    my $pp = '0';
                    my $pp2= '0';
                    my $part2 = substr $seed_input_new2, -$overlap-8;

                    while ($pp < length($part2)-$overlap)
                    {
                        my $seed_check = substr $part2, $pp, $overlap;
                        my $seed_check_reverse = reverse($seed_check);
                        $seed_check_reverse =~ tr/ACTG/TGAC/;

                        if (exists($hash2b{$seed_check}))
                        {
                            $pp2++;
                        }
                        if (exists($hash2c{$seed_check}))
                        {
                            $pp2++;
                        }
                        if (exists($hash2b{$seed_check_reverse}))
                        {
                            $pp2++;
                        }
                        if (exists($hash2c{$seed_check_reverse}))
                        {
                            $pp2++;
                        }
                        if ($pp2 > 4)
                        {
                            $seed{$seed_input_id} = $seed_input_new2;
                            $seeds_check{$seed_input_id} = undef;
                            $read1{$seed_input_new2} =undef;
                            $contig_count{$seed_input_id} = '0';
                            $position{$seed_input_id} = length($seed{$seed_input_id});
                            $bad_read{$seed_input_id} = "yes";
                            print "...OK\n";
                            print "\nInitial read retrieved successfully: ".$seed_input_new2."\n";
                            print OUTPUT4 "\nInitial read retrieved successfully: ".$seed_input_new2."\n";
                            print OUTPUT5 "\nInitial read retrieved successfully: ".$seed_input_new2."\n";
                            if ($bad_read eq "yes" && keys %contigs)
                            {
                                $noback{$seed_input_id} = "stop";
                                $tree{$id} = $seed_input_id;
                            }
                            else
                            {
                                $tree{"START"} = $seed_input_id;
                                $first_contig_id = $seed_input_id;
                            }
                            if (exists($old_id{$id}))
                            {
                                $old_id{$old_id{$id}} = $seed_input_id;
                            }
                            goto REF2;
                        }  
                        $pp++;
                    }                    
                }
            }
        }
        $n++;
        if ($build eq "yes2")
        {
            $n += 4;
        }
    }
    if ($build eq "")
    {
        $build = "yes";
        goto REF0;
    }
    elsif ($build eq "yes")
    {
        $build = "yes2";
        goto REF0;
    }
    else
    {
        print "\n\nINVALID SEED, PLEASE TRY AGAIN WITH A NEW ONE\n\n";
        print OUTPUT4 "\n\nINVALID SEED, PLEASE TRY AGAIN WITH A NEW ONE\n\n";
        exit;
    }
}
else
{
REF:foreach my $line (keys %hash)
    {
        my $count = '0';
        my $n = '5';
        my @line_read = split /,/, $hash{$line};
        my $line_read = $line_read[0];

        while ($n < length($line_read) - $overlap)
        {
            my $line_part = substr $line_read, $n, $overlap;
            if (exists($hash2b{$line_part}))
            {
                my @results = split /,/, $hash2b{$line_part};
                foreach (@results)
                {
                    $count++;
                }
            }
            if (exists($hash2c{$line_part}))
            {
                my @results = split /,/, $hash2c{$line_part};
                foreach (@results)
                {
                    $count++;
                }
            }
            if ($count > 10)
            {
                my $line_read2 = $line_read[1];
                my $count2 = '0';
                my $n2 = '5';
                
                while ($n2 < length($line_read2) - $overlap)
                {
                    my $line_part2 = substr $line_read2, $n2, $overlap;
                    if (exists($hash2b{$line_part2}))
                    {
                        my @results = split /,/, $hash2b{$line_part2};
                        foreach (@results)
                        {
                            $count2++;
                        }
                    }
                    if (exists($hash2c{$line_part2}))
                    {
                        my @results = split /,/, $hash2c{$line_part2};
                        foreach (@results)
                        {
                            $count2++;
                        }
                    }
                    if ($count2 > 10)
                    {                        
                        $seed{$line."1"} = $line_read;
                        $seeds_check{$line."1"} = undef;
                        print "...OK\n";
                        print "\nInitial read retrieved successfully: ".$line_read."\n";
                        print OUTPUT4 "\nInitial read retrieved successfully: ".$line_read."\n";
                        last REF;
                    }          
                    $n2++;   
                }
            }          
            $n++;   
        }
    }
}
close INPUT3;
REF2:
my $seed = "";
my $seed_id = "";
my $read_short_end = "";
my %read_short_end;
my @read_short_end;
my $read_short_end2 = "";
my %read_short_end2;
my $read_short_zone = "";
my $read_short_start = "";
my %read_short_start;
my @read_short_start;
my $read_short_start2 = "";
my %read_short_start2;
my $read_end = "";
my @read_end;
my @read_end_chars;
my $read_end_b = "";
my @read_end_b;
my %read_end;
my %read_end_b;
my %read_start;
my %read_start_b;
my $read_start = "";
my $read_start_b = "";
my @read_start;
my @read_start_b;
my @read_start_chars;
my $read_short_zone_start = "";
my $match_pair = "";
my $match_pair2 = "";
my $extension = "";
my $extension_match = "";
my $match = "";
my $pair = "";
my $regex = "";
my @match_pair;
my %merged_match;
my %merged_match_back;
my %merged_match_pos;
my %merged_match_back_pos;
my @matches;
my @matches1;
my @matches2;
my $read_new = "";
my $id_test = "";
my $best_extension = "";
my $best_match2 = "";
my $best_extension1 = "";
my $best_extension2 = "";
my $best_extension3 = "";
my $best_extension4 = "";
my $use_regex = "";
my $use_regex_back = "";
my %extensions;
my %extensions_original;
my %extensions1;
my %extensions2;
my %extensions1b;
my %extensions2b;
my %extensionsb;
my @extensions;
my @extensions1;
my @extensions2;
my $new_best = "";
my $extra_seed = "";
my $extra_regex = "";
my $position = length($seed_input_new2);
my %insert_size2;
my $insert_size2 = $insert_size;
my $position_back = '0';
my $match_end = "";
my $match_start = "";
my $match_end2 = "";
my $match_start2 = "";
my $noback = "";
my $noforward = "";
my $split = "";
my $best_extension_split;
my $merge = "";
my $X = '0';
my $last_chance = "";
my %last_chance;
my $last_chance_back = "";
my %last_chance_back;
my $circle;
my $containX_short_end2 = "";
my $contain_dot_short_end2 = "";
my $containX_short_start2 = "";
my $contain_dot_short_start2 = "";
my $read_test = "";
my $position_adjust = "";
my $position_adjust_back = "";
my $AT_rich;
my $id_old;
my $enough;
my $enough_back;
my $read_new1;
my $delete_first;
my $delete_second;
my $delete_third;
my $correct_after_split;
my $sc = '0';
my $super_best_extension;
my $still_first_seed = "yes";
my $ext_before;
my %rep_return;
my %rep_return_back;

foreach (keys %seed)
{
    $sc++;    
}
if ($bad_read ne "yes" || $sc eq '1')
{
    foreach $seed_id (keys %seed)
    {
    }
    if ($y > $startprint2)
    {
        print "\nStart Assembly...\n\n";
        print OUTPUT4 "\nStart Assembly...\n\n";
    }
}
ITERATION: while ($y < $iterations)
{
    if ($y > $startprint2)
    {
        print OUTPUT5 "\n".$y."\n\n";
    }
    if (!%seed)
    {
        last ITERATION;
    }
    my $length_other_contig = '0';  
    foreach my $contig_tmp (keys %contigs)
    {
        $length_other_contig += length($contigs{$contig_tmp});
    }

    $|=1;
    my $progress = length($read)+$length_other_contig." bp assembled";
    
    print "\b" x length($progress_before);
    print ' ' x length($progress_before);
    print "\b" x length($progress_before);
    $progress_before = $progress;
    print $progress;
    
    if ($still_first_seed ne "yes")
    {
        $still_first_seed = "yes2";
    }
    else
    {
        $still_first_seed = "";
    }

SEED: foreach $seed_id (keys %seed)
{
    
    if ($still_first_seed ne "yes2")
    {
        my $test_first_seed;
        if ($seed_id =~ m/.*_(\d+)$/)
        {
            $test_first_seed = $1;
        }
        else
        {
            $test_first_seed = $seed_id;
        }
        if ($test_first_seed ne $first_contig_id)
        {
            next SEED;
        }
        $still_first_seed = "yes";
    }
    
    $merge = "";
    $split = "";
    $circle = "";
    $AT_rich = "";

    undef %extensions;
    undef %extensions_original;
    undef %extensions1;
    undef %extensions2;
    undef %extensions1b;
    undef %extensions2b;
    undef @extensions;
    undef @extensions1;
    undef @extensions2;
    undef @matches;
    undef @matches1;
    undef @matches2;    
    undef %merged_match;
    undef %merged_match_pos;
    undef %merged_match_back;
    undef %merged_match_back_pos;
    undef @read_end;
    undef @read_end_b;
    undef %read_end;
    undef %read_end_b;
    undef %read_start;
    undef %read_start_b;
            undef %read_end_tmp;
            undef %read_end_b_tmp;
            undef %read_short_end_tmp;
            undef %read_short_zone_tmp;
    undef %read_short_start_tmp;
    undef %read_short_zone_start_tmp;
    undef %read_start_tmp;
    undef %read_start_b_tmp;
    undef %SNR_length;
    undef %match_rep;
    undef %count_rep;
    undef @read_start;
    undef @read_start_b;
    undef @read_end_chars;
    undef @read_short_end;
    undef @read_start_chars;
    undef @read_short_start;
    undef %filter_before1;
    undef %filter_before2;
    undef %filter_before3;
    undef %filter_before4;
    undef %filter_before1_pair;
    undef %filter_before2_pair;
    undef %filter_before3_pair;
    undef %filter_before4_pair;
    $containX_short_end2 = '0';
    $contain_dot_short_end2 = '0';
    $containX_short_start2 = '0';
    $contain_dot_short_start2 = '0';
    $enough = "";
    $enough_back = "";
    $delete_first = "";
    $delete_second = "";
    $delete_third = "";
    my $merge_read = "";
    my $merge_read_pair = "";
    my $merge_read_length = '0';
    $SNR_read = "";
    $SNR_read2 = "";
    $SNR_read_back = "";
    $SNR_read_back2 = "";
    $split_forward = "";
    $deletion = "";
    $deletion_back = "";
    $reference_guided = "";
    $reference_guided_back = "";
    $repetitive_detect = "";
    $repetitive_detect_back = "";
    $repetitive_detect2 = "";
    $repetitive_detect_back2 = "";
    $contig_end = "";
    $contig_end_check = "";
    $repetitive = "";
    $before_repetitive = "";
    $before_repetitive_short = "";
    $CP_check = "";
    $before_extension1 = "";
    $before_extension2 = "";
    $before_extension_back1 = "";
    $before_extension_back2 = "";
    $id_split1 = "";
    $id_split2 = "";
    $id_split3 = "";
    $super_best_extension = "";
    $no_contig_id2 = "";
    $no_contig_id1 = "";
    $no_contig_id3 = "";
    $no_contig_id4 = "";
    $rep_detect2 = "";
    $no_next_seed = "";
    $count_split = "";
    $yuyu_option_A = "";
    $yuyu_option_C = "";
    $yuyu_option_T = "";
    $yuyu_option_G = "";
    $yuyu_option_A_back = "";
    $yuyu_option_C_back = "";
    $yuyu_option_T_back = "";
    $yuyu_option_G_back = "";
    $extensions_before = "";
    $SNR_next_seed = "";
    $ext_before = "";
    $AT_rich_before = "";
    $AT_rich_before_back = "";
    $SNR_nucleo = "";
    $SNR_nucleo_back = "";
    $insert_range_shorter = "";
    $merge_now = "";
    
    if (exists($indel_split{$seed_id}))
    {
        $indel_split = $indel_split{$seed_id};
        $insert_range = $insert_range_c;
        if ($y > $startprint2)
        {
            print OUTPUT5 "\n".$indel_split." INDEL_SPLIT\n";
        }
    }
    else
    {
        $indel_split = '0';
        $insert_range = $insert_range_b;
    }
    if (exists($indel_split_back{$seed_id}))
    {
        $indel_split_back = $indel_split_back{$seed_id};
        $insert_range_back = $insert_range_c;
        if ($y > $startprint2)
        {
            print OUTPUT5 "\n".$indel_split_back." INDEL_SPLIT_BACK\n";
        }
    }
    else
    {
        $indel_split_back = '0';
        $insert_range_back = $insert_range_b;
    }
    
    if (exists($insert_size2{$seed_id}))
    {
        $insert_size2 = $insert_size2{$seed_id};
    }
    
    if (exists($seed{$seed_id}))
    {
        $seed = $seed{$seed_id};
           
        if ($y > $startprint2)
        {
            print OUTPUT5 "\n".$seed_id." SEED_exists\n\n";
            print OUTPUT5 length($seed)." READ_LENGTH\n";
        }
        
        $id = $seed_id;
        $read = $seed;
        

        my $SNR_end0 = substr $read, -20, 20;
        my $SNR_end0t = substr $read, -$overlap, $overlap;
        $SNR_end0 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
        $SNR_end0t =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
        my @SNR_end0 = split //, $SNR_end0;
        my $u0 = length($SNR_end0);
        my $v0 = '1';
        my $other = "";
        my $Spn0 = '-1';
                
ALREADY_X0:  while ($v0 < $u0)
        {              
            if ($SNR_end0[$u0-$v0-1] eq $SNR_end0[$u0-$v0] || $SNR_end0[$u0-$v0-1] eq ".")
            {
                $v0++;
                if(($v0 > 7 && $other eq "") || $v0 > 10)
                {
                    $SNR_read = "yes";
                    $SNR_nucleo = substr $read, $Spn0, 1;
                    if ($SNR_end0[$u0-1] eq $SNR_nucleo)
                    {
                        $SNR_read2 = "yes";
                    }
                    last ALREADY_X0;
                }
            }
            elsif($SNR_end0[$u0-$v0-1] eq "X")
            {
                if ($y > $startprint2)
                {
                    print OUTPUT5 "ALREADY_X\n";
                }  
                $SNR_read = "X";
                last ALREADY_X0;
            }
            elsif(($v0 > 7 && $other eq "") || $v0 > 10)
            {     
                $SNR_read = "yes";
                if (length($best_extension_prev{$id}) > 1)
                {
                    $SNR_nucleo = substr $read, $Spn0, 1;
                }
                else
                {
                    $SNR_nucleo = substr $read, $Spn0, 1;
                }  
                if ($SNR_end0[$u0-1] eq $SNR_nucleo)
                {
                    $SNR_read2 = "yes";
                }
                print OUTPUT5 $SNR_end0[$u0-1]." ALREADY_X\n";
                last ALREADY_X0;
            }
            elsif($other eq "")
            {     
                $other = "yes";
                if ($v0 eq '1')
                {
                    $Spn0 = '-3';
                }     
                $v0++;
                $v0++;
            }
            else
            {     
                $SNR_read = "";
                last ALREADY_X0;
            }
        }
        if ($SNR_read eq "")
        {
            my $SNR_check = $SNR_end0t =~ s/AAAAAAAA|CCCCCCCC|GGGGGGGG|TTTTTTTT|TATATATATATATATA//;
            if ($SNR_check > 0)
            {
                $SNR_read = "yes"; 
            }
        }
        
        my $SNR_end0b = substr $read, 0, 20;
        my $SNR_end0bt = substr $read, 0, $overlap;
        $SNR_end0b =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
        $SNR_end0bt =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
        my @SNR_end0b = split //, $SNR_end0b;
        my $u0b = length($SNR_end0b);
        my $v0b = '1';
        my $otherb = "";
        my $Spn0b = '0';
                
ALREADY_X0b:  while ($v0b < $u0b)
        {              
            if ($SNR_end0b[$v0b-1] eq $SNR_end0b[$v0b] || $SNR_end0b[$v0b] eq ".")
            {
                $v0b++;
                if(($v0b > 7 && $otherb eq "") || $v0b > 9)
                {     
                    $SNR_read_back = "yes";
                    $SNR_nucleo_back = substr $read, $Spn0b, 1; 
                    if ($SNR_end0b[0] eq $SNR_nucleo_back)
                    {
                        $SNR_read_back2 = "yes";
                    }
                    last ALREADY_X0b;
                }
            }
            elsif($SNR_end0b[$v0b] eq "X")
            {
                if ($y > $startprint2)
                {
                    print OUTPUT5 "ALREADY_Xb\n";
                }  
                $SNR_read_back = "X";
                last ALREADY_X0b;
            }
            elsif(($v0b > 7 && $otherb eq "") || $v0b > 9)
            {     
                $SNR_read_back = "yes";
                if (length($best_extension_back_prev{$id}) > 1)
                {
                    $SNR_nucleo_back = substr $read, $Spn0b, 1;
                }
                else
                {
                    $SNR_nucleo_back = substr $read, $Spn0b, 1;
                }  
                if ($SNR_end0b[0] eq $SNR_nucleo_back)
                {
                    $SNR_read_back2 = "yes";
                }
                last ALREADY_X0b;
            }
            elsif($otherb eq "")
            {     
                $otherb = "yes";
                if ($v0b eq '1')
                {
                    $Spn0b = '2';
                }  
                $v0b++;
                $v0b++;
            }
            else
            {     
                $SNR_read_back = "";
                last ALREADY_X0b;
            }
        }
        if ($SNR_read_back eq "")
        {
            my $SNR_check = $SNR_end0bt =~ s/AAAAAAAA|CCCCCCCC|GGGGGGGG|TTTTTTTT|TATATATATATATATA//;
            if ($SNR_check > 0)
            {
                $SNR_read_back = "yes";  
            }
        }
        
        $contig_count = $contig_count{$id};
        
        if (exists($position{$id}))
        {
            $position = $position{$id};
        }
        if (exists($position_back{$id}))
        {
            $position_back = $position_back{$id};
        }
        
        $read_short_end = substr $read, -($insert_size*$insert_range)+(($read_length-$overlap-8)/2), ((($insert_size*$insert_range)-$insert_size)*2)+$overlap+10+8;
        $read_short_end2 = substr $read, -350;
        $read_short_start = substr $read, ($insert_size*$insert_range_back)-(($read_length-$overlap-8)/2) - (((($insert_size*$insert_range_back)-$insert_size)*2)+$overlap+10+8), ((($insert_size*$insert_range_back)-$insert_size)*2)+$overlap+10+8;
        $read_short_start2 = substr $read, 0, 350;
        
        $read_short_end =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
        $read_short_end2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
        $read_short_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
        $read_short_start2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
               
        $containX_short_end2 = $read_short_end2 =~ tr/X|\*//;
        $contain_dot_short_end2 = $read_short_end2 =~ tr/\.//;
        $containX_short_start2 = $read_short_start2 =~ tr/X|\*//;
        $contain_dot_short_start2 = $read_short_start2 =~ tr/\.//;
        
        if ($y eq '1' || exists($old_id2{$id}) || exists($bad_read{$id}))
        {
            delete $bad_read{$id};

            $position{$id} = length($read);
            $position = length($read);
            $position_back = '0';
            $position_back{$id} = '0';

            if ($y eq '1' || ($correct_after_split eq "yesssss" && length($read) <= $read_length+1))
            {   
                $read = correct ($read);
                $read =~ s/\s+$//;
                $read =~ s/\t+$//;
            }
            delete $old_id2{$id};
        }
        if ($y eq '2')
        {
            my $start_point = '25'; 
            $first_contig_start_reverse = substr $read, $start_point, $overlap;

            my $check_start = $first_contig_start_reverse =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
            while ($check_start > 0)
            {
                $start_point += 5;
                $first_contig_start_reverse = substr $read, $start_point, $overlap;
                $check_start = $first_contig_start_reverse =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
            }
            $first_contig_start_reverse = reverse($first_contig_start_reverse);
            $first_contig_start_reverse =~ tr/ATCG/TAGC/;
        }     
            
        if ($y > $startprint2)
        {
            print OUTPUT5 $insert_size." INSERT_SIZE\n\n";
            print OUTPUT5 $position." POSITION\n";
            print OUTPUT5 $position_back." POSITION_BACK\n";
        }
        if ($use_regex eq "")
        {
            my $read_end_dot = substr $read_short_end2, -($read_length-($overlap+$left+1-$containX_short_end2-$containX_short_end2));
            my $read_end_dot_check = $read_end_dot =~ tr/\./\./;
            if ($read_end_dot_check > 0)
            {
                $use_regex = "yes2";
            }
        }
        if ($use_regex_back eq "")
        {
            my $read_start_dot = substr $read_short_start2, 0, ($read_length-($overlap+$left+1-$containX_short_end2-$containX_short_end2));
            my $read_start_dot_check = $read_start_dot =~ tr/\./\./;
            if ($read_start_dot_check > 0)
            {
                $use_regex_back = "yes2";
            }
        }        

        if ($SNR_read ne "")
        {
            if ($y > $startprint2)
            {
                print OUTPUT5 $SNR_nucleo." SNR_READ\n";
            }
        }
        if ($SNR_read2 ne "")
        {
            if ($y > $startprint2)
            {
                print OUTPUT5 "SNR_READ2\n";
            }
        }
        if ($SNR_read_back ne "")
        {
            if ($y > $startprint2)
            {
                print OUTPUT5 "SNR_BACK_READ\n";
            }
        }
        if (exists($regex{$id}))
        {
            $use_regex = $regex{$id};
        }
        elsif ($use_regex ne "yes2")
        {
            $use_regex = "";
        }
        if (exists($before{$id}))
        {
            $before = $before{$id};
            if ($y > $startprint2)
            {
                print OUTPUT5 "BEFORE\n";
            } 
        }
        else
        {
            $before = "";
        }
        if (exists($before_back{$id}))
        {
            $before_back = $before_back{$id};
            if ($y > $startprint2)
            {
                print OUTPUT5 "BEFORE_BACK\n";
            }
        }
        else
        {
            $before_back = "";
        }
        if (exists($regex_back{$id}))
        {
            $use_regex_back = $regex_back{$id};
        }
        elsif ($use_regex_back ne "yes2")
        {
            $use_regex_back = "";
        }
        if (exists($noforward{$id}))
        {
            $noforward = "stop";
        }
        else
        {
            $noforward = "";
        }
        my $check_old_id = "";
        if (exists($old_id{$id}))
        {
            $check_old_id = "yes";
        }
        if (exists($noback{$id}) || (length($read) < $insert_size && $noforward eq "" && $check_old_id eq ""))
        {
            $noback = "stop";
        }
        else
        {
            $noback = "";
        }         
        if (exists($last_chance{$id}) || length($read) < $insert_size)
        {
            $last_chance = "yes";
            if ($y > $startprint2)
            {
                print OUTPUT5 "LAST_CHANCE\n";
            }
        }
        else
        {
            $last_chance = "";
        }
        if (exists($last_chance_back{$id}))
        {
            $last_chance_back = $last_chance_back{$id};
            if ($y > $startprint2)
            {
                print OUTPUT5 "LAST_CHANCE_BACK\n";
            }
        }
        else
        {
            $last_chance_back = "";
        }
        if (exists($indel_split_skip{$id}))
        {
            $indel_split_skip = $indel_split_skip{$id};
            if ($y > $startprint2)
            {
                print OUTPUT5 "INDEL_SPLIT_SKIP\n";
            }
        }
        else
        {
            $indel_split_skip = "";
        }
        if (exists($indel_split_skip_back{$id}))
        {
            $indel_split_skip_back = $indel_split_skip_back{$id};
            if ($y > $startprint2)
            {
                print OUTPUT5 "INDEL_SPLIT_SKIP_BACK\n";
            }
        }
        else
        {
            $indel_split_skip_back = "";
        }
        if (exists($position_adjust{$id}))
        {
            $position_adjust = $position_adjust{$id};
        }
        else
        {
            $position_adjust = "";
        }
        if (exists($SNP_active{$id}))
        {
            $SNP_active = "yes";
            if ($y > $startprint2)
            {
                print OUTPUT5 "SNP_ACTIVE\n";
            }
        }
        else
        {
            $SNP_active = "";
        }
        if (exists($SNP_active_back{$id}))
        {
            $SNP_active_back = "yes";
            if ($y > $startprint2)
            {
                print OUTPUT5 "SNP_ACTIVE_BACK\n";
            }
        }
        else
        {
            $SNP_active_back = "";
        }
        if (exists($nosecond{$id}))
        {
            $nosecond = "yes";
        }
        else
        {
            $nosecond = "";
        }
        if (exists($before_shorter_skip{$id}))
        {
            $before_shorter_skip = "yes";
            if ($y > $startprint2)
            {
                print OUTPUT5 "BEFORE SHORTER\n";
            }
        }
        else
        {
            $before_shorter_skip = "";
        }
        if (exists($before_shorter_skip_back{$id}))
        {
            $before_shorter_skip_back = "yes";
            if ($y > $startprint2)
            {
                print OUTPUT5 "BEFORE SHORTER BACK\n";
            }
        }
        else
        {
            $before_shorter_skip_back = "";
        }
        if (exists($jump_rep{$id}))
        {
            $jump_rep = "yes";
            if ($y > $startprint2)
            {
                print OUTPUT5 "JUMP REP\n";
            }
        }
        else
        {
            $jump_rep = "";
        }
        if (exists($jump_rep_back{$id}))
        {
            $jump_rep_back = "yes";
            if ($y > $startprint2)
            {
                print OUTPUT5 "JUMP REP BACK\n";
            }
        }
        else
        {
            $jump_rep_back = "";
        }

        if (exists($yuyu_option{$id.'A'}))
        {
            $yuyu_option_A = "A";                        
        }
        if (exists($yuyu_option{$id.'C'}))
        {
            $yuyu_option_C = "C";         
        }
        if (exists($yuyu_option{$id.'G'}))
        {
            $yuyu_option_G = "G";         
        }
        if (exists($yuyu_option{$id.'T'}))
        {
            $yuyu_option_T = "T";         
        }
        
        if (exists($yuyu_option_back{$id.'A'}))
        {
            $yuyu_option_A_back = "A";                        
        }
        if (exists($yuyu_option_back{$id.'C'}))
        {
            $yuyu_option_C_back = "C";         
        }
        if (exists($yuyu_option_back{$id.'G'}))
        {
            $yuyu_option_G_back = "G";         
        }
        if (exists($yuyu_option_back{$id.'T'}))
        {
            $yuyu_option_T_back = "T";         
        }

        if ($noback ne "stop")
        {
            my $start_repetitive = substr $read, $overlap, $insert_size+100;
            my $repetitive_test = substr $read_short_start2, 0, 15;
            my $SNR_skip = "yes";                                        
            my $u = '0';
            while ($SNR_skip eq "yes")
            { 
                $repetitive_test = substr $read_short_start2, $u, 15;
                my $SNR_checkA = $repetitive_test =~ tr/A/A/;
                my $SNR_checkC = $repetitive_test =~ tr/C/C/;
                my $SNR_checkT = $repetitive_test =~ tr/T/T/;
                my $SNR_checkG = $repetitive_test =~ tr/G/G/;
                my $SNR_checkdot = $repetitive_test =~ tr/\./\./;
                if ($SNR_checkA+$SNR_checkdot > 12 || $SNR_checkC+$SNR_checkdot > 12 || $SNR_checkG+$SNR_checkdot > 12 || $SNR_checkT+$SNR_checkdot > 12)
                {
                    $SNR_skip = "yes";
                    $u += 10;
                }
                else
                {
                    $SNR_skip = "";
                }
            }
            $repetitive_test =~ tr/\*//d;
            my $check_repetitive = $start_repetitive =~ s/$repetitive_test/$repetitive_test/g;
            
            if (exists ($rep_return_back{$id}))
            {
                $repetitive_detect_back = "yes2";
                if ($rep_return_back{$id} < $position_back-(2*$read_length))
                {
                    delete $rep_return_back{$id};
                }
                if ($y > $startprint2)
                {
                    print OUTPUT5 "DETECT_REPETITIVE_BACK_RETURN\n";
                    print OUTPUT5 $start_repetitive." START_READ\n";
                }
            }
            if ($check_repetitive > 0)
            {
                if (exists ($rep_return_back{$id}))
                {
                    
                }
                elsif (length($read) > $insert_size+100)
                {
                    substr $read, 0, $read_length, "";
                    $seed{$id} = $read;
                    $seeds_check{$id} = undef;
                    $position_back{$id} = $position_back-$read_length;
                    $rep_return_back{$id} = $position_back;
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "RETURN_REPETITIVE_BACK\n";
                    }
                    goto SEED;
                }
                $repetitive_detect_back = "yes";
                if ($y > $startprint2)
                {
                    print OUTPUT5 "DETECT_REPETITIVE_back\n";
                    print OUTPUT5 $start_repetitive." START_READ\n";
                }
                if ($check_repetitive > 1)
                {
                    my $start_repetitive1 = substr $read, $insert_size+50, 300;
                    my $check_repetitive1= $start_repetitive1 =~ s/$repetitive_test/$repetitive_test/g;
                    if ($check_repetitive > 1)
                    {
                        $repetitive_detect_back2 = "yes";
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "DETECT_REPETITIVE_back2\n";
                        }
                    }
                }
            }
        }
        if ($noforward ne "stop")
        {
            my $end_repetitive = substr $read, -$insert_size-100,-$overlap;
            my $repetitive_test2 = substr $read_short_end2, -15;
            my $SNR_skip = "yes";                                        
            my $u = '15';
            while ($SNR_skip eq "yes")
            { 
                $repetitive_test2 = substr $read_short_end2, -$u, 15;
                my $SNR_checkA = $repetitive_test2 =~ tr/A/A/;
                my $SNR_checkC = $repetitive_test2 =~ tr/C/C/;
                my $SNR_checkT = $repetitive_test2 =~ tr/T/T/;
                my $SNR_checkG = $repetitive_test2 =~ tr/G/G/;
                my $SNR_checkdot = $repetitive_test2 =~ tr/\./\./;
                if ($SNR_checkA+$SNR_checkdot > 12 || $SNR_checkC+$SNR_checkdot > 12 || $SNR_checkG+$SNR_checkdot > 12 || $SNR_checkT+$SNR_checkdot > 12)
                {
                    $SNR_skip = "yes";
                    $u += 10;
                }
                else
                {
                    $SNR_skip = "";
                }
            }
            $repetitive_test2 =~ tr/\*//d;
            my $check_repetitive2 = $end_repetitive =~ s/$repetitive_test2/$repetitive_test2/g;

            if (exists ($rep_return{$id}))
            {
                $repetitive_detect = "yes2";
                if ($y > $startprint2)
                {
                    print OUTPUT5 "DETECT_REPETITIVE_RETURN\n";
                    print OUTPUT5 $end_repetitive." END_READ\n";
                }
                if ($rep_return{$id} < $position-(2*$read_length))
                {
                    delete $rep_return{$id};
                }
            }
            if ($check_repetitive2 > 0)
            {
                if (exists ($rep_return{$id}))
                {
                    
                }
                elsif (length($read) > $insert_size+100)
                {
                    substr $read, -$read_length, $read_length, "";
                    $seed{$id} = $read;
                    $seeds_check{$id} = undef;
                    $position{$id} = $position-$read_length;
                    $rep_return{$id} = $position;
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "RETURN_REPETITIVE\n";
                    }
                    goto SEED;
                }
                $repetitive_detect = "yes";
                if ($y > $startprint2)
                {
                    print OUTPUT5 "DETECT_REPETITIVE\n";
                    print OUTPUT5 $end_repetitive." END_READ\n";
                }
                if ($check_repetitive2 > 1)
                {
                    my $end_repetitive1 = substr $read, -$insert_size-350, 300;
                    my $check_repetitive21 = $end_repetitive1 =~ s/$repetitive_test2/$repetitive_test2/g;
                    if ($check_repetitive21 > 1)
                    {
                        $repetitive_detect2 = "yes";
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "DETECT_REPETITIVE2\n";
                        }
                    }
                }
            }
        }
MERGE:
        my $merge_extra = '0';
        if ($y > $startprint2)
        {    
            if (exists($old_id{$id}))
            {
                print OUTPUT5 $old_id{$id}." OLD_ID\n";
            }
            if (exists($old_rep{$id}))
            {
                print OUTPUT5 "OLD_REP\n";
            }
        }
        if (exists $old_id{$id} && exists $old_rep{$id})
        {
            if ($position_back > $read_length && $position_back < ($insert_size*3) && $repetitive_detect_back eq "" && $noback ne "stop")
            {
                my $read_oldie = $seed_old{$old_id{$id}};
                my $read_newest = $read;
                                                            
                my $start_seq = substr $read, 0, $insert_size*3;
                my $start_seq1 = substr $read, 0, 39;
                my $start_seq2 = substr $read, 30, 39;
                my $end_seq = substr $read_oldie, -$insert_size*3;
                $merge_read_length = length ($read);
                $start_seq =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                $start_seq1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                $start_seq2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                $end_seq =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
      
                if ($end_seq =~ m/.*.$start_seq1(.*)$/)
                {
                    my $r = length($1);
                    my $read_temp = $read;
                    my $read1 = substr $read_oldie, 0, -39-$r;         
    
                    $read = $read1.$read_temp;                                                      
                    $merge_read = "yes";
                }
                elsif ($end_seq =~ m/.*.$start_seq2(.*)$/)
                {
                    my $r = length($1);
                    my $read_temp = $read;
                    my $read1 = substr $read_oldie, 0, -39-$r;         
                    my $read2 = substr $read_temp, 30;
                    
                    $read = $read1.$read2;                                      
                    $merge_read = "yes";
                }
                if ($merge_read eq "yes")
                {            
                    $seed{$id} = $read;
                    $seeds_check{$id} = undef;
                    delete $tree{$old_id{$id}};
                    
                    foreach my $tree_tmp (keys %tree)
                    {
                        my $old = $old_id{$id};
                        my $tree2 = $tree{$tree_tmp};
                        my $tree3 = $tree{$tree_tmp};

                        if ($old_id{$id} =~ m/.*_(\d+)$/)
                        {
                            $old = $1;
                        }
                        my $id_tmp = $id;
                        if ($id =~ m/.*_(\d+)$/)
                        {
                            $id_tmp = $1;
                        }
                        my @ids_split = split /\*/, $tree2;
                        foreach my $id_split (@ids_split)
                        {
                            if ($id_split  =~ m/^$old(REP)*$/)
                            {
                                if ($tree2 =~ m/^(.*\*)*$old(REP)*(\*.*)*$/)
                                {
                                    if (defined($1))
                                    {
                                        $tree3 = $1.$id_tmp;
                                    }
                                    else
                                    {
                                        $tree3 = $id_tmp;
                                    }
                                    if (defined($2))
                                    {
                                        $tree3 = $tree3."REP";
                                    }
                                    if (defined($3))
                                    {
                                        $tree3 = $tree3.$3;
                                    }
                                }
                            }
                        }
                        delete $tree{$tree_tmp};
                        $tree{$tree_tmp} = $tree3;
                        foreach my $contigs_end (keys %contigs_end)
                        {
                            if ($contigs_end{$contigs_end} eq $old)
                            {
                                delete $contigs_end{$contigs_end};
                                $contigs_end{$contigs_end} = $id_tmp;
                            }
                        }
                    }
                    delete $old_id{$id};
                    delete $old_rep{$id};
                    delete $old_rep_old{$old_id{$id}};
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "Merged both contigs rep!\n";
                        print OUTPUT5 ">".$read_newest."\n";
                    }            
                    $noback = "stop";
                    $noback{$id} = "stop";
                    $merge_extra = $position;
                }
            }
            elsif ((($noback eq "stop" || $position_back >= ($insert_size*3)) && $position > $insert_size+200) || $merge_now ne "")
            {
                $tree{$old_id{$id}} = $id."REP";
                $noback{$id} = "stop";
                $noback = "stop";
                $merge_extra = $position;
                print OUTPUT5 ">Remove rep merge!\n";
                 
                if ($contig_num eq '1')
                {
                    $contigs{$contig_num."+".$old_id{$id}} = $seed_old{$old_id{$id}};
                    my $start_point = '500';
                    my $check_repetitive = '3';
                                                                                
                    while ($check_repetitive > 2)
                    {
                        my $repetitive = substr $seed_old{$old_id{$id}}, $start_point, 15;
                        $repetitive =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                        my $read_short_area = substr $seed_old{$old_id{$id}}, $start_point -170, 340;
                        $check_repetitive = $read_short_area =~ s/$repetitive/$repetitive/g;
                        if ($check_repetitive > 2)
                        {
                            $start_point += 20;
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "DETECT_REPETITIVE_IN_START_SEQUENCE\n";
                            }
                        } 
                    }
                    $first_contig_start = substr $seed_old{$old_id{$id}}, $start_point, $overlap;
                    my $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                    while ($check_start > 0)
                    {
                        $start_point += 10;
                        $first_contig_start = substr $seed_old{$old_id{$id}}, $start_point, $overlap;
                        $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                    }
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $first_contig_start." CONTIG_START2\n";
                    }
                    foreach my $seedie (keys %seed)
                    {
                        my $seedie_part = substr $seed{$seedie}, 0, 1000;
                        my $check_seedie = $seedie_part =~ s/$first_contig_start/$first_contig_start/;
                        if ($check_seedie > 0)
                        {
                            delete $seed{$seedie};
                        }
                    }
                }
                else
                {
                    $contigs{$contig_num."+".$old_id{$id}} = $seed_old{$old_id{$id}};
                }
                $contig_num++;
                delete $old_id{$id};
                delete $old_rep{$id};
                delete $old_rep_old{$old_id{$id}};
                if ($noback eq "stop" && $noforward eq "stop" && $merge_now eq "yes2")
                {
                    $contigs{$contig_num."+".$id} = $read;
                    $contig_num++;
                }
            }
        }
        else
        {
            if (exists $old_id{$id} && $position_back > 25 && $position_back < ($insert_size*3) && $noback ne "stop")
            {                                    
                my $read_oldie = $seed_old{$old_id{$id}};
                my $read_newest = $read;
                                                            
                my $start_seq = substr $read, 0, $insert_size+200;
                my $start_seq1 = substr $read, 0, 39;
                my $start_seq2 = substr $read, 30, 39;
                my $end_seq = substr $read_oldie, -$insert_size+200;
                my $end_seq1 = substr $read_oldie, -39, 39;
                my $end_seq2 = substr $read_oldie, -72, 39;
                $merge_read_length = length ($read);
                $start_seq =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                $start_seq1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                $start_seq2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                $end_seq =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                $end_seq1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                $end_seq2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                
      
                if ($end_seq =~ m/.*.$start_seq1(.*)$/)
                {
                    my $r = length($1);
                    my $read_temp = $read;
                    my $read1 = substr $read_oldie, 0, -39-$r;         
    
                    $read = $read1.$read_temp;                                                      
                    $merge_read = "yes";
                }
                elsif ($end_seq =~ m/.*.$start_seq2(.*)$/)
                {
                    my $r = length($1);
                    my $read_temp = $read;
                    my $read1 = substr $read_oldie, 0, -39-$r;         
                    my $read2 = substr $read_temp, 30;
                    
                    $read = $read1.$read2;                                      
                    $merge_read = "yes";
                }
                elsif ($start_seq =~ m/(.*.)$end_seq1.*$/)
                {
                    my $r = length($1);
                    my $read_temp = $read;
                    my $read1 = substr $read_oldie, 0, -39-$r;         
                    
                    $read = $read1.$read_temp;                                      
                    $merge_read = "yes";
                }
                elsif ($start_seq =~ m/(.*.)$end_seq2.*$/)
                {
                    my $r = length($1);
                    my $read_temp = $read;
                    my $read1 = substr $read_oldie, 0, -39-$r;         
                    my $read2 = substr $read_temp, 30;
                    
                    $read = $read1.$read2;                                      
                    $merge_read = "yes";
                }
                
                if ($merge_read eq "yes")
                {            
                    $seed{$id} = $read;
                    $seeds_check{$id} = undef;
                    delete $tree{$old_id{$id}};
                    foreach my $tree_tmp (keys %tree)
                    {
                        my $old = $old_id{$id};
                        my $tree2 = $tree{$tree_tmp};
                        my $tree3 = $tree{$tree_tmp};
                        if ($old_id{$id} =~ m/.*_(\d+)$/)
                        {
                            $old = $1;
                        }
                        my $id_tmp = $id;
                        if ($id =~ m/.*_(\d+)$/)
                        {
                            $id_tmp = $1;
                        }
                        my @ids_split = split /\*/, $tree2;
                        foreach my $id_split (@ids_split)
                        {
                            if ($id_split  =~ m/^$old(REP)*$/)
                            {
                                if ($tree2 =~ m/^(.*\*)*$old(REP)*(\*.*)*$/)
                                {
                                    if (defined($1))
                                    {
                                        $tree3 = $1.$id_tmp;
                                    }
                                    else
                                    {
                                        $tree3 = $id_tmp;
                                    }
                                    if (defined($2))
                                    {
                                        $tree3 = $tree3."REP";
                                    }
                                    if (defined($3))
                                    {
                                        $tree3 = $tree3.$3;
                                    }
                                }
                            }
                        }
                        delete $tree{$tree_tmp};
                        $tree{$tree_tmp} = $tree3;
                        foreach my $contigs_end (keys %contigs_end)
                        {
                            if ($contigs_end{$contigs_end} eq $old)
                            {
                                delete $contigs_end{$contigs_end};
                                $contigs_end{$contigs_end} = $id_tmp;
                            }
                        }
                    }
                    delete $old_id{$id};
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "Merged both contigs!\n";
                        print OUTPUT5 ">".$read_newest."\n";
                    }            
                    $noback = "stop";
                    $noback{$id} = "stop";
                    $merge_extra = $position;
                }
            }
            if (((exists $old_id{$id} && ($noback eq "stop" || $position_back >= ($insert_size*3)) && $position > 1000) || $merge_now ne "") && $merge_read eq "")
            {                                    
                $merge_read_length = length($read);
                $merge_read = "yes"; 
                $read = $seed_old{$old_id{$id}} ."LLLLLLLLLLLLLLL".$read;
                $seed{$id} = $read;
                $seeds_check{$id} = undef;
                $hasL = "yes";
                foreach my $tree_tmp (keys %tree)
                {
                    my $old = $old_id{$id};
                    my $tree2 = $tree{$tree_tmp};
                    my $tree3 = $tree{$tree_tmp};
                    if ($old_id{$id} =~ m/.*_(\d+)$/)
                    {
                        $old = $1;
                    }
                    my $id_tmp = $id;
                    if ($id =~ m/.*_(\d+)$/)
                    {
                        $id_tmp = $1;
                    }
                    my @ids_split = split /\*/, $tree2;
                    foreach my $id_split (@ids_split)
                    {
                        if ($id_split  =~ m/^$old(REP)*$/)
                        {
                            if ($tree2 =~ m/^(.*\*)*$old(REP)*(\*.*)*$/)
                            {
                                if (defined($1))
                                {
                                    $tree3 = $1.$id_tmp;
                                }
                                else
                                {
                                    $tree3 = $id_tmp;
                                }
                                if (defined($2))
                                {
                                    $tree3 = $tree3."REP";
                                }
                                if (defined($3))
                                {
                                    $tree3 = $tree3.$3;
                                }
                            }
                        }
                    }
                    delete $tree{$tree_tmp};
                    $tree{$tree_tmp} = $tree3;
                    foreach my $contigs_end (keys %contigs_end)
                    {
                        if ($contigs_end{$contigs_end} eq $old)
                        {
                            delete $contigs_end{$contigs_end};
                            $contigs_end{$contigs_end} = $id_tmp;
                        }
                    }
                }
                delete $old_id{$id};
                $noback{$id} = "stop";
                $noback = "stop";
                $merge_extra = $position;
                if ($y > $startprint2)
                {
                    print OUTPUT5 ">Merged contigs with LLLLLLLLLLL!\n";
                }
                $contig_gap_min{$id."_".$contig_count} = ($contig_gap_min{$id."_".$contig_count}-$position_back);
                $contig_gap_max{$id."_".$contig_count} = ($contig_gap_max{$id."_".$contig_count}-$position_back);
                if ($noback eq "stop" && $noforward eq "stop" && $merge_now eq "yes2")
                {
                    $contigs{$contig_num."+".$id} = $read;
                    $contig_num++;
                }
            }
        }
        if ($merge_now ne "")
        {
            if ($noforward eq "stop")
            {  
                print OUTPUT5 ">MERGE_NOW ".$id."\n";
                print OUTPUT5 $read."\n";
                delete $seed{$id};
                if (!keys %seed)
                {
                    $circle = "contigs";
                    goto FINISH;
                }
                else
                {
                    goto ITERATION;
                }
            }
            else
            {
                print OUTPUT5 ">\n";
                $noback = "stop";
                $noback{$id} = "stop";
                $seed_split{$id} = undef;
                $best_extension = "";
                goto FINISH;
            }
        }
        

    if (exists $old_id{$id})
    {
            
    }
    elsif (keys %contigs && $first_contig_start ne "")
    {
        my $check_start0 = '0';
        if ($merge_extra ne "")
        {
            my $long_end = substr $read, -1250-$merge_extra;
            $check_start0 = $long_end =~ s/$first_contig_start/$first_contig_start/;
        }     
        my $check_start1 = $read_short_end2 =~ s/$first_contig_start/$first_contig_start/;
        my $check_start2 = $read_short_end2 =~ s/$first_contig_start_reverse/$first_contig_start_reverse/;
        if ($check_start1 > 0 || $check_start0 > 0)
        {
            $tree{$id} = "END";
        }
        if ($check_start2 > 0)
        {
            $tree{$id} = "END_REVERSE";
        }
        if ($check_start1 > 0 || $check_start2 > 0 || $check_start0 > 0)
        {
            $noforward{$id} = "stop";
            $noforward = "stop";
            if ($y > $startprint2)
            {                                                                                 
                if ($check_start1 > 0 || $check_start0 > 0)
                {
                    print OUTPUT5 "\nSTOP_CONTIG, encouter start sequence\n\n";
                }
                if ($check_start2 > 0)
                {
                    print OUTPUT5 "\nSTOP_CONTIG, encouter start sequence reverse\n\n";
                }
                print OUTPUT5 ">".$id."\n";
                print OUTPUT5 $read."\n";
            }

            delete $seed{$id};                                         
            if ($check_start1 > 0 || $check_start0 > 0)
            {
                $contigs{$contig_num."+".$id} = $read;
                $contig_num++;
            } 
            goto SEED;
        }
    }
        
REPEAT:        
        $read_new = $read;
        $read_new1 = $read;
        
        if ($y > $startprint2)
        {
            if ($use_regex eq "yes")
            {
                print OUTPUT5 "USE_REGEX\n";
            }
            if ($use_regex eq "yes2")
            {
                print OUTPUT5 "USE_REGEX2\n";
            }
        }
        if ($y > $startprint2)
        {
            if ($use_regex_back eq "yes2")
            {
                print OUTPUT5 "USE_REGEX_BACK2\n";
            }
            if ($use_regex_back eq "yes")
            {
                print OUTPUT5 "USE_REGEX_BACK2\n";
            }
        }


                                                
                                                if (length($read) > $genome_range_low)
                                                {
                                                    my $start_seq = substr $read_new, 0, 200;
                                                    my $start_seq1 = substr $read_new, 30, 42;
                                                    my $start_seq2 = substr $read_new, 60, 42;
                                                    my $start_seq3 = substr $read_new, 90, 42;
                                                    my $end_seq = substr $read_new, -200;
                                                    my $end_seq1 = substr $read_new, -72, 42;
                                                    my $end_seq2 = substr $read_new, -102, 42;
                                                    my $end_seq3 = substr $read_new, -132, 42;
                                                    my $end_seq1_merge = "";
                                                    my $end_seq2_merge = "";
                                                    my $end_seq3_merge = "";
                                                    $start_seq =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    $start_seq1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    $start_seq2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    $start_seq3 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    $end_seq =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    $end_seq1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    $end_seq2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    $end_seq3 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    
                                                    if ($merge_read eq "yes")
                                                    {
                                                        $end_seq1_merge = substr $read_new, -$merge_read_length+20, 42;
                                                        $end_seq2_merge = substr $read_new, -$merge_read_length+100, 42;
                                                        $end_seq3_merge = substr $read_new, -$merge_read_length+170, 42;
                                                        $end_seq1_merge =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                        $end_seq2_merge =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                        $end_seq3_merge =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    }
                                                    
                                                    if ($end_seq =~ m/.*.$start_seq1(.*)$/)
                                                    {
                                                        my $r = length($1);
                                                        my $read_new_temp = $read_new;
                                                        if ($r eq '0')
                                                        {
                                                            $read_new = substr $read_new_temp, 30+42;             
                                                        }
                                                        else
                                                        {
                                                            $read_new = substr $read_new_temp, 30+42, -$r;                                                       
                                                        }
                                                        $read = $read_new;
                                                        
                                                        $circle = "yes";
                                                        $noback = "stop";
                                                    }
                                                    elsif ($end_seq =~ m/.*.$start_seq2(.*)$/)
                                                    {
                                                        my $r = length($1);
                                                        my $read_new_temp = $read_new;
                                                        if ($r eq '0')
                                                        {
                                                            $read_new = substr $read_new_temp, 60+42;     
                                                        }
                                                        else
                                                        {
                                                            $read_new = substr $read_new_temp, 60+42, -$r;                                             
                                                        }
                                                        $read = $read_new;
                                                        
                                                        $circle = "yes";
                                                        $noback = "stop";
                                                    }
                                                    elsif ($end_seq =~ m/.*.$start_seq3(.*)$/)
                                                    {
                                                        my $r = length($1);
                                                        my $read_new_temp = $read_new;
                                                        if ($r eq '0')
                                                        {
                                                            $read_new = substr $read_new_temp, 90+42;     
                                                        }
                                                        else
                                                        {
                                                            $read_new = substr $read_new_temp, 90+42, -$r;                                             
                                                        }
                                                        $read = $read_new;
                                                        
                                                        $circle = "yes";
                                                        $noback = "stop";
                                                    }
                                                    elsif ($start_seq =~ m/(.*.)$end_seq1.*$/)
                                                    {
                                                        my $r = length($1);
                                                        my $read_new_temp = $read_new;
                                                        $read_new = substr $read_new_temp, $r, -(42+30);
                                                        $read = $read_new;
                                                        
                                                        $circle = "yes";
                                                        $noback = "stop";
                                                    }
                                                    elsif ($start_seq =~ m/(.*.)$end_seq2.*$/)
                                                    {
                                                        my $r = length($1);
                                                        my $read_new_temp = $read_new;
                                                        $read_new = substr $read_new_temp, $r, -(42+60);
                                                        $read = $read_new;
                                                        
                                                        $circle = "yes";
                                                        $noback = "stop";
                                                    }
                                                    elsif ($start_seq =~ m/(.*.)$end_seq3.*$/)
                                                    {
                                                        my $r = length($1);
                                                        my $read_new_temp = $read_new;
                                                        $read_new = substr $read_new_temp, $r, -(42+90);
                                                        $read = $read_new;
                                                        
                                                        $circle = "yes";
                                                        $noback = "stop";
                                                    }
                                                    if ($merge_read eq "yes" && $circle ne "yes")
                                                    {
                                                        if ($start_seq =~ m/(.*.)$end_seq1_merge.*$/)
                                                        {
                                                            my $r = length($1);
                                                            my $read_new_temp = $read_new;
                                                            $read_new = substr $read_new_temp, $merge_read_length-20, -$r;
                                                            $read = $read_new;
                                                            
                                                            $circle = "yes";
                                                            $noback = "stop";
                                                        }
                                                        elsif ($start_seq =~ m/(.*.)$end_seq2_merge.*$/)
                                                        {
                                                            my $r = length($1);
                                                            my $read_new_temp = $read_new;
                                                            $read_new = substr $read_new_temp, $merge_read_length-100, -$r;
                                                            $read = $read_new;
                                                            
                                                            $circle = "yes";
                                                            $noback = "stop";
                                                        }
                                                        elsif ($start_seq =~ m/(.*.)$end_seq3_merge.*$/)
                                                        {
                                                            my $r = length($1);
                                                            my $read_new_temp = $read_new;
                                                            $read_new = substr $read_new_temp, $merge_read_length-170, -$r;
                                                            $read = $read_new;
                                                            
                                                            $circle = "yes";
                                                            $noback = "stop";
                                                        }
                                                    }
                                                }
                                                if (keys %contigs)
                                                {
                                                    my $total_length = length($read);
                                                    foreach my $contig_tmp (keys %contigs)
                                                    {
                                                        $total_length = $total_length + length($contigs{$contig_tmp});
                                                    }
                                                    if ($total_length > $genome_range_high + $genome_range_high/2)
                                                    {
                                                        $circle = "contigs";
                                                        $noback = "stop";
                                                        $noforward = "stop";
                                                        $contigs{$contig_num."+".$id} = $read;
                                                        $contig_num++;
                                                        delete $seed{$id};
                                                        if ($y > $startprint2)
                                                        {
                                                            print OUTPUT5 "TOTAL_LENGTH_ABOVE_UP_LIMIT\n";
                                                        }
                                                        goto FINISH2;
                                                    }
                                                }
                                                if (length($read) > $genome_range_high)
                                                {
                                                    $circle = "no";
                                                    $noforward = "stop";
                                                    $noback = "stop";

                                                    $contigs{$contig_num."+".$id} = $read;
                                                    $contig_num++;
                                                    delete $seed{$id};
                                                    if ($y > $startprint2)
                                                    {
                                                        print OUTPUT5 "READ_ABOVE_UP_LIMIT\n";
                                                    }
                                                    goto FINISH2;
                                                }
                                chomp $read;
                           
                                if ($position > $insert_size2 - $read_length + 300 || $noback eq "")
                                {
                                    my $read_end_AT = substr $read_short_end2, -$read_length+$right-3;
                                    my $A_rich_test = $read_end_AT =~ tr/A\.//;
                                    my $T_rich_test = $read_end_AT =~ tr/T\.//;
                                    my $G_rich_test = $read_end_AT =~ tr/G\.//;
                                    my $C_rich_test = $read_end_AT =~ tr/C\.//;
                                    if ($A_rich_test > $read_length-$right || $T_rich_test > $read_length-$right || $G_rich_test > $read_length-$right || $C_rich_test > $read_length-$right)
                                    {
                                        $AT_rich = "yes";
                                        goto FINISH;
                                    }
                                    my $read_end_AT2 = substr $read_short_end2, -$overlap-3;
                                    my $A_rich_test2 = $read_end_AT2 =~ tr/A\.//;
                                    my $T_rich_test2 = $read_end_AT2 =~ tr/T\.//;
                                    my $G_rich_test2 = $read_end_AT2 =~ tr/G\.//;
                                    my $C_rich_test2 = $read_end_AT2 =~ tr/C\.//;
                                    if ($A_rich_test2 > $overlap-5 || $T_rich_test2 > $overlap-5 || $G_rich_test2 > $overlap-5 || $C_rich_test2 > $overlap-5)
                                    {
                                        $AT_rich_before = "yes";
                                        $last_chance = ""
                                    }
                                    if ($noback eq "")
                                    {
                                        my $read_start_AT = substr $read_short_start2, 0, $read_length-$right+5;
                                        my $A_rich_test3 = $read_start_AT =~ tr/A/A/;
                                        my $T_rich_test3 = $read_start_AT =~ tr/T/T/;
                                        my $G_rich_test3 = $read_start_AT =~ tr/G/G/;
                                        my $C_rich_test3 = $read_start_AT =~ tr/C/C/;
                                        my $dot_rich_test3 = $read_start_AT =~ tr/\./\./;
                                        if ($A_rich_test3+$dot_rich_test3 > $read_length-$right || $T_rich_test3+$dot_rich_test3 > $read_length-$right || $G_rich_test3+$dot_rich_test3 > $read_length-$right || $C_rich_test3+$dot_rich_test3 > $read_length-$right)
                                        {
                                            $noback = "stop";
                                            $noback{$id} = "stop";
                                            goto MERGE;
                                        }
                                        my $read_start_AT2 = substr $read_short_start2, 0, $overlap+3;
                                        my $A_rich_test32 = $read_start_AT2 =~ tr/A/A/;
                                        my $T_rich_test32 = $read_start_AT2 =~ tr/T/T/;
                                        my $G_rich_test32 = $read_start_AT2 =~ tr/G/G/;
                                        my $C_rich_test32 = $read_start_AT2 =~ tr/C/C/;
                                        my $dot_rich_test32 = $read_start_AT2 =~ tr/\./\./;
                                        if ($A_rich_test32+$dot_rich_test32 > $overlap-5 || $T_rich_test32+$dot_rich_test32 > $overlap-5 || $G_rich_test32+$dot_rich_test32 > $overlap-5 || $C_rich_test32+$dot_rich_test32 > $overlap-5)
                                        {
                                            $AT_rich_before_back = "yes";
                                            print OUTPUT5 $AT_rich_before_back." AT_RICH_BACK\n";
                                            $last_chance_back = "";
                                        }
                                    }
                                                                    
                                    my $s = '0';
                                    my $e = '0';
                                    
                                    while ($s < $read_length-($overlap+$left+1-$containX_short_end2-$containX_short_end2)  && $e < $read_length-($overlap+$left+1-$containX_short_start2-$containX_short_start2))
                                    {
                                        if ($high_coverage eq "yes")
                                        {
                                            my $xd = '0';
                                            foreach (keys %merged_match)
                                            {
                                                $xd++;
                                            }
                                            if ($xd > $coverage_cut_off)
                                            {
                                                $enough = "yes";
                                            }
                                            my $ebrp = '0';
                                            foreach (keys %merged_match_back)
                                            {
                                                $ebrp++;
                                            }
                                            if ($ebrp > $coverage_cut_off)
                                            {
                                                $enough_back = "yes";
                                            }           
                                        }                                    
                                        my $read_end_d = substr $read_short_end2, -($s+$overlap), $overlap;
                                        my $read_start_t = substr $read_short_start2, $e, $overlap;
                                        
                                        if ($containX_short_end2 > 0)
                                        {
                                            my $star = $read_end_d =~ tr/\*//;

                                            $read_end_d = substr $read_short_end2, -($s+$overlap+($star*2)), $overlap+($star*2);
                                            my $star2 = $read_end_d =~ tr/\*//;                                                
                                            while ($star2 > $star)
                                            {
                                                $read_end_d = substr $read_short_end2, -($s+$overlap+($star*2)+(($star2-$star)*2)), $overlap+($star*2)+(($star2-$star)*2);
                                                $star = $star2;
                                                $star2 = $read_end_d =~ tr/\*//;
                                            }   
                                        }
                                        if ($containX_short_start2 > 0)
                                        {
                                            my $star = $read_start_t =~ tr/\*//;
                
                                            $read_start_t = substr $read_short_start2, $e, $overlap+($star*2);
                                            my $star2 = $read_start_t =~ tr/\*//;                                                
                                            while ($star2 > $star)
                                            {
                                                $read_start_t = substr $read_short_start2, $e, $overlap+($star*2)+(($star2-$star)*2);
                                                $star = $star2;
                                                $star2 = $read_start_t =~ tr/\*//;
                                            }
                                        }     
                                        if ($s eq 0)
                                        {                      
                                            my $X = "";
                                            my $Xb = "";
                                            $read_end = $read_end_d;
                                            $read_start = $read_start_t;
                                            $read_end =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                            $read_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                            
                                            my $read_end_c = $read_end;
                                            $read_end_c =~ tr/ATCG/TAGC/;
                                            $read_end_b = reverse($read_end_c);                                                
                                            
                                            my $read_start_c = $read_start;
                                            $read_start_c =~ tr/ATCG/TAGC/;
                                            $read_start_b = reverse($read_start_c);
  
                                            if ($containX_short_end2 > 0)
                                            {
                                                $X = $read_end =~ tr/\*//;
                                                if ($X > 0)
                                                {
                                                    %read_end = build_partial3c $read_end, "";
                                                    %read_end_b = build_partial3c ($read_end_b, "reverse");
                                                }
                                                else
                                                {
                                                    $read_end{$read_end} = undef;
                                                    $read_end_b{$read_end_b} = undef;
                                                }
                                            }
                                            else
                                            {
                                                $read_end{$read_end} = undef;
                                                $read_end_b{$read_end_b} = undef;
                                            }                                                   
                                            if ($containX_short_start2 > 0)
                                            {
                                                $Xb = $read_start =~ tr/\*//;
                                                if ($Xb > 0)
                                                {
                                                    %read_start = build_partial3c $read_start, "reverse";
                                                    %read_start_b = build_partial3c ($read_start_b, "");
                                                }
                                                else
                                                {
                                                    $read_start{$read_start} = undef;
                                                    $read_start_b{$read_start_b} = undef;
                                                }
                                            }
                                            else
                                            {
                                                $read_start{$read_start} = undef;
                                                $read_start_b{$read_start_b} = undef;
                                            }
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $read_end." READ_END\n";
                                                if ($X > 0)
                                                {
                                                   foreach (keys %read_end)
                                                    {
                                                        print OUTPUT5 $_."\n";
                                                    } 
                                                }
                                                print OUTPUT5 $read_start." READ_START\n";
                                                if ($Xb > 0)
                                                {
                                                    foreach (keys %read_start)
                                                    {
                                                        print OUTPUT5 $_."\n";
                                                    } 
                                                }
                                            }
                                        }
                                        if ($enough ne "yes" && $noforward ne "stop") 
                                        {
                                            my $read_end_c = $read_end_d;
                                            
                                            $read_end_c =~ tr/ATCG/TAGC/;
                                            my $read_end_e = reverse($read_end_c);
                                                                                     
   
                                            my %read_end_e;
                                            undef %read_end_e;
                                            
                                            if ($contain_dot_short_end2 > 0 || $containX_short_end2 > 0)
                                            {
                                                my $dot = $read_end_e =~ tr/\.//;
                                                my $X = $read_end_e =~ tr/\*//;
                                                if ($X > 0 || $dot > 0)
                                                {
                                                    %read_end_e = build_partial3c ($read_end_e, "reverse");
                                                }
                                                else
                                                {
                                                    $read_end_e{$read_end_e} = undef;
                                                }
                                            }
                                            else
                                            {
                                                $read_end_e{$read_end_e} = undef;
                                            }
                                            if ($use_regex ne "" || $last_chance eq "yes")
                                            {
                                                my %read_end_e = build_partial3b ($read_end_e, "reverse");
                                                my $X_test = $read_end_e =~ tr/\./\./;
                                                my %list;
                                                undef %list;
                                                if ($X_test < 2 && $use_regex ne "yes2")
                                                {
                                                    %list = build_partial2b %read_end_e;
                                                    %read_end_e = %list;
                                                }

                                                foreach my $list (keys %read_end_e) 
                                                {
                                                    if (exists($hash2c{$list}))
                                                    {                       
                                                        my $search = $hash2c{$list};
                                                                
                                                        $search = substr $search, 1;
                                                        my @search = split /,/,$search;
                                                                        
                                                        foreach my $search (@search)
                                                        {
                                                            my $search_tmp = substr $search, 0, -1;
                                                            my $search_end = substr $search, -1;
                                                            if (exists($hash{$search_tmp}))
                                                            {
                                                                my @search_tmp = split /,/,$hash{$search_tmp};
                                                                my $found;
                                                                if ($search_end eq "1")
                                                                {
                                                                    $found = $search_tmp[0];
                                                                }
                                                                elsif ($search_end eq "2")
                                                                {
                                                                    $found = $search_tmp[1];
                                                                }
                                                                if ($encrypt eq "yes")
                                                                {
                                                                    $found = decrypt $found;
                                                                }
                                                                $merged_match{$search} = $found;
                                                                $merged_match_pos{$search} = $s;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                foreach my $read_end_e2 (keys %read_end_e)
                                                {                                                   
                                                    if (exists($hash2c{$read_end_e2}))
                                                    {                       
                                                        my $search = $hash2c{$read_end_e2};
                                                        $search = substr $search, 1;
                                                        my @search = split /,/,$search;
                                            
                                                        foreach my $search (@search)
                                                        {
                                                            my $search_tmp = substr $search, 0, -1;
                                                            my $search_end = substr $search, -1;
                                                            if (exists($hash{$search_tmp}))
                                                            {
                                                                my @search_tmp = split /,/,$hash{$search_tmp};
                                                                my $found;
                                                                if ($search_end eq "1")
                                                                {
                                                                    $found = $search_tmp[0];
                                                                }
                                                                elsif ($search_end eq "2")
                                                                {
                                                                    $found = $search_tmp[1];
                                                                }
                                                                if ($encrypt eq "yes")
                                                                {
                                                                    $found = decrypt $found;
                                                                }
                                                                $merged_match{$search} = $found;
                                                                $merged_match_pos{$search} = $s;
                                                            }
                                                        }
                                                    }         
                                                }                                                                               
                                            }
                                            if ($last_chance eq "yes")
                                            {
                                                my %read_end_d;
                                                undef %read_end_d;
                                            
                                                if ($contain_dot_short_end2 > 0 || $containX_short_end2 > 0)
                                                {
                                                    my $dot = $read_end_d =~ tr/\.//;
                                                    my $X = $read_end_d =~ tr/\*//;
                                                    if ($X > 0 || $dot > 0)
                                                    {
                                                        %read_end_d = build_partial3b $read_end_d, "";
                                                    }
                                                    else
                                                    {
                                                        $read_end_d{$read_end_d} = undef;
                                                    }
                                                }
                                                else
                                                {
                                                    $read_end_d{$read_end_d} = undef;
                                                }
                                                
                                                my $X_test = $read_end_d =~ tr/\./\./;
                                                my %list;
                                                undef %list;
                                                if ($X_test < 2 && $use_regex ne "yes2")
                                                {
                                                    %list = build_partial2b %read_end_d;
                                                }
                                                else
                                                {
                                                    %list = %read_end_d;
                                                }
                                                foreach my $list (keys %list) 
                                                {
                                                    if (exists($hash2b{$list}))
                                                    {
                                                        my $search = $hash2b{$list};
                                                        
                                                        $search = substr $search, 1;
                                                        my @search = split /,/,$search;
                                                                    
                                                        foreach my $search (@search)
                                                        {
                                                            my $search_tmp = substr $search, 0, -1;
                                                            my $search_end = substr $search, -1;
                                                            if (exists($hash{$search_tmp}))
                                                            {
                                                                my @search_tmp = split /,/,$hash{$search_tmp};
                                                                my $found;
                                                                if ($search_end eq "1")
                                                                {
                                                                    $found = $search_tmp[0];
                                                                }
                                                                elsif ($search_end eq "2")
                                                                {
                                                                    $found = $search_tmp[1];
                                                                }
                                                                if ($encrypt eq "yes")
                                                                {
                                                                    $found = decrypt $found;
                                                                }
                                                                $merged_match{$search} = $found;
                                                                $merged_match_pos{$search} = $s;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if ($noback ne "stop" && $enough_back ne "yes")
                                        {
                                            my %read_start_t;
                                            undef %read_start_t;
                                            
                                            if ($contain_dot_short_start2 > 0 || $containX_short_start2 > 0)
                                            {
                                                my $dot = $read_start_t =~ tr/\.//;
                                                my $X = $read_start_t =~ tr/X|\*//;
                                                if ($X > 0 || $dot > 0)
                                                {
                                                    %read_start_t = build_partial3c ($read_start_t, "reverse");
                                                }
                                                else
                                                {
                                                    $read_start_t{$read_start_t} = undef;
                                                }
                                            }
                                            else
                                            {
                                                $read_start_t{$read_start_t} = undef;
                                            }                                                                          
                                            if ($last_chance_back eq "yes" || $use_regex_back ne "")
                                            {
                                                my %read_start_t = build_partial3b ($read_start_t, "reverse");
                                                my $X_test = $read_start_t =~ tr/\./\./;
                                                my %list;
                                                undef %list;
                                                if ($X_test < 2 && $use_regex_back ne "yes2")
                                                {
                                                    %list = build_partial2b %read_start_t;
                                                    %read_start_t = %list;
                                                }
                                                foreach my $list (keys %read_start_t) 
                                                {
                                                    if (exists($hash2c{$list}))
                                                    {                       
                                                        my $search = $hash2c{$list};
                                                                                
                                                        $search = substr $search, 1;
                                                        my @search = split /,/,$search;
                                                                                        
                                                        foreach my $search (@search)
                                                        {
                                                            my $search_tmp = substr $search, 0, -1;
                                                            my $search_end = substr $search, -1;
                                                            if (exists($hash{$search_tmp}))
                                                            {
                                                                my @search_tmp = split /,/,$hash{$search_tmp};
                                                                my $found;
                                                                if ($search_end eq "1")
                                                                {
                                                                    $found = $search_tmp[0];
                                                                }
                                                                elsif ($search_end eq "2")
                                                                {
                                                                    $found = $search_tmp[1];
                                                                }
                                                                if ($encrypt eq "yes")
                                                                {
                                                                    $found = decrypt $found;
                                                                }
                                                                $merged_match_back{$search} = $found;
                                                                $merged_match_back_pos{$search} = $e;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                foreach my $read_start_t2 (keys %read_start_t)
                                                {
                                                    if (exists($hash2c{$read_start_t2}))
                                                    {                       
                                                        my $search = $hash2c{$read_start_t2};
                                                        $search = substr $search, 1;
                                                        my @search = split /,/,$search;
                                                                                        
                                                        foreach my $search (@search)
                                                        {
                                                            my $search_tmp = substr $search, 0, -1;
                                                            my $search_end = substr $search, -1;
                                                            if (exists($hash{$search_tmp}))
                                                            {
                                                                my @search_tmp = split /,/,$hash{$search_tmp};
                                                                my $found;
                                                                if ($search_end eq "1")
                                                                {
                                                                    $found = $search_tmp[0];
                                                                }
                                                                elsif ($search_end eq "2")
                                                                {
                                                                    $found = $search_tmp[1];
                                                                }
                                                                if ($encrypt eq "yes")
                                                                {
                                                                    $found = decrypt $found;
                                                                }
                                                                $merged_match_back{$search} = $found;
                                                                $merged_match_back_pos{$search} = $e;
                                                            }
                                                        }
                                                    }
                                                }         
                                            }
                                            if ($last_chance_back eq "yes")
                                            {
                                                my $read_start_c = $read_start_t;
                                                $read_start_c =~ tr/ATCG/TAGC/;
                                                my $read_start_e = reverse($read_start_c);
                                                
                                                my %read_start_e;
                                                undef %read_start_e;
                                            
                                                if ($contain_dot_short_start2 > 0 || $containX_short_start2 > 0)
                                                {
                                                    my $dot = $read_start_e =~ tr/\.//;
                                                    my $X = $read_start_e =~ tr/\*//;
                                                    if ($X > 0 || $dot > 0)
                                                    {
                                                        %read_start_e = build_partial3c ($read_start_e, "");
                                                    }
                                                    else
                                                    {
                                                        $read_start_e{$read_start_e} = undef;
                                                    }
                                                }
                                                else
                                                {
                                                    $read_start_e{$read_start_e} = undef;
                                                }
                                                    %read_start_e = build_partial3b ($read_start_e, "");
                                                    my $X_test = $read_start_e =~ tr/\./\./;
                                                    my %list;
                                                    undef %list;
                                                    if ($X_test < 2 && $use_regex_back ne "yes2")
                                                    {
                                                        %list = build_partial2b %read_start_e;
                                                        %read_start_e = %list;
                                                    }
                                                    foreach my $list (keys %read_start_e) 
                                                    {
                                                        if (exists($hash2b{$list}))
                                                        {                       
                                                            my $search = $hash2b{$list};
                                                                                    
                                                            $search = substr $search, 1;
                                                            my @search = split /,/,$search;
                                                                                            
                                                            foreach my $search (@search)
                                                            {
                                                                my $search_tmp = substr $search, 0, -1;
                                                                my $search_end = substr $search, -1;
                                                                if (exists($hash{$search_tmp}))
                                                                {
                                                                    my @search_tmp = split /,/,$hash{$search_tmp};
                                                                    my $found;
                                                                    if ($search_end eq "1")
                                                                    {
                                                                        $found = $search_tmp[0];
                                                                    }
                                                                    elsif ($search_end eq "2")
                                                                    {
                                                                        $found = $search_tmp[1];
                                                                    }
                                                                    if ($encrypt eq "yes")
                                                                    {
                                                                        $found = decrypt $found;
                                                                    }
                                                                    $merged_match_back{$search} = $found;
                                                                }
                                                            }
                                                        }
                                                    }
                                            }
                                        }  
                                        $s++;
                                        $e++;
                                    }
                                }                                     
                                else
                                { 
                                    my $s = '0';
                                    while ($s < $read_length-($overlap+$right))
                                    {                          
                                            my $read_end_d = substr $read, -($s+$overlap), $overlap;
                                            $read_end_d =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
  
                                                if ($s eq 0)
                                                {
                                                    $read_end = substr $read, -$overlap, $overlap;
                                                    $read_end =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;

                                                    my $read_end_c = $read_end;
                                                    $read_end_c =~ tr/ATCG/TAGC/;
                                                    $read_end_b = reverse($read_end_c);
                                                    
                                                    $read_end{$read_end} = undef;
                                                    $read_end_b{$read_end_b} = undef;
                                                    if ($y > $startprint2)
                                                    {
                                                        print OUTPUT5 $read_end." READ_END\n";
                                                        foreach (keys %read_end)
                                                        {
                                                            print OUTPUT5 $_."\n";
                                                        }
                                                    }
                                                }
                                                my $read_end_c = $read_end_d;
                                                $read_end_c =~ tr/ATCG/TAGC/;
                                                my $read_end_e = reverse($read_end_c);
                                                
                                                if ($use_regex eq "yes" || $use_regex eq "yes2")
                                                {
                                                    my %list = build_partial3b $read_end_e;
                                                
                                                    foreach my $list (keys %list)
                                                    {
                                                        if (exists($hash2c{$list}))
                                                        {                       
                                                            my $search = $hash2c{$list};
                                                            $search = substr $search, 1;
                                                            my @search = split /,/,$search;
                                                                
                                                            foreach my $search (@search)
                                                            {
                                                                my $search_tmp = substr $search, 0, -1;
                                                                my $search_end = substr $search, -1;
                                                                if (exists($hash{$search_tmp}))
                                                                {
                                                                    my @search_tmp = split /,/,$hash{$search_tmp};
                                                                    my $found;
                                                                    if ($search_end eq "1")
                                                                    {
                                                                        $found = $search_tmp[0];
                                                                    }
                                                                    elsif ($search_end eq "2")
                                                                    {
                                                                        $found = $search_tmp[1];
                                                                    }
                                                                    if ($encrypt eq "yes")
                                                                    {
                                                                        $found = decrypt $found;
                                                                    }
                                                                    $merged_match{$search} = $found;
                                                                    $merged_match_pos{$search} = $s;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (exists($hash2c{$read_end_e}))
                                                    {                       
                                                        my $search = $hash2c{$read_end_e};
                                                        $search = substr $search, 1;
                                                        my @search = split /,/,$search;
                                                            
                                                        foreach my $search (@search)
                                                        {
                                                            my $search_tmp = substr $search, 0, -1;
                                                            my $search_end = substr $search, -1;
                                                            if (exists($hash{$search_tmp}))
                                                            {
                                                                my @search_tmp = split /,/,$hash{$search_tmp};
                                                                my $found;
                                                                if ($search_end eq "1")
                                                                {
                                                                    $found = $search_tmp[0];
                                                                }
                                                                elsif ($search_end eq "2")
                                                                {
                                                                    $found = $search_tmp[1];
                                                                }
                                                                if ($encrypt eq "yes")
                                                                {
                                                                    $found = decrypt $found;
                                                                }
                                                                $merged_match{$search} = $found;
                                                                $merged_match_pos{$search} = $s;
                                                            }
                                                        }
                                                    } 
                                                }
                                                if ($last_chance eq "yes")
                                                {
                                                    my %read_end_d;
                                                    undef %read_end_d;
                                                
                                                    if ($contain_dot_short_end2 > 0 || $containX_short_end2 > 0)
                                                    {
                                                        my $dot = $read_end_d =~ tr/\.//;
                                                        my $X = $read_end_d =~ tr/\*//;
                                                        if ($X > 0 || $dot > 0)
                                                        {
                                                            %read_end_d = build_partial3b $read_end_d, "";
                                                        }
                                                        else
                                                        {
                                                            $read_end_d{$read_end_d} = undef;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        $read_end_d{$read_end_d} = undef;
                                                    }
                                                    
                                                    my $X_test = $read_end_d =~ tr/\./\./;
                                                    my %list;
                                                    undef %list;
                                                    if ($X_test < 2 && $use_regex ne "yes2")
                                                    {
                                                        %list = build_partial2b %read_end_d;
                                                    }
                                                    else
                                                    {
                                                        %list = %read_end_d;
                                                    }
                                                    foreach my $list (keys %list) 
                                                    {
                                                        if (exists($hash2b{$list}))
                                                        {
                                                            my $search = $hash2b{$list};
                                                            
                                                            $search = substr $search, 1;
                                                            my @search = split /,/,$search;
                                                                        
                                                            foreach my $search (@search)
                                                            {
                                                                my $search_tmp = substr $search, 0, -1;
                                                                my $search_end = substr $search, -1;
                                                                if (exists($hash{$search_tmp}))
                                                                {
                                                                    my @search_tmp = split /,/,$hash{$search_tmp};
                                                                    my $found;
                                                                    if ($search_end eq "1")
                                                                    {
                                                                        $found = $search_tmp[0];
                                                                    }
                                                                    elsif ($search_end eq "2")
                                                                    {
                                                                        $found = $search_tmp[1];
                                                                    }
                                                                    if ($encrypt eq "yes")
                                                                    {
                                                                        $found = decrypt $found;
                                                                    }
                                                                    $merged_match{$search} = $found;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            $s++;
                                        }
                                    }
            my $mmr = '0';
            my $mmbr = '0';
            foreach (keys %merged_match)
            {
                $mmr++;
            }               
            foreach (keys %merged_match_back)
            {
                $mmbr++;
            }
            if ($y > $startprint2)
            {
                print OUTPUT5 $mmr." MATCH_ARRAY_READ\n";
                print OUTPUT5 $mmbr." MATCH_ARRAY_BACK_READ\n";
            }
            
            foreach my $add_read2 (keys %merged_match)
            {
                my $add_read = substr $add_read2, 0, -1;
                $count_reads_all{$add_read} = undef;
            }
            foreach my $add_read2 (keys %merged_match_back)
            {
                my $add_read = substr $add_read2, 0, -1;
                $count_reads_all{$add_read} = undef;
            }
REGEX:
            my $X4 = $read_short_end =~ tr/\*//;
        
            if ($X4 > 0)
            {
                %read_short_end_tmp = build_partial3c ($read_short_end, "");
            }
            else
            {
               $read_short_end_tmp{$read_short_end} = undef;
            }

            my $read_count = '0';
            my $read_ex = '0';
            my $read2_ex = '0';
            
            if ($y > $startprint2)
            {
                if ($extra_regex eq "yes" || $use_regex ne "")
                {
                    print OUTPUT5 "USE_REGEX_REVERSE\n";
                }
            }
            if ($noforward eq "stop")
            {
                goto BACK;
            }
            
NO_MATCH:   foreach my $ln (keys %merged_match)
            {              
                $match = $merged_match{$ln};
                $id_match = $ln;
                chomp $id_match;
                chomp $match;
                
                my $id_match_end1 = substr $id_match, -1, 1;
                my $id_match_tmp1 = substr $id_match, 0, -1;

                if ($id_match_end1 eq "1")
                {
                    $id_pair_match = $id_match_tmp1."2";
                }
                elsif ($id_match_end1 eq "2")
                {
                    $id_pair_match = $id_match_tmp1."1";
                }

                my %match;
                my %match2;
                undef %match;
                undef %match2;
        
                my $n = '0';
                while ($n < length($match)-($overlap+1))
                {
                    my $match_part = substr $match, $n, $overlap;
                    my $match_part2 = substr $match, 0, $n;
                    my $match_part3 = substr $match, $n+$overlap;
                    if (exists($match{$match_part}))
                    {
                        $match{$match_part} = $match_part2;
                    }
                    else
                    {
                        $match{$match_part} = $match_part2;
                        $match2{$match_part} = $match_part3;
                    }
                    $n++;     
                }
                    if ($last_chance eq "yes")
                    {                             
                        my $forward = "";
                        foreach my $line (keys %read_end_b)
                        {
                            my @read_end_b_sub;
                            if ($position > $insert_size+200)
                            {
                                @read_end_b_sub = build_partialb $line;
                            }
                            else
                            {
                                push @read_end_b_sub, $line;
                            }
                                                       
                            my $found_seq = '0';
                            my $match2 = $match;
                                                                
                            foreach my $read_end_b_sub (@read_end_b_sub)
                            {
                                $found_seq = $match2 =~ s/(.)$read_end_b_sub/$1+/;
                                if ($found_seq > 1)
                                {
                                    my $pos = $merged_match_pos{$ln};
                                    my $match4b = substr $match, 0, -$pos;
                                    $match4b =~ s/(.+)$line/$1+/;
                                    my @ext = split /\+/, $match4b;
                                    my $extension5 = $ext[0];
                                    $extension5 =~ tr/ATCG/TAGC/;
                                    $extension = reverse ($extension5);                           
    
                                    $read_count++;
                                    goto LAST1;     
                                }
                                elsif ($found_seq > 0)
                                {
                                    my @ext = split /\+/, $match2;
                                    my $extension5 = $ext[0];
                                    $extension5 =~ tr/ATCG/TAGC/;
                                    $extension = reverse ($extension5);                           
    
                                    $read_count++;
                                    goto LAST1;     
                                }
                            }
                        }

                        my @read_end_sub;
                        if (length($read) > $insert_size+200)
                        {
                            @read_end_sub = build_partialb $read_end;
                        }
                        else
                        {
                            push @read_end_sub, $read_end;
                        }                            
                        my $found_seq = '0';
                        my $match2 = $match;
                                                            
                        foreach my $read_end_sub (@read_end_sub)
                        {
                            $found_seq = $match2 =~ s/.$read_end_sub/+/;
                            if ($found_seq > 0)
                            {
                                my @ext = split /\+/, $match2;
                                $extension = $ext[1];                                                 
                                $read_count++;
                                $forward = "yes";
                                goto LAST1;
                            }
                        }

                        next NO_MATCH;
LAST1:                  my $id_match_end = substr $id_match, -1, 1;
                        my $id_match_tmp = substr $id_match, 0, -1,;
                                   
                        if (index ($id_match_tmp, $id) eq "-1")
                        {
                            my $nuc_exlude = "test";
                            if ($yuyu_option_A eq "A" || $yuyu_option_C eq "C" || $yuyu_option_T eq "T" || $yuyu_option_G eq "G")
                            {
                                $nuc_exlude = substr $extension, 0, 1;
                            }
                            if ($extension ne " " && $extension ne "" && $nuc_exlude ne $yuyu_option_A && $nuc_exlude ne $yuyu_option_C && $nuc_exlude ne $yuyu_option_T && $nuc_exlude ne $yuyu_option_G)
                            {
                                push @matches, $id_match.",".$extension.","."".",".$match.","."";
                                if ($forward eq "yes")
                                {
                                    $extensions2b{$id_match} = $extension;
                                }
                                else
                                {
                                    $extensions1b{$id_match} = $extension;
                                }
                                $extensions2{$extension} = $id_match;
                                push @extensions2, $extension;
                            }
                        }
                        next NO_MATCH;
                    }
                    elsif ($extra_regex eq "yes" || $use_regex eq "yes")
                    {  
                        if ($repetitive_detect eq "")
                        {
                            foreach my $line (keys %read_end_b)
                            {
                                if (exists($match{$line}))
                                {
                                    my $extension5 = $match{$line};
                                    $extension5 =~ tr/ATCG/TAGC/;
                                    $extension = reverse ($extension5);                           
                                                           
                                    $read_count++;
                                    goto FOUND;
                                }
                            }
                        }
                        foreach my $line (keys %read_end_b)
                        {
                            my $found_seq = '0';
                            my $match4 = $match;
                            my @read_end_b_sub = build_partialb $line;
                            
                            foreach my $read_end_b_subc (@read_end_b_sub)
                            {
                                $found_seq = $match4 =~ s/(.)$read_end_b_subc/$1+/g;
                                if ($found_seq > 1)
                                {
                                    my $pos = $merged_match_pos{$ln};
                                    my $match4b = substr $match, 0, -$pos;
                                    $match4b =~ s/(.+)$line/$1+/;
                                    my @ext = split /\+/, $match4b;
                                    my $extension5 = $ext[0];
                                    $extension5 =~ tr/ATCG/TAGC/;
                                    $extension = reverse ($extension5);                           
        
                                    $read_count++;
                                    goto FOUND;     
                                }
                                elsif ($found_seq > 0)
                                {
                                    my @ext = split /\+/, $match4;
                                    my $extension5 = $ext[0];
                                    $extension5 =~ tr/ATCG/TAGC/;
                                    $extension = reverse ($extension5); 
                                                                                           
                                    $read_count++;
                                    goto FOUND;     
                                }
                            }
                        }
                        my @match = split //, $match;
                        
                        foreach my $line (keys %read_end_b)
                        {
                            my @line = split //,$line;
                            my $gh = length($match)-($overlap);
                            my $th = '0';
         
REGEXORNOT:                 while ($gh > 1)
                            {
                                my $d = '0';
                                my $next = "";
                                                                    
                                while ($d < $overlap)
                                {
                                    $th = $d + $gh;
                                    if ($match[$th] eq $line[$d])
                                    {
                                    }
                                    elsif ($line[$d] eq ".")
                                    {
                                    }
                                    elsif ($next eq "")
                                    {
                                        $next = "yes";
                                    }
                                    elsif ($next eq "yes")
                                    {
                                        $next = "yes2";
                                    }
                                    elsif ($next eq "yes2")
                                    {
                                        $next = "yes3";
                                    }
                                    else
                                    {
                                        $gh--;
                                        goto REGEXORNOT;
                                    }
                                    $d++    
                                }
                                my $extension5 = substr $match,0, $gh;
                                $extension5 =~ tr/ATCG/TAGC/;
                                $extension = reverse ($extension5);                           
                                                                                       
                                $read_count++;
                                goto FOUND;     
                            }
                        }
                        next NO_MATCH;
                    }
                    else
                    {
                        if ($repetitive_detect eq "")
                        {
                            foreach my $line (keys %read_end_b)
                            {
                                if (exists($match{$line}))
                                {
                                    my $extension5 = $match{$line};
                                    $extension5 =~ tr/ATCG/TAGC/;
                                    $extension = reverse ($extension5);                           
                                                           
                                    $read_count++;
                                    goto FOUND;
                                }
                            }
                        }
                        foreach my $line (keys %read_end_b)
                        {
                            my $found_seq = '0';
                            my $match4 = $match;
                            
                            $found_seq = $match4 =~ s/(.)$line/$1+/g;
                            if ($found_seq > 1)
                            {
                                my $pos = $merged_match_pos{$ln};
                                my $match4b = substr $match, 0, -$pos;
                                $match4b =~ s/(.+)$line/$1+/;
                                my @ext = split /\+/, $match4b;
                                my $extension5 = $ext[0];
                                $extension5 =~ tr/ATCG/TAGC/;
                                $extension = reverse ($extension5);                           

                                $read_count++;
                                goto FOUND;     
                            }
                            elsif ($found_seq > 0)
                            {
                                my @ext = split /\+/, $match4;
                                my $extension5 = $ext[0];
                                $extension5 =~ tr/ATCG/TAGC/;
                                $extension = reverse ($extension5);                           

                                $read_count++;
                                goto FOUND;     
                            }
                        }
                        next NO_MATCH;
                    }

FOUND:              if ($last_chance eq "yes")
                    {           
                        next NO_MATCH;
                    }         
                    if ($extension ne "NOOO")
                    {
                        my $id_match_b = $id_match;
                        my $id_match_end = substr $id_match_b, -1, 1,"",;

                        if (exists($hash{$id_match_b}))
                        {
                            my @id_match_b = split /,/, $hash{$id_match_b};
                            
                            if ($id_match_end eq "1")
                            {
                                $match_pair = $id_match_b[1];
                            }
                            elsif ($id_match_end eq "2")
                            {
                                $match_pair = $id_match_b[0];
                            }
                            else
                            {
                                next NO_MATCH;
                            }                                
                            chomp($match_pair);
                            if ($encrypt eq "yes")
                            {
                                $match_pair = decrypt $match_pair;
                            }
                  
                            $match_end = substr $match_pair, -($overlap+$right), $overlap;
                            $match_start = substr $match_pair, $left, $overlap;
                            
                            my $is = $overlap;
                            if ($indel_split ne '0')
                            {
                                $is = length($match_pair2)-40;
                                $is = $overlap+$indel_split;
                            }                 
                                   
                            if ($extra_regex eq "yes")
                            {                                                   
                                my $match_pair_middle = substr $match_pair, 8, $is+10;
                                $match_pair_middle =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                my @match_pair_middle_sub = split //, $match_pair_middle;

                                my $read_short_end_tempie = substr $read, -($insert_size*$insert_range)+length($extension), ((($insert_size*$insert_range)-$insert_size)*2)+$is+10+8;
                                $read_short_end_tempie =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                undef %read_short_end_tmp;
                                %read_short_end_tmp = build_partial3c $read_short_end_tempie;
                                                                     
                                foreach my $line (keys %read_short_end_tmp)
                                {
                                    my @line = split //,$line;
                                    my $gh = '0';
                                    my $th = '0';
                                    $line =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                 
CHECK_PAIR:                         while ($gh < length($line)-length($match_pair_middle))
                                    {
                                        my $d = '0';
                                        my $next = "";
                                                                            
                                        while ($d < length($match_pair_middle))
                                        {
                                            $th = $d + $gh;
                                            if ($match_pair_middle_sub[$d] eq $line[$th])
                                            {
                                            }
                                            elsif ($line[$th] eq ".")
                                            {
                                            }
                                            elsif ($match_pair_middle_sub[$d] eq ".")
                                            {
                                            }
                                            elsif ($next eq "")
                                            {
                                                $next = "yes";
                                            }
                                            elsif ($next eq "yes")
                                            {
                                                $next = "yes2";
                                            }
                                            elsif ($next eq "yes2")
                                            {
                                                $next = "yes3";
                                            }
                                            else
                                            {
                                                $gh++;
                                                goto CHECK_PAIR;
                                            }
                                            $d++    
                                        }
                                        
                                        $counttest1++;
                                        if ($insert_size_correct eq "yes")
                                        {
                                            my $line_tmp_length = length($line) - length($match_pair_middle);
                                            my $cal = -($insert_size*$insert_range) + (($read_length-$overlap-8)/2) + 15+($overlap/2) - (($overlap-38)/2) + ($insert_size*(($insert_range-1)*2)) + ($overlap-38);
                                            if ($cal > 0)
                                            {
                                                $cal = '0';
                                            }
                                            my $insert_size_tmp = length($extension) - $cal + 8 + length($match_pair_middle) + $line_tmp_length;

                                            push @insert_size, $insert_size_tmp;
                                        }
                                        $extension_match = "";
                                        goto SKIP3;    
                                    }
                                }   
                                $extension_match = "NOOO";
                            }
                            else
                            {
                                my $match_pair_middle = substr $match_pair, 8, $is+10;

                                my $read_short_end_tempie = substr $read, -($insert_size*$insert_range)+length($extension), ((($insert_size*$insert_range)-$insert_size)*2)+$is+10+8;
                                $read_short_end_tempie =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                undef %read_short_end_tmp;
                                my $test_dot = $read_short_end_tempie =~ tr/\./\./;
                                
                                if ($test_dot > 1)
                                {
                                    %read_short_end_tmp = build_partial3c $read_short_end_tempie;
                                    my $ff = '0';
                                    my %hash_read_short_end;
                                    undef %hash_read_short_end;
                                    foreach my $read_short_end_tempie (keys %read_short_end_tmp)
                                    {
                                        $ff = '0';
                                        while ($ff < (length($read_short_end_tempie)-length($match_pair_middle)-1))
                                        {
                                            my $read_short_end_part = substr $read_short_end_tempie, $ff, length($match_pair_middle);
                                                                           
                                            $hash_read_short_end{$read_short_end_part} = undef;
                                            $ff++;
                                        }
                                    }
                                    foreach my $line (keys %hash_read_short_end)
                                    {
                                        my $found_seq = '0';
                                        my $star_first = substr $line, 0, 1;
                                        if ($star_first eq "*")
                                        {
                                            $line = substr $line, 1;
                                        }  
    
                                        $found_seq = $match_pair_middle =~ s/$line/$line/;
                                        if ($found_seq > 0)
                                        {
                                            if ($insert_size_correct eq "yes")
                                            {
                                                my $line_tmp = $line;
                                                $line_tmp =~ s/^.*$match_pair_middle//;
                                                my $cal = -($insert_size*$insert_range) + (($read_length-$overlap-8)/2) + 15+($overlap/2) - (($overlap-38)/2) + ($insert_size*(($insert_range-1)*2)) + ($overlap-38);
                                                if ($cal > 0)
                                                {
                                                    $cal = '0';
                                                }
                                                my $insert_size_tmp = length($extension) - $cal + 8 + length($match_pair_middle) + length($line_tmp);
    
                                                push @insert_size, $insert_size_tmp;
                                            }
                                            $extension_match = "";
                                            goto SKIP3;
                                        }
                                    }
                                    $extension_match = "NOOO";
                                }
                                else
                                {
                                    %read_short_end_tmp = build_partial3b $read_short_end_tempie;
                                    foreach my $line (keys %read_short_end_tmp)
                                    {
                                        my $found_seq = '0';
    
                                        $found_seq = $line =~ s/$match_pair_middle/$match_pair_middle/;
                                        if ($found_seq > 0)
                                        {
                                            if ($insert_size_correct eq "yes")
                                            {
                                                my $line_tmp = $line;
                                                $line_tmp =~ s/^.*$match_pair_middle//;
                                                my $cal = -($insert_size*$insert_range) + (($read_length-$overlap-8)/2) + 15+($overlap/2) - (($overlap-38)/2) + ($insert_size*(($insert_range-1)*2)) + ($overlap-38);
                                                if ($cal > 0)
                                                {
                                                    $cal = '0';
                                                }
                                                my $insert_size_tmp = length($extension) - $cal + 8 + length($match_pair_middle) + length($line_tmp);
    
                                                push @insert_size, $insert_size_tmp;
                                            }
                                            $extension_match = "";
                                            goto SKIP3;
                                        }
                                    }
                                    $extension_match = "NOOO";
                                }                              
                            }
SKIP3:
                            my $nuc_exlude = "test";
                            if ($yuyu_option_A eq "A" || $yuyu_option_C eq "C" || $yuyu_option_T eq "T" || $yuyu_option_G eq "G")
                            {
                                $nuc_exlude = substr $extension, 0, 1;
                            }
                            if ($extension_match ne "NOOO" && ($extension ne " " && $extension ne "") && $nuc_exlude ne $yuyu_option_A && $nuc_exlude ne $yuyu_option_C && $nuc_exlude ne $yuyu_option_T && $nuc_exlude ne $yuyu_option_G)
                            {
                                $read_ex++;
                                push @matches, $id_match.",".$extension.",".$id_pair_match.",".$match.",".$match_pair;
                                    
                                if ($extension ne " " && $extension ne "")
                                {
                                    $extensions1{$extension} = $id_match;
                                    $extensions1b{$id_match} = $extension;
                                    push @extensions1, $extension;
                                }                                                            
                            }                                                                                           
                    }                                  
                }
            }

            %extensions = (%extensions1, %extensions2);
            %extensions_original = %extensions;
            %extensionsb = (%extensions1b, %extensions2b);
            @extensions = (@extensions1, @extensions2);

            my $ext = '0';
            my $ext_total = '0';
            foreach (@extensions)
            {
                $ext++;
            }
            $ext_total = $ext;
                        
            if ($y > $startprint2)
            {
                print OUTPUT5 "\n".$read_count ." READ_COUNT\n";
                print OUTPUT5 $read_ex ." READ_EX\n";
                print OUTPUT5 $ext ." EXTENSIONS\n";
            }
            
            
            if ($ext < 8 && $extra_regex eq "" && $last_chance ne "yes" && $indel_split eq '0' && $ext < $mmr && $repetitive_detect eq "")
            {
                undef %extensions1;
                undef %extensions2;
                undef %extensions;
                undef %extensions_original;
                undef %extensions1b;
                undef %extensions2b;
                undef %extensionsb;
                undef @extensions1;
                undef @extensions2;
                undef @extensions;
                undef @matches;

                $extra_regex = "yes";
                goto REGEX;
            }          
            $extra_regex = "";
            if ($y > $startprint && $print_log eq '2')
            {
                foreach my $matches (@matches)
                {
                    my @matchesb;
                    undef @matchesb;
                    @matchesb = split /,/, $matches;
                    my $m_reverse = reverse($matchesb[3]);
                    $m_reverse =~ tr/ACTG/TGAC/;
                    my $mp_reverse = reverse($matchesb[4]);
                    $mp_reverse =~ tr/ACTG/TGAC/;
                    print OUTPUT5 $matchesb[0].",".$matchesb[1]."\n";
                    print OUTPUT5 $m_reverse." MATCH_REVERSE\n";
                }               
            }
            
            my $id_original = $id;
            
            my @extensions_group1;
            my @extensions_group2;
            my %extensions_group1;
            my %extensions_group2;
            undef @extensions_group1;
            undef @extensions_group2;
            undef %extensions_group1;
            undef %extensions_group2;
            my @extensions_group3;
            my @extensions_group4;
            my %extensions_group3;
            my %extensions_group4;
            undef @extensions_group3;
            undef @extensions_group4;
            undef %extensions_group3;
            undef %extensions_group4;
            
            my $no_SNR = "";
            my %extensions_backup = %extensions;
            my @extensions_backup = @extensions;
            
SPLIT:
            if ($split eq "yes")
            {
                @extensions = @extensions_group2;
                %extensions = %extensions_group2;
                $split = "yes2";
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            elsif ($split eq "yes2")
            {
                @extensions = @extensions_group3;
                %extensions = %extensions_group3;
                $split = "yes3";
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            elsif ($split eq "yes3")
            {
                @extensions = @extensions_group4;
                %extensions = %extensions_group4;
                $split = "yes4";
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            elsif ($split eq "yes4")
            {
                @extensions = @extensions_group1;
                %extensions = %extensions_group1;
                $split = "yes5";
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            $id = $id_original;
            $position = $position{$id};       
            my $l = '0';
            my $best_extension = "";
            my $SNP = "";
            my $A_SNP = '0';
            my $C_SNP = '0';
            my $T_SNP = '0';
            my $G_SNP = '0';
            my $position_SNP = $position;
            my $pos_SNP = '0';
            
            my $A_SNP2 = '0';
            my $C_SNP2 = '0';
            my $T_SNP2 = '0';
            my $G_SNP2 = '0';
            my $position_SNP2 = $position;
            my $pos_SNP2 = '0';
            
            my $A_SNP3 = '0';
            my $C_SNP3 = '0';
            my $T_SNP3 = '0';
            my $G_SNP3 = '0';
            my $position_SNP3 = $position;
            my $pos_SNP3 = '0';
            
            my %SNR_count;
            my %extensions_new;
            my @extensions_new;
            undef %SNR_count;
            undef %extensions_new;
            undef @extensions_new;
            my $SNR_test = "";
            my $most_SNR = '0';
            my $most_SNR2 = '0';
                      
            if ($SNR_read ne "" && $split eq "" && $SNR_read2 ne "")
            {
                $SNR_test = "yes2";
                if ($SNR_read eq "yes")
                {
                    $SNR_test = "yes2";
                    my $G = '0';
                    my $G2 = '0';
                    my $no_SNR1 = "";
                    my $second_round;
SNR1:               foreach my $extensions (@extensions)
                    { 
                        my @chars = split//, $extensions;
                        my $e = '0';
                        my $check_only_SNR = $extensions =~ s/$SNR_nucleo/$SNR_nucleo/g;                       
                        if (length($extensions) eq $check_only_SNR)
                        {
                            $G2++;
                        }
                        $G++;
                        if ($second_round eq "")
                        {
                            $no_SNR1 = "yes";
                        }
                        while ($SNR_nucleo eq $chars[$e] || $no_SNR1 eq "")
                        {
                            if ($SNR_nucleo ne $chars[$e])
                            {
                                $no_SNR1 = "yes";
                            }
                            $e++;
                        }
                        if ($e < length($extensions))
                        {                      
                            $SNR_count{$extensions} = $e;
                            $SNR_length{$e} .= exists $SNR_length{$e} ? ",$extensions" : $extensions;
                        }
                        $no_SNR1 = "";
                    }
                    my $SNR_length_count2 = '0';
                    my $SNR_length_reads = "";
                    my $first = "";
                    foreach my $SNR_length (keys %SNR_length)
                    {
                        my $SNR_length_count = $SNR_length{$SNR_length} =~ tr/,/,/;
                        if ($SNR_length_count > $SNR_length_count2)
                        {
                            $SNR_length_count2 = $SNR_length_count;
                            $SNR_length_reads = $SNR_length{$SNR_length};
                            $most_SNR = $SNR_length;
                            $first = $SNR_length_count;
                        }
                    }
                    my @SNR_length = split/,/, $SNR_length_reads;
                    my $repetitive_test = substr $read_short_end2, -10, 10;
                    my $SNR_checkSNR = $repetitive_test =~ s/$SNR_nucleo/$SNR_nucleo/g;
                    
                    $SNR_length_count2 = '0';
                    foreach my $SNR_length (keys %SNR_length)
                    {
                        my $SNR_length_count = $SNR_length{$SNR_length} =~ tr/,/,/;
                        if ($SNR_length_count > $SNR_length_count2 && $SNR_length_count ne $first)
                        {
                            $SNR_length_count2 = $SNR_length_count;
                            $most_SNR2 = $SNR_length;
                        }
                    }
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $first." 1 ".$SNR_length_count2." 2 ".$G." G\n";
                    }
                    if ($first < 0.8*$G && $second_round eq "" && $first > 0.35*$G && $SNR_length_count2 > 0.35*($G-$G2))
                    {                
                        $no_SNR = "yes";
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "NO_SNR\n";
                        }
                        goto NUCLEO0;
                    }
                    elsif ($first < 0.8*$G && $second_round eq "")
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 $most_SNR." MOST_SNR SECOND ROUND\n";
                        }
                        $second_round = "yes";
                        $no_SNR1 = "";
                        $most_SNR = '0';
                        $G = '0';
                        undef %SNR_count;
                        undef %SNR_length;
                        goto SNR1;
                    }
                    else
                    {
                        foreach my $SNRie (@SNR_length)
                        {
                            if (exists($extensions{$SNRie}))
                            {
                                $extensions_new{$SNRie} = $extensions{$SNRie};
                                push @extensions_new, $SNRie;
                            }
                        }
                        %extensions = %extensions_new;
                        @extensions = @extensions_new;
                    }            
                }
                if ($SNR{$id} eq "yes2_double")
                {
                    $SNR_test = "yes2_double";
                    foreach my $extensions (@extensions)
                    { 
                        my @chars = split//, $extensions;
                        my $e = '0';
                        if ($SNR_nucleo eq $chars[$e].$chars[$e+1] )
                        {
                            while ($SNR_nucleo eq $chars[$e].$chars[$e+1])
                            {
                                my $tempie = reverse $extensions;
                                chop $tempie;
                                chop $tempie;
                                $extensions = reverse $tempie;
                                $e++;
                                $e++;
                            }
                        }
                        else
                        {  
                            while ($SNR_nucleo eq $chars[$e+1].$chars[$e])
                            {
                                my $tempie = reverse $extensions;
                                chop $tempie;
                                chop $tempie;
                                $extensions = reverse $tempie;
                                $e++;
                                $e++;
                            }
                        }
                        $extensions_new{$extensions} = $extensions{$extensions};
                        push @extensions_new, $extensions;
                        if ($e < length($extensions))
                        {                      
                            $SNR_count{$extensions} = $e;
                        }
                    }
                    %extensions = %extensions_new;
                    @extensions = @extensions_new;
                }
                delete $SNR{$id};
            }
NUCLEO0:            
            if ($SNR_read ne "")
            {
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            my $extra_l = '0';

NUCLEO:     while ($l < $read_length - ($overlap+$left-1) + $extra_l)
            {
                my $A = '0';
                my $C = '0';
                my $T = '0';
                my $G = '0';
                
                if ($SNR_read ne "" && $l > 0)
                {
                    my $best_extension_tmp8 = $best_extension;
                    $best_extension_tmp8 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    
                    my $last_nuc = substr $best_extension_tmp8, -1;
                    my $arrSize1 = @extensions;
                    if ($last_nuc ne '.' && $arrSize1 > 4)
                    {
                        my @extensions_tmp;
                        undef @extensions_tmp;
                        foreach my $extensions (@extensions)
                        {
                            my @chars = split//, $extensions;
                            if ($chars[$l-1] eq $last_nuc || length($extensions) < $l)
                            {
                                push @extensions_tmp, $extensions;
                            }
                        }
                        my $arrSize2 = @extensions_tmp;

                        if ($arrSize1 ne $arrSize2)
                        {
                            undef @extensions;
                            @extensions = @extensions_tmp;
                                                       
                            my $best_extension_dot = $best_extension_tmp8 =~ tr/\./\./;
                            if ($best_extension_dot > 0)
                            {
                                $l = 0;
                                $SNP = "";
                                $best_extension = "";
                                goto NUCLEO;
                            }
                        }
                    }               
                }
                
                foreach my $extensions (@extensions)
                {                                
                    my @chars = split//, $extensions;
                    
                    if ($chars[$l] eq "A")
                    {
                        $A++;
                    }
                    elsif ($chars[$l] eq "C")
                    {
                        $C++;
                    }
                    elsif ($chars[$l] eq "T")
                    {
                        $T++;
                    }
                    elsif ($chars[$l] eq "G")
                    {
                        $G++;
                    }
                }
                my $c = '2.8';
                my $q = '3';
                
                if ($ext > 22 && $SNR_read eq "")
                {
                    $c = '3.7';
                }
                if ($ext > 6 && $SNR_read eq "" && ($type ne "chloro" || $extensions_before ne ""))
                {
                    $c = '5';
                }
                if ($ext > 22 && $SNR_read eq "" && $type ne "chloro")
                {
                    $c = '6.5';
                }
                if ($ext > 38 && $SNR_read eq "" && $type ne "chloro")
                {
                    $c = '8.4';
                }
                if ($repetitive_detect ne "" && $ext < 23 && $SNR_read eq "")
                {
                   $c = '7';
                }
                if ($repetitive_detect ne "" && $ext >= 23 && $SNR_read eq "")
                {
                   $c = '9';
                }
                if ($repetitive_detect2 eq "yes" || ($repetitive_detect ne "" && $ext > 150))
                {
                   $c = '15';
                }
                if ($extensions_before eq "yes")
                {
                   $c = '6.3';
                }
                if ($type eq "mito_plant")
                {
                    $c += 2;
                }
                my $v = '6';
                my $s = '3';
                my $z = '1';
                if ($split ne "")
                {
                    $v = '10';
                    $z = '0';
                    $s = '2';
                }

                if ($A > ($C + $T + $G)*$c && (($A > $s && ($ext)/($A + $T + $G + $C) < $q) || ($A > $z && $l < $v && ($C + $T + $G) eq 0)))
                {
                    $best_extension = $best_extension."A";
                }
                elsif ($C > ($A + $T + $G)*$c && (($C > $s && ($ext)/($A + $T + $G + $C) < $q) || ($C > $z && $l < $v && ($A + $T + $G) eq 0)))
                {
                    $best_extension = $best_extension."C"; 
                }
                elsif ($T > ($A + $C + $G)*$c && (($T > $s && ($ext)/($A + $T + $G + $C) < $q) || ($T > $z && $l < $v && ($C + $A + $G) eq 0)))
                {
                    $best_extension = $best_extension."T";  
                }
                elsif ($G > ($C + $T + $A)*$c && (($G > $s && ($ext)/($A + $T + $G + $C) < $q) || ($G > $z && $l < $v && ($C + $T + $A) eq 0)))
                {
                    $best_extension = $best_extension."G";
                }
                elsif (($SNP_active eq "yes" || ($SNR_read ne "" && $l > 0) || ($extensions_before eq "yes" && $ext_before ne "yes")) && $SNP eq "" && ($A + $T + $G + $C) > 4 && $l < 15 && ($ext)/($A + $T + $G + $C) < $q && $split eq "") 
                {
                    delete $SNP_active{$id};
                    $SNP = "yes";
                    $A_SNP = $A;
                    $C_SNP = $C;
                    $T_SNP = $T;
                    $G_SNP = $G;
                    $position_SNP += $l;
                    $pos_SNP = $l;
                    
                    my $IUPAC = IUPAC($A,$C,$T,$G);
                    $best_extension = $best_extension.$IUPAC;
                }
                elsif ($SNP eq "yes" && ($A + $T + $G + $C) > 4  && $l < 15)
                {
                    $SNP = "yes2";
                    $A_SNP2 = $A;
                    $C_SNP2 = $C;
                    $T_SNP2 = $T;
                    $G_SNP2 = $G;
                    $position_SNP2 += $l;
                    $pos_SNP2 = $l;
                    
                    my $IUPAC = IUPAC($A,$C,$T,$G);
                    $best_extension = $best_extension.$IUPAC;
                }
                elsif ($SNP eq "yes2" && ($A + $T + $G + $C) > 4  && $l < 15)
                {
                    $SNP = "yes3";
                    $A_SNP3 = $A;
                    $C_SNP3 = $C;
                    $T_SNP3 = $T;
                    $G_SNP3 = $G;
                    $position_SNP3 += $l;
                    $pos_SNP3 = $l;
                    
                    my $IUPAC = IUPAC($A,$C,$T,$G);
                    $best_extension = $best_extension.$IUPAC;
                }
                elsif ($SNP eq "yes3" && ($pos_SNP ne 0 || ($pos_SNP3 > $pos_SNP+12 && $l > 11 ) || ($extensions_before eq "yes" && $l > 9)) && $split eq "")
                {
                    $SNP = "yes4";
                    my $g = $l;
                    my $pos_SNP_tmp = $pos_SNP;
                    if ($pos_SNP3 > $pos_SNP+12)
                    {
                        $pos_SNP_tmp = $pos_SNP3;
                    }
                    if ($extensions_before ne "yes")
                    {
                        while ($g > $pos_SNP_tmp)
                        {                                         
                            chop($best_extension);
                            $g--;
                        }
                    }
                    last  NUCLEO;
                }
                

                elsif (($SNP eq "yes3" && $pos_SNP eq 0 && $l <= 15) || ($indel_split_skip ne "yes" && $l eq 0 && $ext > 4))
                {
                    print OUTPUT5 $SNP." SNP\n";
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $SNP." SNP\n";
                    }
                    if ($SNR_test ne "" && $no_SNR ne "yes")
                    {
                        $l = '0';
                        $best_extension = "";
                        $no_SNR = "yes";
                        %extensions = %extensions_backup;
                        @extensions = @extensions_backup;
                        
                        goto SPLIT;
                    }
                    if ($SNP eq "")
                    {
                        $A_SNP = $A;
                        $C_SNP = $C;
                        $T_SNP = $T;
                        $G_SNP = $G;
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 $best_extension." BEST_EXTENSIONll\n";
                            print OUTPUT5 $A_SNP." A\n";
                            print OUTPUT5 $C_SNP." C\n";
                            print OUTPUT5 $T_SNP." T\n";
                            print OUTPUT5 $G_SNP." G\n";
                        }
                    }
                    $best_extension = "";
                    $split = "yes";
                    undef @firstSNP_max;

                    if ($A_SNP >= ($C_SNP+$A_SNP+$T_SNP+$G_SNP)*0.035 && $A_SNP > 0)
                    {
                        if (exists($yuyu_option{$id.'A'}))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "YUYU_A_EXISTS\n";
                            }
                        }
                        else
                        {
                            push @firstSNP_max, "A";
                        }  
                    }
                    if ($C_SNP >= ($C_SNP+$A_SNP+$T_SNP+$G_SNP)*0.035 && $C_SNP > 0)
                    {
                        if (exists($yuyu_option{$id.'C'}))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "YUYU_C_EXISTS\n";
                            }
                        }
                        else
                        {
                            push @firstSNP_max, "C";
                        }  
                    }
                    if ($T_SNP >= ($C_SNP+$A_SNP+$T_SNP+$G_SNP)*0.035 && $T_SNP > 0)
                    {
                        if (exists($yuyu_option{$id.'T'}))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "YUYU_T_EXISTS\n";
                            }
                        }
                        else
                        {
                            push @firstSNP_max, "T";
                        }  
                    }
                    if ($G_SNP >= ($C_SNP+$A_SNP+$T_SNP+$G_SNP)*0.035 && $G_SNP > 0)
                    {
                        if (exists($yuyu_option{$id.'G'}))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "YUYU_G_EXISTS\n";
                            }
                        }
                        else
                        {
                            push @firstSNP_max, "G";
                        }  
                    }
                    
                    $count_split = @firstSNP_max;
                    my $count_split_tmp = '0';
                    
                    if ($count_split eq '2')
                    {
                        $delete_third = "yes";
                        $delete_second = "yes";
                    }
                    if ($count_split eq '3')
                    {
                        $delete_third = "yes";
                    }
                    
                    foreach my $firstSNP_max (@firstSNP_max)
                    {                                            
                        foreach my $extensions_tmp (@extensions)
                        {                       
                            my @chars = split//, $extensions_tmp;
                            if ($chars[0] eq $firstSNP_max && $count_split_tmp eq '0')
                            {
                                $extensions_group1{$extensions_tmp} = $extensions{$extensions_tmp};
                                push @extensions_group1, $extensions_tmp;
                            }
                            elsif ($chars[0] eq $firstSNP_max && $count_split_tmp eq '1')
                            {
                                $extensions_group2{$extensions_tmp} = $extensions{$extensions_tmp};
                                push @extensions_group2, $extensions_tmp;
                            }
                            elsif ($chars[0] eq $firstSNP_max && $count_split_tmp eq '2')
                            {
                                $extensions_group3{$extensions_tmp} = $extensions{$extensions_tmp};
                                push @extensions_group3, $extensions_tmp;
                            }
                            elsif ($chars[0] eq $firstSNP_max && $count_split_tmp eq '3')
                            {
                                $extensions_group4{$extensions_tmp} = $extensions{$extensions_tmp};
                                push @extensions_group4, $extensions_tmp;
                            }
                        }
                        $count_split_tmp++;
                    }
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $count_split." COUNT_SPLIT\n";
                    }
                    if ($count_split eq '1')
                    {
                        $delete_third = "yes";
                        $delete_second = "yes";
                    }
                    goto SPLIT;
                }
                else
                { 
                    last  NUCLEO;
                }
                $l++;
            }

            if($split eq "" && $rep_detect2 eq "yes")
            {
                
            }
            if ($SNP eq "yes3" && $split eq "" && length($best_extension) < 15)
            {
                $SNP = "yes4";
                my $g = $l;
                my $pos_SNP_tmp = $pos_SNP;
                    
                while ($g > $pos_SNP_tmp)
                {                                         
                    chop($best_extension);
                    $g--;
                }
                if ($y > $startprint2)
                {
                    print OUTPUT5 $best_extension." BEST_EXTENSION_chopped\n";
                }
            }
            my $last_nucleo = substr $best_extension, -1;
            $last_nucleo =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
            while ($last_nucleo eq '.')
            {
                chop($best_extension);
                $last_nucleo = substr $best_extension, -1;
                $last_nucleo =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
            }
            if ($split eq "yes2" || ($split eq "yes3" && $count_split > 2) || ($split eq "yes4" && $count_split > 3))
            {               
                if ($split eq "yes2")
                {
                    $best_extension2 = $best_extension;
                }
                elsif ($split eq "yes3")
                {
                    $best_extension3 = $best_extension;
                }
                elsif ($split eq "yes4")
                {
                    $best_extension4 = $best_extension;
                }
                my $best_extension_tmp = $best_extension;
                $best_extension_tmp =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;

                my $contig_id2_prev = $id;
                
                if ($split eq "yes2")
                {
                    $contig_id2 = "2".length($read);
                    $contig_id2{$contig_id2_prev} = $contig_id2;
                    $contig_id2_prev = $contig_id2;
                    $contig_read2 = substr $read, -$read_length;
                    $contig_read2 = $contig_read2.$best_extension_tmp;
                }
                elsif ($split eq "yes3")
                {
                    $contig_id3 = "3".length($read);
                    $contig_id3{$contig_id2_prev} = $contig_id3;
                    $contig_id2_prev = $contig_id3;
                    $contig_read2 = substr $read_short_end2, -$read_length;
                    $contig_read2 = $contig_read2.$best_extension_tmp;
                }
                elsif ($split eq "yes4")
                {
                    $contig_id4 = "4".length($read);
                    $contig_id4{$contig_id2_prev} = $contig_id4;
                    $contig_id2_prev = $contig_id4;
                    $contig_read2 = substr $read_short_end2, -$read_length;
                    $contig_read2 = $contig_read2.$best_extension_tmp;
                }                
                
                if ($y > $startprint2)
                {
                    if ($split eq "yes2")
                    {
                        print OUTPUT5 "GROUP2\n";
                        foreach my $extensions_tmp (@extensions_group2)
                        {  
                            print OUTPUT5 $extensions_tmp."\n";                        
                        }
                        print OUTPUT5 $best_extension_tmp." BEST_EXTENSION2\n\n";
                    }
                    elsif ($split eq "yes3")
                    {
                         print OUTPUT5 "GROUP3\n";
                        foreach my $extensions_tmp (@extensions_group3)
                        {  
                            print OUTPUT5 $extensions_tmp."\n";                        
                        }
                        print OUTPUT5 $best_extension_tmp." BEST_EXTENSION3\n\n";
                    }
                    elsif ($split eq "yes4")
                    {
                         print OUTPUT5 "GROUP4\n";
                        foreach my $extensions_tmp (@extensions_group4)
                        {  
                            print OUTPUT5 $extensions_tmp."\n";                        
                        }
                        print OUTPUT5 $best_extension_tmp." BEST_EXTENSION4\n\n";
                    }    
                }
                if ($before eq "yes")
                {
                    if ((length($best_extension2) < 3 || (length($best_extension2) < 6 && ($ext > 15 || $SNR_read ne ""))) && $repetitive_detect eq "" && $before eq "yesss")
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nDELETE BEST EXTENSION2/3/4\n\n";
                        }
                        if ($split eq "yes2")
                        {
                            $delete_first = "yes";
                            if ($count_split eq '2')
                            {
                                $split = "yes4";
                            }
                        }
                        if ($split eq "yes3")
                        {
                            $delete_second = "yes";
                            if ($count_split eq '3')
                            {
                                $split = "yes4";
                            }
                        }
                        if ($split eq "yes4")
                        {
                            $delete_third = "yes";
                        }
                        goto SPLIT;               
                    }
                    if ($type eq "chloro222" && length($best_extension2) > 20)
                    {
                        my $best_extension2_reverse2 = $best_extension2;
                        $best_extension2_reverse2 =~ tr/ATCG/TAGC/;
                        my $best_extension2_reverse = reverse($best_extension2_reverse2);
                        
                        $best_extension2_reverse =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        
                        my $read_cp = $read;
                        my $read_end_rev_tmp = reverse($read_end);
                        $read_end_rev_tmp =~ tr/ATCG/TAGC/;
                        if (length($read) > $genome_range_low)
                        {
                            $read_cp = substr $read, 5000;
                            $read_cp =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        }
                        
                        my $found_seq_cp = $read_cp =~ s/$best_extension2_reverse/$best_extension2_reverse/;
                        my $found_seq_cp3 = $read_cp =~ s/$read_end_rev_tmp/$read_end_rev_tmp/;
                        
                        if ($found_seq_cp > 0 && $found_seq_cp3 > 0)
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nDELETE BEST EXTENSION 2 (CP)\n\n";
                            }
                            if ($split eq "yes2")
                            {
                                $delete_first = "yes";
                                if ($count_split eq '2')
                                {
                                    $split = "yes4";
                                }
                            }
                            if ($split eq "yes3")
                            {
                                $delete_second = "yes";
                                if ($count_split eq '3')
                                {
                                    $split = "yes4";
                                }
                            }
                            if ($split eq "yes4")
                            {
                                $delete_third = "yes";
                            }
                            goto SPLIT;                      
                        }
                    }
                    my $end_SNR = substr $read_end, -4;
                    my $GGGG = $end_SNR =~ tr/G/G/;
                    my $TTTT = $end_SNR =~ tr/T/T/;
                    my $CCCC = $end_SNR =~ tr/C/C/;
                    my $AAAA = $end_SNR =~ tr/A/A/;
                    if ($GGGG eq '7' || $TTTT eq '7' || $CCCC eq '7' || $AAAA eq '7')
                    {
                        $GGGG = $best_extension2 =~ tr/G/G/;
                        $TTTT = $best_extension2 =~ tr/T/T/;
                        $CCCC = $best_extension2 =~ tr/C/C/;
                        $AAAA = $best_extension2 =~ tr/A/A/;
                        if ($GGGG eq length($best_extension2) || $TTTT eq length($best_extension2) || $CCCC eq length($best_extension2) || $AAAA eq length($best_extension2))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nDELETE BEST EXTENSION 2/3/4 (SNR)\n\n";
                            }
                            if ($split eq "yes2")
                            {
                                $delete_first = "yes";
                                if ($count_split eq '2')
                                {
                                    $split = "yes4";
                                }
                            }
                            if ($split eq "yes3")
                            {
                                $delete_second = "yes";
                                if ($count_split eq '3')
                                {
                                    $split = "yes4";
                                }
                            }
                            if ($split eq "yes4")
                            {
                                $delete_third = "yes";
                            }
                            goto SPLIT;    
                        }             
                    }

                    my $contigs_end2 = substr $best_extension2, 0, 7;
                    $contigs_end2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    my $contigs_end0 = substr $read_end, -5;
    
                    if (exists($contigs_end{$contigs_end0.$contigs_end2}))
                    {
                        $tree{$id} = $contigs_end{$contigs_end0.$contigs_end2};
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nDELETE BEST EXTENSION 2 (CONTIG_END)\n\n";
                        }
                        $contig_end_check ="yes";
                        if ($split eq "yes2")
                        {
                            $delete_first = "yes";
                            if ($count_split eq '2')
                            {
                                $split = "yes4";
                            }
                        }
                        if ($split eq "yes3")
                        {
                            $delete_second = "yes";
                            if ($count_split eq '3')
                            {
                                $split = "yes4";
                            }
                        }
                        if ($split eq "yes4")
                        {
                            $delete_third = "yes";
                        }
                        goto SPLIT;  
                    }
                 
                    my $hasdot = $best_extension_tmp =~ tr/\./\./;
                    if (length($best_extension_tmp) > 9 && $hasdot < 2 && $hasdot eq "yesssssjeij")
                    {
                        my $end_tmp = substr $read_end, 10;
                        if (length($best_extension_tmp) < 15)
                        {
                            $end_tmp = substr $read_end, length($best_extension_tmp)-5;
                        }
                        $end_tmp = $end_tmp.$best_extension_tmp;
                        my $s = '0';
                        my $foundit = "";
                        while ($s < length($end_tmp)-$overlap)
                        {
                            my $end_tmp_d = substr $end_tmp, -($s+$overlap), $overlap;
                            
                            if ($containX_short_end2 > 0)
                            {
                                my $star = $end_tmp_d =~ tr/\*//;
                        
                                $end_tmp_d = substr $end_tmp, -($s+$overlap+($star*2)), $overlap+($star*2);
                                my $star2 = $end_tmp_d =~ tr/\*//;                                                
                                while ($star2 > $star)
                                {
                                    $end_tmp_d = substr $end_tmp, -($s+$overlap+($star*2)), $overlap+($star*2);
                                    $star = $star2;
                                    $star2 = $end_tmp_d =~ tr/\*//;
                                }
                            }    
                            my %end_tmp_d = build_partial3b $end_tmp_d;
                            foreach my $end_tmp_d (keys %end_tmp_d)
                            {
                                if (exists($hash2b{$end_tmp_d}))
                                {
                                    $foundit = "yes";
                                }
                                elsif (exists($hash2c{$end_tmp_d}))
                                {
                                    $foundit = "yes";
                                }
                            }
                            $s++;
                        }
                        if ($foundit ne "yes")
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nDELETE BEST EXTENSION 2/3/4 (no reverse match)\n\n";
                            }
                            if ($split eq "yes2")
                            {
                                $delete_first = "yes";
                                if ($count_split eq '2')
                                {
                                    $split = "yes4";
                                }
                            }
                            if ($split eq "yes3")
                            {
                                $delete_second = "yes";
                                if ($count_split eq '3')
                                {
                                    $split = "yes4";
                                }
                            }
                            if ($split eq "yes4")
                            {
                                $delete_third = "yes";
                            }
                            goto SPLIT;    
                        
                        }
                    }
                }              
                if ($split eq "yes2")
                {
                    delete $seed{$id};
                    $id      = "c".$position."_".$id;
                    $id_split1 = $id;
                }
                elsif ($split eq "yes3")
                {
                    $id      = "d".$position."_".$id;
                    $id_split2 = $id;
                }
                elsif ($split eq "yes4")
                {
                    $id      = "e".$position."_".$id;
                    $id_split3 = $id;
                }
                $position{$id} = $position;
                $position_back{$id} = $position_back;
                $contig_count{$id} = $contig_count{$id_original};
                my $count_contig_tmp = $contig_count;
                while ($count_contig_tmp > 0)
                {
                    $contig_gap_min{$id."_".$count_contig_tmp} = $contig_gap_min{$id_original."_".$count_contig_tmp};
                    $contig_gap_max{$id."_".$count_contig_tmp} = $contig_gap_max{$id_original."_".$count_contig_tmp};
                    $count_contig_tmp--;    
                } 
                if (exists($noback{$id_original}))
                {
                    $noback{$id} = $noback;
                }
                if (exists($old_id{$id_original}))
                {
                    $old_id{$id} = $old_id{$id_original};
                }
                if (exists($old_rep{$id_original}))
                {
                    $old_rep{$id} = $old_rep{$id_original};
                }
                if (exists($old_rep_old{$id_original}))
                {
                    $old_rep_old{$id} = $old_rep_old{$id_original};
                }
                if (exists($nosecond{$id_original}))
                {
                    $nosecond{$id} = undef;
                }
                if (exists($seed_split{$id_original}))
                {
                    $seed_split{$id} = undef;
                }
                if (exists($rep_return{$id_original}))
                {
                    $rep_return{$id} = $rep_return{$id_original};
                }
                if (exists($rep_return_back{$id_original}))
                {
                    $rep_return_back{$id} = $rep_return_back{$id_original};
                }
            }
            elsif ($split eq "yes5")
            {
                $best_extension1 = $best_extension;
                
                
                my $contig_id1_prev = $id;
                my $best_extension1_tmp = $best_extension1;
                $best_extension1_tmp =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;

                $contig_id1 = "1".length($read);
                $contig_id1{$contig_id1_prev} = $contig_id1;
                $contig_id1_prev = $contig_id1;
                $contig_read1 = substr $read, -$insert_size;
                $contig_read1 = $contig_read1.$best_extension1_tmp;

                $contig_id1 = $contig_id1{$id};
                
                if ($y > $startprint2)
                {
                    print OUTPUT5 "GROUP1\n";
                    foreach my $extensions_tmp (@extensions_group1)
                    {  
                        print OUTPUT5 $extensions_tmp."\n";
                    }
                    print OUTPUT5 $best_extension1." BEST_EXTENSION1\n\n";
                }
                my $contigs_end0;
                if ($before eq "yes")
                {
                    if ((length($best_extension1) < 3 || (length($best_extension1) < 6 && ($ext > 15 || $SNR_read ne ""))) && $repetitive_detect eq "" && $repetitive_detect eq "yesff")
                    {               
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nDELETE BEST EXTENSION1\n\n";
                        }
                        delete $seed{$id};
                        if ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes" && $last_chance eq "yes" && $noback eq "stop")
                        {
                            $delete_first = "yes2";
                            $delete_second = "yes2";
                            $delete_third = "yes2";
                            goto FINISH;                
                        }
                        elsif ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes")
                        {
                            $best_extension = "";
                            $delete_first = "yes2";
                            $delete_second = "yes2";
                            $delete_third = "yes2";
                            goto AFTER_EXT;
                        }
                        elsif (($delete_first eq "yes" && $delete_second eq "yes") || ($delete_first eq "yes" && $delete_third eq "yes") || ($delete_second eq "yes" && $delete_third eq "yes"))
                        {
                            goto SEED;  
                        }
                        else
                        {
                            $best_extension = "";
                        }               
                    }
                
                    if ($type eq "chloro222" && length($best_extension1) > 20 )
                    {
                        my $read_cp = $read;
                        my $read_end_rev_tmpb = $read_end.$best_extension1;
                        my $read_end_rev_tmp = reverse($read_end_rev_tmpb);
                        $read_end_rev_tmp =~ tr/ATCG/TAGC/;
                        $read_end_rev_tmp =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        if (length($read) > 80000)
                        {
                            $read_cp = substr $read, -80000;
                        }
                        
                        my $found_seq_cp3 = $read_cp =~ s/$read_end_rev_tmp/$read_end_rev_tmp/;
                        
                        if (($found_seq_cp3 > 0) && $before eq "yes")
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nDELETE BEST EXTENSION1 (CP)\n\n";
                            }
                            delete $seed{$id};
                            if ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes" && $last_chance eq "yes" && $noback eq "stop")
                            {
                                $delete_first = "yes2";
                                $delete_second = "yes2";
                                $delete_third = "yes2";
                                goto FINISH;                
                            }
                            elsif ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes")
                            {
                                $best_extension = "";
                                $delete_first = "yes2";
                                $delete_second = "yes2";
                                $delete_third = "yes2";
                                goto AFTER_EXT;
                            }
                            elsif (($delete_first eq "yes" && $delete_second eq "yes") || ($delete_first eq "yes" && $delete_third eq "yes") || ($delete_second eq "yes" && $delete_third eq "yes"))
                            {
                                goto SEED;  
                            }
                            else
                            {
                                $best_extension = "";
                            }                        
                        }                    
                    }
                    my $end_SNR = substr $read_end, -4;
                    my $GGGG = $end_SNR =~ tr/G/G/;
                    my $TTTT = $end_SNR =~ tr/T/T/;
                    my $CCCC = $end_SNR =~ tr/C/C/;
                    my $AAAA = $end_SNR =~ tr/A/A/;
                    if ($GGGG eq '7' || $TTTT eq '7' || $CCCC eq '7' || $AAAA eq '7' && $before eq "yessss")
                    {
                        $GGGG = $best_extension1 =~ tr/G/G/;
                        $TTTT = $best_extension1 =~ tr/T/T/;
                        $CCCC = $best_extension1 =~ tr/C/C/;
                        $AAAA = $best_extension1 =~ tr/A/A/;
                        if ($GGGG eq length($best_extension1) || $TTTT eq length($best_extension1) || $CCCC eq length($best_extension1) || $AAAA eq length($best_extension1))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nDELETE BEST EXTENSION 1 (SNR)\n\n";
                            }
                            delete $seed{$id};
                            if ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes" && $last_chance eq "yes" && $noback eq "stop")
                            {
                                $delete_first = "yes2";
                                $delete_second = "yes2";
                                $delete_third = "yes2";
                                goto FINISH;                
                            }
                            elsif ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes")
                            {
                                $best_extension = "";
                                $delete_first = "yes2";
                                $delete_second = "yes2";
                                $delete_third = "yes2";
                                goto AFTER_EXT;
                            }
                            elsif (($delete_first eq "yes" && $delete_second eq "yes") || ($delete_first eq "yes" && $delete_third eq "yes") || ($delete_second eq "yes" && $delete_third eq "yes"))
                            {
                                goto SEED;  
                            }
                            else
                            {
                                $best_extension = "";
                            }
                        }             
                    }

                    my $contigs_end1 = substr $best_extension1, 0, 7;
                    $contigs_end1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    my $contigs_end0 = substr $read_end, -5;
    
                    if (exists($contigs_end{$contigs_end0.$contigs_end1}))
                    {
                        $tree{$id} = $contigs_end{$contigs_end0.$contigs_end1};
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nDELETE BEST EXTENSION 1 (CONTIG_END)\n\n";
                        }
                        my $contigs_end1 = substr $best_extension1, 0, 10;
                        $contigs_end1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        my $contigs_end0 = substr $read_end, -15;
                        my $repetitive_test = $contigs_end0.$contigs_end1;
                        my $end_repetitive = $read_short_end2;
                        my $check_repetitive = $end_repetitive =~ s/$repetitive_test/$repetitive_test/g;
                        if ($check_repetitive > 1)
                        {
                            $delete_first = "yes";
                            $delete_second = "yes";
                            $delete_third = "yes";
                        }
                        if ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes")
                        {
                            $contig_end = "yes";
                            $delete_first = "yes2";
                            $delete_second = "yes2";
                            $delete_third = "yes2";
                            goto INDEL;
                        }
                        elsif (($delete_first eq "yes" && $delete_second eq "yes") || ($delete_first eq "yes" && $delete_third eq "yes") || ($delete_second eq "yes" && $delete_third eq "yes"))
                        {
                            delete $seed{$id};
                            goto SEED;  
                        }
                        else
                        {
                            delete $seed{$id};
                            $best_extension = "";
                        } 
                    }
                
                    my $hasdot = $best_extension1 =~ tr/\./\./;
                    if (length($best_extension1) > 9 && $hasdot < 2 && $before eq "yesss")
                    {
                        my $end_tmp = substr $read_end, 10;
                        if (length($best_extension1) < 15)
                        {
                            $end_tmp = substr $read_end, length($best_extension1)-5;
                        }
                        $end_tmp = $end_tmp.$best_extension1;
                        my $s = '0';
                        my $foundit = "";
                        while ($s < length($end_tmp)-$overlap)
                        {
                            my $end_tmp_d = substr $end_tmp, -($s+$overlap), $overlap;
                            
                            if ($containX_short_end2 > 0)
                            {
                                my $star = $end_tmp_d =~ tr/\*//;
                        
                                $end_tmp_d = substr $end_tmp, -($s+$overlap+($star*2)), $overlap+($star*2);
                                my $star2 = $end_tmp_d =~ tr/\*//;                                                
                                while ($star2 > $star)
                                {
                                    $end_tmp_d = substr $end_tmp, -($s+$overlap+($star*2)), $overlap+($star*2);
                                    $star = $star2;
                                    $star2 = $end_tmp_d =~ tr/\*//;
                                }
                            }    
                            my %end_tmp_d = build_partial3b $end_tmp_d;
                            foreach my $end_tmp_d (keys %end_tmp_d)
                            {
                                if (exists($hash2b{$end_tmp_d}))
                                {
                                    $foundit = "yes";
                                }
                                elsif (exists($hash2c{$end_tmp_d}))
                                {
                                    $foundit = "yes";
                                }
                            }
                            $s++;
                        }
                        if ($foundit ne "yes")
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nDELETE BEST EXTENSION 1 (no reverse match)\n\n";
                            }
                            delete $seed{$id};
                            if ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes" && $last_chance eq "yes" && $noback eq "stop")
                            {
                                $delete_first = "yes2";
                                $delete_second = "yes2";
                                $delete_third = "yes2";
                                goto FINISH;                
                            }
                            elsif ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes")
                            {
                                $best_extension = "";
                                $delete_first = "yes2";
                                $delete_second = "yes2";
                                $delete_third = "yes2";
                                goto AFTER_EXT;
                            }
                            elsif (($delete_first eq "yes" && $delete_second eq "yes") || ($delete_first eq "yes" && $delete_third eq "yes") || ($delete_second eq "yes" && $delete_third eq "yes"))
                            {
                                goto SEED;  
                            }
                            else
                            {
                                $best_extension = "";
                            }   
                        }  
                    }
                    if ($repetitive_detect ne "" && $contig_end_check eq "yes")
                    {
                        my $contigs_end1 = substr $best_extension1, 0, 10;
                         $contigs_end1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        my $contigs_end0 = substr $read_end, -15;
                        my $repetitive_test = $contigs_end0.$contigs_end1;
                        my $end_repetitive = $read_short_end2;
                        my $check_repetitive = $end_repetitive =~ s/$repetitive_test/$repetitive_test/g;
                        if ($check_repetitive > 1)
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nDELETE BEST EXTENSION 1 (CONTIG_END+REPETITIVE)\n\n";
                            }
                            if ($delete_first eq "yes" && $delete_second eq "yes" && $delete_third eq "yes")
                            {
                                $delete_first = "yes2";
                                $delete_second = "yes2";
                                $delete_third = "yes2";
                                goto INDEL;
                            }
                            elsif (($delete_first eq "yes" && $delete_second eq "yes") || ($delete_first eq "yes" && $delete_third eq "yes") || ($delete_second eq "yes" && $delete_third eq "yes"))
                            {
                                delete $seed{$id};
                                goto SEED;  
                            }
                            else
                            {
                                delete $seed{$id};
                                $best_extension = "";
                            }
                        }
                    }
                }
                
                my $read_part_tmp = substr $read_short_end2, -$read_length;
                my $star_check = $read_part_tmp =~ tr/\*/\*/;
                if ($SNP_active eq "yes" && ($before eq "yes" || $extensions_before eq "yes") && $count_split eq '2' && $star_check eq 0)
                {
                    my @chars;
                    my @chars2;
                    undef @chars;
                    undef @chars2;
                    my $p = '0';
                    my $amatch = '0';
                    my $nomatch = '0';
                    my $best_extension_short = "";
                    my $best_extension_long = "";
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "CHECK_DELETION\n";
                    }
                    if (length($best_extension1) >= length($best_extension2))
                    {
                        $best_extension_short = $best_extension2;
                        $best_extension_long = $best_extension1;
                        $best_extension_short =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        $best_extension_long =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        @chars = split //, $best_extension_short;
                        @chars2 = split //, $best_extension_long;
                    }
                    elsif (length($best_extension2) > length($best_extension1))
                    { 
                        $best_extension_short = $best_extension1;
                        $best_extension_long = $best_extension2;
                        $best_extension_short =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        $best_extension_long =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        @chars = split //, $best_extension_short;
                        @chars2 = split //, $best_extension_long;
                    }
                    while ($p < length($best_extension_short) && $p < 25)
                    {
                        my $i = '1';
                        $amatch = '0';
                        $nomatch = '0';
                        while ($i < @chars-$p)
                        {
                            if ($chars[$i+$p] eq $chars2[$i-1])
                            {
                                $amatch++;   
                            }
                            elsif ($chars[$i+$p] eq ".")
                            {                            
                            }
                            else
                            {                            
                                $nomatch++;
                            }
                            $i++;
                        }
                        $p++;
                        if ($amatch > ($nomatch*8) && (($amatch > 2 && $p < 3) || $amatch > 7))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "DELETION\n";
                            }
                            $deletion = "yes";
                                                  
                            my $indel = substr $best_extension_short,0, $p;
                            $indel =~ s/A/A*/g;
                            $indel =~ s/T/T*/g;
                            $indel =~ s/G/G*/g;
                            $indel =~ s/C/C*/g;
                            
                            my $after_indel;
                            if ($nomatch > 0)
                            {
                                my $after_indel1 = substr $best_extension_short, $p;
                                my @chars3 = split //, $after_indel1;
                                my $f = '0';
                                
INDEL0:                         while ($f < @chars3)
                                {
                                    if ($chars3[$f] eq $chars2[$f])
                                    {
                                        $after_indel .= $chars3[$f];
                                    }
                                    elsif ($f < 20)
                                    {
                                        $after_indel .= ".";
                                    }
                                    else
                                    {
                                        last INDEL0;
                                    }
                                    $f++;
                                }
                            }
                            else
                            {
                                my $after_indel1 = substr $best_extension_short, $p;
                                my @chars3 = split //, $after_indel1;
                                my $f = '0';
                                
INDEL1:                         while ($f < @chars3)
                                {
                                    if ($chars3[$f] eq $chars2[$f])
                                    {
                                        $after_indel .= $chars3[$f];
                                    }
                                    else
                                    {
                                        last INDEL1;
                                    }
                                    $f++;
                                }
                            }
                            $best_extension = $indel.$after_indel;
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 $best_extension." BEST_EXTENSION\n";
                            }
                            if ($after_indel eq "")
                            {
                                $deletion = "";
                                $best_extension = "";
                            }
                            goto INDELa;
                        }  
                    }
                    $p = '0';
                    while ($p < length($best_extension_short) && $p < 25)
                    {
                        my $i = '1';
                        $amatch = '0';
                        $nomatch = '0';
                        while ($i <= @chars && $i <= @chars2-$p)
                        {
                            if ($chars[$i-1] eq $chars2[$i+$p])
                            {
                                $amatch++;
                            }
                            elsif ($chars[$i-1] eq ".")
                            {                            
                            }
                            else
                            {
                                $nomatch++;
                            }
                            $i++;
                        }
                        $p++;
                        if ($amatch > ($nomatch*8) && (($amatch > 2 && $p < 3) || $amatch > 7))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "DELETION\n";
                            }
                            $deletion = "yes";
                           
                            my $indel = substr $best_extension_long,0, $p;
                            my $indel2 = $indel;
                            $indel =~ s/A/A*/g;
                            $indel =~ s/T/T*/g;
                            $indel =~ s/G/G*/g;
                            $indel =~ s/C/C*/g;
                            
                            my $after_indel;
                            if ($nomatch > 0)
                            {
                                my $after_indel1 = substr $best_extension_long, $p;
                                my @chars3 = split //, $after_indel1;
                                my $f = '0';
                                
INDEL2:                         while ($f < @chars3)
                                {
                                    if ($chars3[$f] eq $chars[$f])
                                    {
                                        $after_indel .= $chars3[$f];
                                    }
                                    elsif ($f < 20)
                                    {
                                        $after_indel .= ".";
                                    }
                                    else
                                    {
                                        last INDEL2;
                                    }
                                    $f++;
                                }
                            }
                            else
                            {
                                my $after_indel1 = substr $best_extension_long, $p;
                                my @chars3 = split //, $after_indel1;
                                my $f = '0';
                                
INDEL3:                         while ($f < @chars3)
                                {
                                    if ($chars3[$f] eq $chars[$f])
                                    {
                                        $after_indel .= $chars3[$f];
                                    }
                                    else
                                    {
                                        last INDEL3;
                                    }
                                    $f++;
                                }
                            }
                            $best_extension = $indel.$after_indel;
                            
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 $best_extension." BEST_EXTENSION\n";
                            }
                            if ($after_indel eq "")
                            {
                                $deletion = "";
                                $best_extension = "";
                            }
                            goto INDELa;
                        }  
                    }
                }
                if (length($best_extension1) > 4 && length($best_extension2) > 4 && $cp_input ne "" && $type eq "mito_plant")
                {
                    my $p = -$overlap;
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "CHECK_REFERENCE\n\n";
                    }                   
                    while ($p > (-$overlap*2))
                    {
                        my $ref_part2 = substr $read_short_end2, $p, $overlap;
                        my %ref_part = build_partial3b $ref_part2, "";
                        
                        foreach my $ref_part (keys %ref_part)
                        {
                            if (exists($hashref{$ref_part}))
                            {                            
                                my $ref_loc = $p + $overlap;
                                
                                my $ref_id3 = $hashref{$ref_part};      
                                my $ref_id2 = substr $ref_id3, 1;
                                my @ref_id3 = split /,/,$ref_id2;
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 $ref_part." EXISTS ".$ref_loc." LOC ".$ref_id3." LOC_REF\n";
                                }
                           
                                foreach my $ref_id (@ref_id3)
                                {
                                    my $prev_loc1 = $ref_id + $overlap - $ref_loc;
                                    
                                    if (exists($hashref2{$prev_loc1}))
                                    {
                                        my $ref_check = $hashref2{$prev_loc1};

                                        print OUTPUT5 $ref_check." EXISTSREF1\n";
   
                                        my $best_extension1_part = substr $best_extension1, 0, 30;
                                        my $best_extension2_part = substr $best_extension2, 0, 30;
                                        my @ref1 = build_partial $best_extension1_part;
                                        
                                        foreach my $best_extension1_part2 (@ref1)
                                        {                                    
                                            my $ref_check_tmp = $ref_check;
                                            if ($ref_check_tmp =~ s/$best_extension1_part2/$best_extension1_part2/)
                                            {    
                                                print OUTPUT5 "REFERNCE_GUIDED\n";
                                                print OUTPUT5 $best_extension2." BEST_EXTENSION2\n\n";
                                                $best_extension = $best_extension2;
                                                $reference_guided = "yes";
                                            }
                                        }
                                        my @ref2 = build_partial $best_extension2_part;
                                        foreach my $best_extension2_part2 (@ref2)
                                        {
                                            my $ref_check_tmp = $ref_check;
                                            if ($ref_check_tmp =~ s/$best_extension2_part2/$best_extension2_part2/)
                                            {                                           
                                                $best_extension = $best_extension1;
                                                if ($reference_guided eq "")
                                                {
                                                    $reference_guided = "yes";
                                                }
                                                else
                                                {
                                                    $reference_guided = "";
                                                }
                                                print OUTPUT5 $reference_guided." REFERNCE_GUIDED\n";
                                                print OUTPUT5 $best_extension1." BEST_EXTENSION1\n\n";
                                            }
                                        }
                                        if ($reference_guided eq "yes")
                                        {
                                            goto INDELa;
                                        }
                                        else
                                        {
                                            $best_extension = "";
                                        }
                                    }
                                }                        
                            }
                        }
                        $p--;
                    }
                }

INDELa:                
                if (($last_chance ne "yes" || length($read) <= $insert_size) && $before eq "yes" && $indel_split_skip ne "yes" && $ext_total > 15 && ($delete_first ne "yes" || (length($best_extension1) < 5 && length($best_extension2) < 5)) && $indel_split < ($read_length-25 -$overlap))
                {
                    delete $seed{$id};
                    delete $seed{$id_split1};
                    delete $seed{$id_split2};
                    delete $seed{$id_split3};
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "\nDELETE SPLIT\n\n";
                    }
                }
                elsif ((($before eq "yes" && $SNP_active eq "") || $extensions_before eq "yes") && ($delete_first ne "yes" ||  (length($best_extension1) < 5 && length($best_extension2) < 5 && length($best_extension3) < 5 && length($best_extension4) < 5)))
                {
                    delete $seed{$id};
                    delete $seed{$id_split1};
                    delete $seed{$id_split2};
                    delete $seed{$id_split3};
                    $SNP_active{$id_original} = "yes";
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "\nDELETE SPLIT2\n\n";
                    }
                }
INDEL:
                $id             = $id_original;
                $split          = "";
                $split_forward  = "yes";
                

                if ($jump_rep eq "yes" || ($extensions_before ne "yes" && $contig_end ne "yes" && $indel_split_skip ne "yes" && ($delete_first eq "" || $delete_second eq "" || $delete_third eq "") && $before ne "yes" && $reference_guided ne "yes"))
                {                   
                    my $before_shorter = "";
                    my $overhangb = 1+($read_length/40);
                    my $count_one = '0';
                    my $count_two = '0';
                    while (($count_one < 2 || $count_two < 2) && $overhangb < 15)
                    {
                        $overhangb++;
                        foreach my $one (@extensions_group1)
                        {
                            if (length($one) < $overhangb && length($one) > 0)
                            {
                                $count_one++;
                            }
                        }
                        foreach my $two (@extensions_group2)
                        {
                            if (length($two) < $overhangb && length($two) > 0)
                            {
                                $count_two++;
                            }
                        }
                    }
                    my $overhang = sprintf("%.0f", $overhangb);
                    delete $jump_rep{$id};
BEFORE:
                    my $s = '0';
                    
                    if ($before_shorter eq "yes")
                    {
                        $overhang += 5;
                    }
                    $before_shorter = "";

                    print OUTPUT5 $overhang." OVERHANG\n";
                    my @extensions_yuyu;
                    undef @extensions_yuyu;
                    my %extensions_yuyu;
                    undef %extensions_yuyu;
                    my %before1F;
                    my %before2F;
                    my %before3F;
                    my %before4F;
                    my %before1B;
                    my %before2B;
                    my %before3B;
                    my %before4B;
                    undef %repetitive_pair;
                    my $skip_overhang = "";
                    if ($AT_rich_before eq "yes")
                    {
                        my $before_split_tmp = substr $read_short_end2, -($read_length-$overhang);
                        my $A_rich_test = $before_split_tmp =~ tr/AX\.//;
                        my $T_rich_test = $before_split_tmp =~ tr/TX\.//;
                        my $G_rich_test = $before_split_tmp =~ tr/GX\.//;
                        my $C_rich_test = $before_split_tmp =~ tr/CX\.//;
                        if ($A_rich_test > ($read_length-$overhang-3) || $T_rich_test > ($read_length-$overhang-3) || $G_rich_test > ($read_length-$overhang-3) || $C_rich_test > ($read_length-$overhang-3))
                        {
                            $skip_overhang = "yes";
                            $overhang -= 5;
                        }
                    }

                    my $before_split = substr $read_short_end2, -($read_length-$right-1), $overhang+$overlap;
                    my $star0 = '0';
                    if ($containX_short_end2 > 0)
                    {
                        my $before_split2 = substr $read_short_end2, -($read_length-$right-1);
                        $star0 = $before_split2 =~ tr/\*//;
                        if ($star0 > 0)
                        {
                            my $star_tmp = $before_split2 =~ tr/\*//;
                            $before_split = substr $read_short_end2, -($read_length-$right-1+($star_tmp*2)), $overhang+$overlap+($star_tmp*2);
                            $before_split2 = substr $read_short_end2, -($read_length-$right-1+($star_tmp*2));
                            my $star2_tmp = $before_split2 =~ tr/\*//;                                                
                            while ($star2_tmp > $star_tmp)
                            {
                                $before_split = substr $read_short_end2, -($read_length-$right-1+($star_tmp*2)+(($star2_tmp-$star_tmp)*2)), $overhang+$overlap+($star_tmp*2)+(($star2_tmp-$star_tmp)*2);
                                $before_split2 = substr $read_short_end2, -($read_length-$right-1+($star_tmp*2)+(($star2_tmp-$star_tmp)*2));
                                $star_tmp = $star2_tmp;
                                $star2_tmp = $before_split2 =~ tr/\*//;
                                $star0 = $star2_tmp;
                            }   
                        }
                    }
                    $before_split =~ tr/ACTG/TGAC/;
                    my $before_split_reverse = reverse($before_split);
                    while ($s <= length($before_split_reverse)-$overlap)
                    {
                        my $line_tmp = substr $before_split_reverse, $s, $overlap;
                        if ($star0 > 0)
                        {
                            my $star = $line_tmp =~ tr/\*//;
                            $line_tmp = substr $before_split_reverse, $s, $overlap+($star*2);
                            my $star2 = $line_tmp =~ tr/\*//;                                                
                            while ($star2 > $star)
                            {
                                $line_tmp = substr $before_split, $s, $overlap+($star*2)+(($star2-$star)*2);
                                $star = $star2;
                                $star2 = $line_tmp =~ tr/\*//;
                            }   
                        }
                        if ($containX_short_end2 > 0)
                        {
                            my $after_before_split_tmp = substr $read_short_end2, -($read_length-$right-1+($star0*2));
                        }
                        my %line_tmp = build_partial3b $line_tmp, "reverse";
                        foreach my $line (keys %line_tmp)
                        {
                            if (exists($hash2c{$line}))
                            {                       
                                my $search = $hash2c{$line};
                                $search = substr $search, 1;
                                my @search = split /,/,$search;
                                                                            
                                foreach my $search (@search)
                                {                             
                                    my $search_tmp = substr $search, 0, -1;
                                    my $search_end = substr $search, -1;
                                    if (exists($hash{$search_tmp}))
                                    {
                                        my @search_tmp = split /,/,$hash{$search_tmp};
                                        my $found;
                                        if ($search_end eq "1")
                                        {
                                            $found = $search_tmp[0];                                            
                                        }
                                        elsif ($search_end eq "2")
                                        {
                                            $found = $search_tmp[1];                                           
                                        }
                                        if ($encrypt eq "yes")
                                        {
                                            $found = decrypt $found;                                          
                                        }
                                        my $found_new;
                                        my $h;
                                        if (($overhang-$s-($star0*2)-5)<0)
                                        {
                                            $h = '0';
                                        }
                                        else
                                        {
                                            $h = ($overhang-$s-($star0*2)-5);
                                        }
                                        my $last_10 = substr $found, $h, 30+($star0*2)+5;
                                        my $first_nuc;
                                        my $found_reverse;
                                        $last_10 =~ tr/ACTG/TGAC/;
                                        my $last_10r = reverse($last_10);
                                        my $last_20_read_end = substr $read_end, -20;
                                        my $check_last10 = $last_10r =~ s/(.)$last_20_read_end/$last_20_read_end/;
                                        if ($check_last10 > 0)
                                        {
                                            $found =~ tr/ACTG/TGAC/;
                                            my $found_reverse2 = reverse($found);
                                            if ($found_reverse2 =~ m/(.*.$last_20_read_end)((.).*)$/)
                                            {
                                                $first_nuc = $3;
                                                $found = $2;
                                                $found_reverse = $1;
                                            }       
                                            my $first_nuc1 = substr $best_extension1, 0, 1;
                                            my $first_nuc2 = substr $best_extension2, 0, 1;
                                            my $first_nuc3 = substr $best_extension3, 0, 1;
                                            my $first_nuc4 = substr $best_extension4, 0, 1;
                                            
                                            my $extension_yuyu = $found;
                                            $extensions_yuyu{$search} = $extension_yuyu;
    
                                            if ($first_nuc eq $first_nuc1 && $check_last10 > 0)
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                my $found_tmp = reverse($found_reverse);
                                                $filter_before1_pair{$search} = $found_tmp;
                                                $before1B{$search} = $found_tmp;                                 
                                            }
                                            elsif ($first_nuc eq $first_nuc2 && $check_last10 > 0)
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                my $found_tmp = reverse($found_reverse);
                                                $filter_before2_pair{$search} = $found_tmp;
                                                $before2B{$search} = $found_tmp;
                                            }
                                            elsif ($first_nuc eq $first_nuc3 && $check_last10 > 0)
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                my $found_tmp = reverse($found_reverse);
                                                $filter_before3_pair{$search} = $found_tmp;
                                                $before3B{$search} = $found_tmp;
                                            }
                                            elsif ($first_nuc eq $first_nuc3 && $check_last10 > 0)
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                my $found_tmp = reverse($found_reverse);
                                                $filter_before4_pair{$search} = $found_tmp;
                                                $before4B{$search} = $found_tmp;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        $s++;
                    }
                    $s = '0';
                    $before_split = substr $read_short_end2, -($read_length-$left-1), $overhang+$overlap;
                    my $star3 = '0';
                    if ($containX_short_end2 > 0)
                    {
                        my $before_split2 = substr $read_short_end2, -($read_length-$left-1);
                        $star3 = $before_split2 =~ tr/\*//;
                        if ($star3 > 0)
                        {
                            $before_split = substr $read_short_end2, -($read_length-$left-1+($star3*2)), $right+$overlap+($star3*2);
                        }
                    }
                    while ($s <= length($before_split)-$overlap)
                    {
                        my $line_tmp = substr $before_split, $s, $overlap;
                        if ($star3 > 0)
                        {
                            my $star = $line_tmp =~ tr/\*//;
                            $line_tmp = substr $before_split, $s, $overlap+($star*2);
                            my $star2 = $line_tmp =~ tr/\*//;                                                
                            while ($star2 > $star)
                            {
                                $line_tmp = substr $before_split, $s, $overlap+($star*2)+(($star2-$star)*2);
                                $star = $star2;
                                $star2 = $line_tmp =~ tr/\*//;
                            }   
                        }
                        my %line_tmp = build_partial3b $line_tmp, "";
                        foreach my $line (keys %line_tmp)
                        {
                            if (exists($hash2b{$line}))
                            {                       
                                my $search = $hash2b{$line};
                                my $search_rev;
                                $search = substr $search, 1;
                                my @search = split /,/,$search;
                                                                            
                                foreach my $search (@search)
                                {                             
                                    my $search_tmp = substr $search, 0, -1;
                                    my $search_end = substr $search, -1;
                                    if (exists($hash{$search_tmp}))
                                    {
                                        my @search_tmp = split /,/,$hash{$search_tmp};
                                        my $found;
                                        my $found_rev;
                                        if ($search_end eq "1")
                                        {
                                            $found = $search_tmp[0];
                                            $found_rev = $search_tmp[1];
                                            $search_rev = $search_tmp."2";
                                        }
                                        elsif ($search_end eq "2")
                                        {
                                            $found = $search_tmp[1];
                                            $found_rev = $search_tmp[0];
                                            $search_rev = $search_tmp."1";
                                        }
                                        if ($encrypt eq "yes")
                                        {
                                            $found = decrypt $found;
                                            $found_rev = decrypt $found_rev;
                                        }
                                        my $found_new;
                                        my $first_nuc;
                                        my $last_10 = substr $found, -(11+$s), 10;
                                        my $last_10b = substr $found,  -(11+$s-($star3)), 10;
                                        
                                        my $last_10_read_end = substr $read_end, -9;
                                        my $check_last10 = $last_10 =~ s/(.)$last_10_read_end/$1$last_10_read_end/;
                                        my $check_last10b = $last_10b =~ s/(.)$last_10_read_end/$1$last_10_read_end/;
                                        if ($check_last10 > 0)
                                        {
                                            $found_new = substr $found, 0, -(1+$s), "";
                                            $first_nuc = substr $found, 0, 1;
                                        }
                                        elsif ($check_last10b > 0)
                                        {
                                            $found_new = substr $found, 0, -(1+$s-($star3)),"";
                                            $first_nuc = substr $found,  0, 1;
                                        }
                                        if ($check_last10 > 0 || $check_last10b > 0)
                                        {
                                            my $first_nuc1 = substr $best_extension1, 0, 1;
                                            my $first_nuc2 = substr $best_extension2, 0, 1;
                                            my $first_nuc3 = substr $best_extension3, 0, 1;
                                            my $first_nuc4 = substr $best_extension4, 0, 1;
                                            my $extension_yuyu = $found;
                                            $extensions_yuyu{$search} = $extension_yuyu;
                                            if ($first_nuc eq $first_nuc1 && ($check_last10 > 0 || $check_last10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                my $found_tmp = reverse($found_new);
                                                $before1F{$search} = $found_tmp;
                                                $found_rev  =~ tr/ACTG/TGAC/;
                                                my $found_rev2 = reverse($found_rev);
                                                $repetitive_pair{$found_rev2} = $search_rev;
                                            }
                                            elsif ($first_nuc eq $first_nuc2 && ($check_last10 > 0 || $check_last10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                my $found_tmp = reverse($found_new);
                                                $before2F{$search} = $found_tmp;
                                                $found_rev  =~ tr/ACTG/TGAC/;
                                                my $found_rev2 = reverse($found_rev);
                                                $repetitive_pair{$found_rev2} = $search_rev;
                                            }
                                            elsif ($first_nuc eq $first_nuc3 && ($check_last10 > 0 || $check_last10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                my $found_tmp = reverse($found_new);
                                                $before3F{$search} = $found_tmp;
                                                $found_rev  =~ tr/ACTG/TGAC/;
                                                my $found_rev2 = reverse($found_rev);
                                                $repetitive_pair{$found_rev2} = $search_rev;
                                            }
                                            elsif ($first_nuc eq $first_nuc4 && ($check_last10 > 0 || $check_last10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                my $found_tmp = reverse($found_new);
                                                $before4F{$search} = $found_tmp;
                                                $found_rev  =~ tr/ACTG/TGAC/;
                                                my $found_rev2 = reverse($found_rev);
                                                $repetitive_pair{$found_rev2} = $search_rev;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        $s++;
                    }
                    my %before_all1 = (%before1B, %before1F);
                    my %before_all2 = (%before2B, %before2F);
                    my %before_all3 = (%before3B, %before3F);
                    my %before_all4 = (%before4B, %before4F);


                    my $end_short_tmp = substr $read_short_end2, -($read_length+20);

BEFORE_EXTRA:       
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $end_short_tmp." END_SHORT_TMP\n\n";
                        print OUTPUT5 reverse($end_short_tmp)." REVERSE_END_SHORT_TMP\n\n";
                    }
                    my $first_yuyu = "";
                    my $second_yuyu = "";
                    my $third_yuyu = "";
                    my $fourth_yuyu = "";
                    my $first_yuyu2 = '0';
                    my $second_yuyu2 = '0';
                    my $third_yuyu2 = '0';
                    my $fourth_yuyu2 = '0';
                    my @filter_dot_before1;
                    my @filter_dot_before2;
                    my @filter_dot_before3;
                    my @filter_dot_before4;
                    undef @filter_dot_before1;
                    undef @filter_dot_before2;
                    undef @filter_dot_before3;
                    undef @filter_dot_before4;
                    my @extensions_before;
                    my @extensions_before1;
                    my @extensions_before2;
                    my @extensions_before3;
                    my @extensions_before4;
                    undef @extensions_before;
                    undef @extensions_before1;
                    undef @extensions_before2;
                    undef @extensions_before3;
                    undef @extensions_before4;
                    
                    my $end_short_tmp_part = substr $end_short_tmp, -($read_length-$overhang-5-1);
                    my $star_first = substr $end_short_tmp_part, 0, 1;
                    if ($star_first eq "*")
                    {
                        $end_short_tmp_part = substr $end_short_tmp, -($read_length-$overhang-5-1+1);
                    }  
                    
                    foreach my $search (keys %before_all2)
                    {
                        my $overhang_tmp = $overhang;
                        my $yuyu0 = $before_all2{$search};
                        my $yuyu2 = reverse($yuyu0);
                        if (length($yuyu2) >= $read_length-$overhang_tmp-5-1)
                        {
                            $yuyu2 =~ tr/N/\./;
                            my $yuyu2_tmp = $yuyu2;
                            my $check_yuyuy2 = $yuyu2_tmp =~ s/$end_short_tmp_part/$end_short_tmp_part/;

                            if ($check_yuyuy2 > 0)
                            {
                                $second_yuyu = "yes";
                                $second_yuyu2++;
                                my $end_short_tmp_part = substr $end_short_tmp, -length($yuyu2)+1;
                                my $end_short_tmp_part2;
                                my $star = $end_short_tmp_part =~ tr/\*//;
                                my $check_yuyuy2b2 = '0';
                                if ($star > 0)
                                {
                                    $end_short_tmp_part = substr $end_short_tmp, -length($yuyu2)+1-($star*2);
                                    $end_short_tmp_part2 = substr $end_short_tmp, -length($yuyu2)+1-($star);
                                    my $yuyu2_tmp2 = $yuyu2;
                                    $check_yuyuy2b2 = $yuyu2_tmp2 =~ s/.$end_short_tmp_part2/$end_short_tmp_part2/;
                                }      
                                my $yuyu2_tmp3 = $yuyu2;
                                my $check_yuyuy2b = $yuyu2_tmp3 =~ s/.$end_short_tmp_part/$end_short_tmp_part/;
                                if ($check_yuyuy2b > 0 || $check_yuyuy2b2 > 0)
                                {
                                    push @filter_dot_before2, $yuyu0;
                                    if (exists($extensions_yuyu{$search}))
                                    {     
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $yuyu0." FOUND2   ";
                                            print OUTPUT5 $extensions_yuyu{$search}." EXTe2\n";
                                        }
                                        push @extensions_before, $extensions_yuyu{$search};
                                        push @extensions_before2, $extensions_yuyu{$search};
                                    }
                                    $filter_before2{$search} = undef;
                                }
                            }
                        }
                    }
                    foreach my $search (keys %before_all1)
                    {
                        my $overhang_tmp = $overhang;
                        my $yuyu0 = $before_all1{$search};
                        my $yuyu1 = reverse($yuyu0);
                        if (length($yuyu1) >= $read_length-$overhang_tmp-5-1)
                        {
                            $yuyu1 =~ tr/N/\./;
                            my $yuyu1_tmp = $yuyu1;
                            my $check_yuyuy = $yuyu1_tmp =~ s/$end_short_tmp_part/$end_short_tmp_part/;
                            
                            if ($check_yuyuy > 0)
                            {
                                $first_yuyu2++;
                                $first_yuyu = "yes";
                                my $end_short_tmp_part = substr $end_short_tmp, -length($yuyu1)+1;
                                my $star = $end_short_tmp_part =~ tr/\*//;
                                my $check_yuyuy1b2 = '0';
                                my $end_short_tmp_part2;
                                if ($star > 0)
                                {
                                    $end_short_tmp_part = substr $end_short_tmp, -length($yuyu1)+1-($star*2);
                                    $end_short_tmp_part2 = substr $end_short_tmp, -length($yuyu1)+1-($star);
                                    my $yuyu1_tmp2 = $yuyu1;
                                    $check_yuyuy1b2 = $yuyu1_tmp2 =~ s/.$end_short_tmp_part2/$end_short_tmp_part2/;
                                }
                                my $yuyu1_tmp3 = $yuyu1;
                                my $check_yuyuy1b = $yuyu1_tmp3 =~ s/.$end_short_tmp_part/$end_short_tmp_part/;
                                if ($check_yuyuy1b > 0 || $check_yuyuy1b2 > 0)
                                {
                                    push @filter_dot_before1, $yuyu0;
                                    if (exists($extensions_yuyu{$search}))
                                    {
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $yuyu0." FOUND1   ";
                                            print OUTPUT5 $extensions_yuyu{$search}." EXTe1\n";
                                        }
                                        push @extensions_before, $extensions_yuyu{$search};
                                        push @extensions_before1, $extensions_yuyu{$search};
                                    }
                                    $filter_before1{$search} = undef;
                                }
                            }
                        }
                    }
                    foreach my $search (keys %before_all3)
                    {
                        my $overhang_tmp = $overhang;
                        my $yuyu0 = $before_all3{$search};
                        my $yuyu3 = reverse($yuyu0);
                        if (length($yuyu3) >= $read_length-$overhang_tmp-5-1)
                        {
                            $yuyu3 =~ tr/N/\./;
                            my $yuyu3_tmp = $yuyu3;
                            my $check_yuyuy = $yuyu3_tmp =~ s/$end_short_tmp_part/$end_short_tmp_part/;
                            
                            if ($check_yuyuy > 0)
                            {
                                $third_yuyu2++;
                                $third_yuyu = "yes";
                                my $end_short_tmp_part = substr $end_short_tmp, -length($yuyu3)+1;
                                my $star = $end_short_tmp_part =~ tr/\*//;
                                my $check_yuyuy3b2 = '0';
                                my $end_short_tmp_part2;
                                if ($star > 0)
                                {
                                    $end_short_tmp_part = substr $end_short_tmp, -length($yuyu3)+1-($star*2);
                                    $end_short_tmp_part2 = substr $end_short_tmp, -length($yuyu3)+1-($star);
                                    my $yuyu3_tmp2 = $yuyu3;
                                    $check_yuyuy3b2 = $yuyu3_tmp2 =~ s/.$end_short_tmp_part2/$end_short_tmp_part2/;
                                }
                                my $yuyu3_tmp3 = $yuyu3;
                                my $check_yuyuy3b = $yuyu3_tmp3 =~ s/.$end_short_tmp_part/$end_short_tmp_part/;
                                if ($check_yuyuy3b > 0 || $check_yuyuy3b2 > 0)
                                {
                                    push @filter_dot_before3, $yuyu0;
                                    if (exists($extensions_yuyu{$search}))
                                    {
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $yuyu0." FOUND3   ";
                                            print OUTPUT5 $extensions_yuyu{$search}." EXTe3\n";
                                        }
                                        push @extensions_before, $extensions_yuyu{$search};
                                        push @extensions_before3, $extensions_yuyu{$search};
                                    }
                                    $filter_before3{$search} = undef;
                                }
                            }
                        }
                    }
                    foreach my $search (keys %before_all4)
                    {
                        my $overhang_tmp = $overhang;
                        my $yuyu0 = $before_all4{$search};
                        my $yuyu4 = reverse($yuyu0);
                        if (length($yuyu4) >= $read_length-$overhang_tmp-5-1)
                        {
                            $yuyu4 =~ tr/N/\./;
                            my $yuyu4_tmp = $yuyu4;
                            my $check_yuyuy = $yuyu4_tmp =~ s/$end_short_tmp_part/$end_short_tmp_part/;
                            
                            if ($check_yuyuy > 0)
                            {
                                $fourth_yuyu2++;
                                $fourth_yuyu = "yes";
                                my $end_short_tmp_part = substr $end_short_tmp, -length($yuyu4)+1;
                                my $star = $end_short_tmp_part =~ tr/\*//;
                                my $check_yuyuy4b2 = '0';
                                my $end_short_tmp_part2;
                                if ($star > 0)
                                {
                                    $end_short_tmp_part = substr $end_short_tmp, -length($yuyu4)+1-($star*2);
                                    $end_short_tmp_part2 = substr $end_short_tmp, -length($yuyu4)+1-($star);
                                    my $yuyu4_tmp2 = $yuyu4;
                                    $check_yuyuy4b2 = $yuyu4_tmp2 =~ s/.$end_short_tmp_part2/$end_short_tmp_part2/;
                                }
                                my $yuyu4_tmp3 = $yuyu4;
                                my $check_yuyuy4b = $yuyu4_tmp3 =~ s/.$end_short_tmp_part/$end_short_tmp_part/;
                                if ($check_yuyuy4b > 0 || $check_yuyuy4b2 > 0)
                                {
                                    push @filter_dot_before4, $yuyu0;
                                    if (exists($extensions_yuyu{$search}))
                                    {
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $yuyu0." FOUND4   ";
                                            print OUTPUT5 $extensions_yuyu{$search}." EXTe4\n";
                                        }
                                        push @extensions_before, $extensions_yuyu{$search};
                                        push @extensions_before4, $extensions_yuyu{$search};
                                    }
                                    $filter_before4{$search} = undef;
                                }
                            }
                        }
                    }

                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $first_yuyu2." FIRST_YUYU\n";
                        print OUTPUT5 $second_yuyu2." SECOND_YUYU\n";
                        print OUTPUT5 $third_yuyu2." THIRD_YUYU\n";
                        print OUTPUT5 $fourth_yuyu2." FOURTH_YUYU\n";
                    }
                    my $correction = '0';
                    my $count_all = '0';
                    $count_all = (@filter_dot_before1) + (@filter_dot_before2) + (@filter_dot_before3) + (@filter_dot_before4);
                    my $end_short_tmp2 = substr $read_short_end2, -($read_length-$overhang-5-1);
                    my $check_dot = $end_short_tmp2 =~ tr/\./\./;
                    if ($check_dot eq "")
                    {
                        $check_dot = '0'
                    }
                    my $check_star = $end_short_tmp2 =~ tr/\*/\*/;
                    my $check_star2 = "";
                    my $deletion = "";
print OUTPUT5 $check_star." CHECK_STAR\n";
print OUTPUT5 $count_all." COUNT_ALL\n";
                    if ($count_all > 2 && ($check_dot > 0 || $check_star > 0))
                    {
                        my @split_dot = split /\./, reverse($end_short_tmp2);
                        my $size = '0';
                        my $size2 = '1';
                        my $length_total = '0';
                        my @filter_dot_before_all;
                        undef @filter_dot_before_all;
                        @filter_dot_before_all = (@filter_dot_before1, @filter_dot_before2, @filter_dot_before3, @filter_dot_before4);
                        
                        if ($check_star > 0)
                        {
                            my @deletion = split //, reverse($end_short_tmp2);
                            my $found_star = "";
                            foreach my $nuc (@deletion)
                            {
                                if ($found_star eq "yes")
                                {
                                    $found_star = "";
                                    $deletion .= $nuc;
                                }
                                if ($nuc eq "*")
                                {
                                    $found_star = "yes";
                                }
                            }
                            print OUTPUT5 $deletion." NUC0\n";
                        }
                        foreach my $dot_split (@split_dot)
                        {
                           $size++; 
                        }
                        foreach my $dot_split (@split_dot)
                        {
                            if ($size2 ne $size)
                            {
                                $length_total += length($dot_split);
                                my $A = '0';
                                my $C = '0';
                                my $T = '0';
                                my $G = '0';
                                my $check_star2_part = "";
                                if ($check_star > 0)
                                {
                                    $check_star2_part = substr reverse($end_short_tmp2), 0, $length_total;
                                    $check_star2 = $check_star2_part =~ tr/\*/\*/;
                                    print OUTPUT5 $check_star2." CHECK_STAR2\n";
                                }
                                foreach my $dot_before2 (@filter_dot_before_all)
                                {  
                                    my $dot = substr $dot_before2, $length_total, 1;
                                    if ($check_star2 > 0)
                                    {
                                        print OUTPUT5 $deletion." DEL_DOT\n";
                                        my $check_star2_part_tmp1 = $check_star2_part;
                                        $check_star2_part_tmp1 =~ tr/\*//d;
                                        my $dot_before2_tmp = $dot_before2;
                                        my $check_star1 = $dot_before2_tmp =~ s/$check_star2_part_tmp1/$check_star2_part_tmp1/;
                                        if ($check_star1 > 0)
                                        {
                                            $dot = substr $dot_before2, $length_total-$check_star2, 1;
                                        }
                                        else
                                        {
                                            $dot = substr $dot_before2, $length_total-($check_star2*2), 1;
                                        }
                                    }
                                    if ($dot eq "A")
                                    {
                                        $A++;
                                    }
                                    elsif ($dot eq "C")
                                    {
                                        $C++;
                                    }
                                    elsif ($dot eq "T")
                                    {
                                        $T++;
                                    }
                                    elsif ($dot eq "G")
                                    {
                                        $G++;
                                    } 
                                }
                                if (($A > 2 && ($C+$T+$G) eq 0) || ($A > ($C+$T+$G)*10 && $repetitive_detect eq ""))
                                {
                                    substr $read, -$length_total-1, 1, "A";
                                    print OUTPUT5 "DOT AAAAAAAAAAAA\n";
                                    $correction++;
                                }
                                if (($C > 2 && ($A+$T+$G) eq 0) || ($C > ($A+$T+$G)*10 && $repetitive_detect eq ""))
                                {
                                    substr $read, -$length_total-1, 1, "C";
                                    print OUTPUT5 "DOT CCCCCCCCCCCC\n";
                                    $correction++;
                                }
                                if (($T > 2 && ($C+$A+$G) eq 0) || ($T > ($C+$A+$G)*10 && $repetitive_detect eq ""))
                                {
                                    substr $read, -$length_total-1, 1, "T";
                                    print OUTPUT5 "DOT TTTTTTTTTTTT\n";
                                    $correction++;
                                }
                                if (($G > 2 && ($C+$T+$A) eq 0) || ($G > ($C+$T+$A)*10 && $repetitive_detect eq ""))
                                {
                                    substr $read, -$length_total-1, 1, "G";
                                    print OUTPUT5 "DOT GGGGGGGGGGGGG\n";
                                    $correction++;
                                }
                                $length_total += 1;
                            }
                            $size2++;
                        }
                    }
              
                    if ($first_yuyu eq "yes" && $second_yuyu ne "yes" && $third_yuyu ne "yes" && $fourth_yuyu ne "yes" && $first_yuyu2 > 3 && ($correction eq $check_dot || $overhang >= $read_length-($overlap+15)))
                    {
                        if (@extensions_before > 3000000)
                        {
                            $ext_before = "yes";
                            @extensions_before = @extensions_before1;
                            goto EXT_BEFORE;
                        }
                        delete $seed{$id_split1};
                        delete $seed{$id_split2};
                        delete $seed{$id_split3};
                        $best_extension = $best_extension1;
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nONLY FIRST\n\n";
                        }
                        delete $SNP_active{$id};
                        delete $before{$id};
                        delete $before_shorter_skip{$id};
                        goto AFTER_EXT;
                    }
                    elsif ($second_yuyu eq "yes" && $first_yuyu ne "yes" && $third_yuyu ne "yes" && $fourth_yuyu ne "yes" && $second_yuyu2 > 3 && ($correction eq $check_dot || $overhang >= $read_length-($overlap+15)))
                    {
                        if (@extensions_before > 3000000)
                        {
                            $ext_before = "yes";
                            @extensions_before = @extensions_before2;
                            goto EXT_BEFORE;
                        }
                        delete $seed{$id};
                        delete $seed{$id_split2};
                        delete $seed{$id_split3};
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nONLY SECOND\n\n";
                        }
                        delete $SNP_active{$id_split1};
                        delete $before{$id_split1};
                        delete $before_shorter_skip{$id_split1};
                        $id = $id_split1;
                        if ($correction > 0)
                        {
                            delete $seed{$id_split1};
                            $best_extension = $best_extension2;
                            goto AFTER_EXT;
                        }
                        else
                        {
                            $read_new = $seed{$id_split1};                       
                            $read_new1 = $read_new;
                            goto BACK;
                        }
                    }
                    elsif ($third_yuyu eq "yes" && $first_yuyu ne "yes" && $second_yuyu ne "yes" && $fourth_yuyu ne "yes" && $third_yuyu2 > 3 && ($correction eq $check_dot || $overhang >= $read_length-($overlap+15)))
                    {
                        if (@extensions_before > 300000)
                        {
                            $ext_before = "yes";
                            @extensions_before = @extensions_before3;
                            goto EXT_BEFORE;
                        }
                        delete $seed{$id};
                        delete $seed{$id_split3};
                        delete $seed{$id_split1};
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nONLY THIRD\n\n";
                        }
                        delete $SNP_active{$id_split2};
                        delete $before{$id_split2};
                        delete $before_shorter_skip{$id_split2};
                        $id = $id_split2;
                        if ($correction > 0)
                        {
                            delete $seed{$id_split2};
                            $best_extension = $best_extension3;
                            goto AFTER_EXT;
                        }
                        else
                        {
                            $read_new = $seed{$id_split2};                       
                            $read_new1 = $read_new;
                            goto BACK;
                        }
                    }
                    elsif ($fourth_yuyu eq "yes" && $first_yuyu ne "yes" && $third_yuyu ne "yes" && $second_yuyu ne "yes" && $fourth_yuyu2 > 3 && ($correction eq $check_dot || $overhang >= $read_length-($overlap+15)))
                    {
                        if (@extensions_before > 3000000)
                        {
                            $ext_before = "yes";
                            @extensions_before = @extensions_before4;
                            goto EXT_BEFORE;
                        }
                        delete $seed{$id};
                        delete $seed{$id_split2};
                        delete $seed{$id_split1};
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nONLY FOURTH\n\n";
                        }
                        delete $SNP_active{$id_split3};
                        delete $before{$id_split3};
                        delete $before_shorter_skip{$id_split3};
                        $id = $id_split3;
                        if ($correction > 0)
                        {
                            $best_extension = $best_extension4;
                            delete $seed{$id_split3};
                            goto AFTER_EXT;
                        }
                        else
                        {
                            $read_new = $seed{$id_split3};                       
                            $read_new1 = $read_new;
                            goto BACK;
                        }
                    }
                    elsif($count_all < 3 && $overhang < $read_length-($overlap+15) && $skip_overhang ne "yes")
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nMAKE BEFORE SHORTER\n\n";
                        }
                        $before_shorter = "yes";
                        goto BEFORE;   
                    }
                    else
                    {
                        my %count1234;
                        undef %count1234;
                        my %count1234b;
                        undef %count1234b;
                        
                        my $count1 = '0';
                        my $count2 = '0';
                        my $count3 = '0';
                        my $count4 = '0';
                        my $count1_pair = '0';
                        my $count2_pair = '0';
                        my $count3_pair = '0';
                        my $count4_pair = '0';
                        foreach my $exb (keys %filter_before1)
                        {
                            my $search_tmp = substr $exb, 0, -1;
                            my $search_end = substr $exb, -1;
                            if (exists($hash{$search_tmp}))
                            {
                                my @search_tmp = split /,/,$hash{$search_tmp};
                                my $found;
                                my $found2;
                                if ($search_end eq "1")
                                {
                                    $found = $search_tmp[1];
                                    $found2 = $search_tmp[0];
                                }
                                elsif ($search_end eq "2")
                                {
                                    $found = $search_tmp[0];
                                    $found2 = $search_tmp[1];
                                }
                                if ($encrypt eq "yes")
                                {
                                    $found = decrypt $found;
                                }
                                $count1++;
                                if (exists($filter_before1{$exb}))
                                {
                                    if (exists($filter_before1_pair{$exb}))
                                    {
                                        $filter_before1{$exb} = $found;
                                        $count1_pair++;
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $found." 1\n";
                                        }
                                    }
                                    else
                                    {
                                        delete $filter_before1{$exb};
                                    }
                                }
                            }
                        }
                        foreach my $exb (keys %filter_before2)
                        {
                            my $search_tmp = substr $exb, 0, -1;
                            my $search_end = substr $exb, -1;
                            if (exists($hash{$search_tmp}))
                            {
                                my @search_tmp = split /,/,$hash{$search_tmp};
                                my $found;
                                my $found2;
                                if ($search_end eq "1")
                                {
                                    $found = $search_tmp[1];
                                    $found2 = $search_tmp[0];
                                }
                                elsif ($search_end eq "2")
                                {
                                    $found = $search_tmp[0];
                                    $found2 = $search_tmp[1];
                                }
                                if ($encrypt eq "yes")
                                {
                                    $found = decrypt $found;
                                }
                                $count2++;
                                if (exists($filter_before2{$exb}))
                                {
                                    if (exists($filter_before2_pair{$exb}))
                                    {
                                        $filter_before2{$exb} = $found;
                                        $count2_pair++;
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $found." 2\n";
                                        }
                                    }
                                    else
                                    {
                                        delete $filter_before2{$exb};
                                    }
                                }
                            }
                        }
                        foreach my $exb (keys %filter_before3)
                        {
                            my $search_tmp = substr $exb, 0, -1;
                            my $search_end = substr $exb, -1;
                            if (exists($hash{$search_tmp}))
                            {
                                my @search_tmp = split /,/,$hash{$search_tmp};
                                my $found;
                                if ($search_end eq "1")
                                {
                                    $found = $search_tmp[1];
                                }
                                elsif ($search_end eq "2")
                                {
                                    $found = $search_tmp[0];
                                }
                                if ($encrypt eq "yes")
                                {
                                    $found = decrypt $found;
                                }
                                $count3++;
                                if (exists($filter_before3{$exb}))
                                {
                                    if (exists($filter_before3_pair{$exb}))
                                    {
                                        $filter_before3{$exb} = $found;
                                        $count3_pair++;
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $found." 3\n";
                                        }
                                    }
                                    else
                                    {
                                        delete $filter_before3{$exb};
                                    }
                                }
                            }
                        }
                        foreach my $exb (keys %filter_before4)
                        {
                            my $search_tmp = substr $exb, 0, -1;
                            my $search_end = substr $exb, -1;
                            if (exists($hash{$search_tmp}))
                            {
                                my @search_tmp = split /,/,$hash{$search_tmp};
                                my $found;
                                if ($search_end eq "1")
                                {
                                    $found = $search_tmp[1];
                                }
                                elsif ($search_end eq "2")
                                {
                                    $found = $search_tmp[0];
                                }
                                if ($encrypt eq "yes")
                                {
                                    $found = decrypt $found;
                                }                               
                                $count4++;
                                if (exists($filter_before4{$exb}))
                                {
                                    if (exists($filter_before4_pair{$exb}))
                                    {
                                        $filter_before4{$exb} = $found;
                                        $count4_pair++;
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $found." 4\n";
                                        }
                                    }
                                    else
                                    {
                                        delete $filter_before4{$exb};
                                    }
                                }
                            }
                        }
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 $count1." COUNT1\n";
                            print OUTPUT5 $count2." COUNT2\n";
                            print OUTPUT5 $count3." COUNT3\n";
                            print OUTPUT5 $count4." COUNT4\n";
                        }
                        $count1234{'1'} = $count1_pair;
                        $count1234{'2'} = $count2_pair;
                        $count1234{'3'} = $count3_pair;
                        $count1234{'4'} = $count4_pair;
                        
                        my $morethan3 = '0';
                        my $difference = "";
                        
                        foreach my $count1234 (keys %count1234)
                        {
                            if ($count1234{$count1234} > 3)
                            {
                                $morethan3++;
                            }
                        }
                        my $best_extension1_tmp = $best_extension1;
                        my $best_extension2_tmp = $best_extension2;
                        my $SNR_check1 = $best_extension1_tmp =~ s/AAAAAAAA|CCCCCCCC|GGGGGGGG|TTTTTTTT//;
                        my $SNR_check2 = $best_extension2_tmp =~ s/AAAAAAAA|CCCCCCCC|GGGGGGGG|TTTTTTTT//;
                        my $GGGG2 = $best_extension2 =~ tr/G/G/;
                        my $TTTT2 = $best_extension2 =~ tr/T/T/;
                        my $CCCC2 = $best_extension2 =~ tr/C/C/;
                        my $AAAA2 = $best_extension2 =~ tr/A/A/;
                        if ($GGGG2 eq length($best_extension2) || $TTTT2 eq length($best_extension2) || $CCCC2 eq length($best_extension2) || $AAAA2 eq length($best_extension2))
                        {
                            $SNR_check2 = '1';
                        }
                        my $GGGG1 = $best_extension1 =~ tr/G/G/;
                        my $TTTT1 = $best_extension1 =~ tr/T/T/;
                        my $CCCC1 = $best_extension1 =~ tr/C/C/;
                        my $AAAA1 = $best_extension1 =~ tr/A/A/;
                        if ($GGGG1 eq length($best_extension1) || $TTTT1 eq length($best_extension1) || $CCCC1 eq length($best_extension1) || $AAAA1 eq length($best_extension1))
                        {
                            $SNR_check1 = '1';
                        }

                        if ((($count1 > 3 && ($count2+$count3+$count4) eq 0) || ($count1 > ($count2+$count3+$count4)*6 && ($count2+$count3+$count4) ne 0)) && $repetitive_detect2 ne "yes" && $correction eq $check_dot)
                        {
                            delete $seed{$id_split1};
                            delete $seed{$id_split2};
                            delete $seed{$id_split3};
                            $best_extension = $best_extension1;
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE9\n\n";
                            }
                            delete $SNP_active{$id};
                            delete $before{$id};
                            delete $before_shorter_skip{$id};
                            goto AFTER_EXT;
                        }
                        elsif ((($count2 > 3 && ($count1+$count3+$count4) eq 0) || ($count2 > ($count1+$count3+$count4)*6 && ($count1+$count3+$count4) ne 0)) && $repetitive_detect2 ne "yes" && $correction eq $check_dot)
                        {
                            delete $seed{$id};
                            delete $seed{$id_split1};
                            delete $seed{$id_split2};
                            delete $seed{$id_split3};
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE10\n\n";
                            }
                            delete $SNP_active{$id_split1};
                            delete $before{$id_split1};
                            delete $before_shorter_skip{$id_split1};
                            $id = $id_split1;

                            $best_extension = $best_extension2;
                            goto AFTER_EXT;
                        }
                        elsif ((($count3 > 3 && ($count1+$count2+$count4) eq 0) || ($count3 > ($count1+$count2+$count4)*6 && ($count1+$count2+$count4) ne 0)) && $repetitive_detect2 ne "yes" && $correction eq $check_dot)
                        {
                            delete $seed{$id};
                            delete $seed{$id_split3};
                            delete $seed{$id_split1};
                            delete $seed{$id_split2};
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE11\n\n";
                            }
                            delete $SNP_active{$id_split2};
                            delete $before{$id_split2};
                            delete $before_shorter_skip{$id_split2};
                            $id = $id_split2;
                            $best_extension = $best_extension3;
                            goto AFTER_EXT;
                        }
                        elsif ((($count4 > 3 && ($count1+$count3+$count2) eq 0) || ($count4 > ($count1+$count3+$count2)*6 && ($count1+$count3+$count2) ne 0)) && $repetitive_detect2 ne "yes" && $correction eq $check_dot)
                        {
                            delete $seed{$id};
                            delete $seed{$id_split2};
                            delete $seed{$id_split1};
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE12\n\n";
                            }
                            delete $SNP_active{$id_split3};
                            delete $before{$id_split3};
                            delete $before_shorter_skip{$id_split3};
                            $id = $id_split3;

                            $best_extension = $best_extension4;
                            goto AFTER_EXT;
                        }
                        elsif((($morethan3 < 2 && $overhang < $read_length-($overlap*1.5)) || ($morethan3 > 1 && length($read) <= $insert_size && $overhang < 30)) && $skip_overhang ne "yes")
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nMAKE BEFORE SHORTER2\n\n";
                            }
                            $before_shorter = "yes";
                            goto BEFORE;
                        }
                        elsif ($morethan3 > 1 && $last_chance ne "yes" && (($SNR_check1 eq "" && $SNR_check2 eq "") || $overhang < 30))
                        {
                            my $count1b = '0';
                            my $count2b = '0';
                            my $count3b = '0';
                            my $count4b = '0';
                                                                                                   
                            my $size = keys %read_short_end_tmp;
                            undef @extensions_before;
                            undef @extensions_before1;
                            undef @extensions_before2;
                            undef @extensions_before3;
                            undef @extensions_before4;
                            
                            my $insert_range_before = 1.65;
                            if ($insert_range_shorter eq "yes")
                            {
                                $insert_range_before = 1.45;
                            }

                            my $read_short_end_tempie = substr $read, -($insert_size-7-($read_length/2))-(($insert_size*$insert_range_before)-$insert_size)-(($read_length-$right-$left)/2), ((($insert_size*$insert_range_before)-$insert_size)*2)+($read_length-$right-$left);
                            $read_short_end_tempie =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                               
                            undef %read_short_end_tmp;

                            %read_short_end_tmp = build_partial3c $read_short_end_tempie;
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 $read_short_end_tempie." READ_SHORT\n";
                            }
                        
                           my $ff = '0';
                            my %hash_read_short_end;
                            undef %hash_read_short_end;
                            foreach my $read_short_end_tempie (keys %read_short_end_tmp)
                            {
                                $ff = '0';
                                while ($ff < (length($read_short_end_tempie)-$read_length+$left+$right+2))
                                {
                                    my $read_short_end_part = substr $read_short_end_tempie, $ff, $read_length-$left-$right-1;
                                                                   
                                    $hash_read_short_end{$read_short_end_part} = undef;
                                    $ff++;
                                }
                            }

                            foreach my $exb0 (keys %filter_before1)
                            { 
                                my $exb = $filter_before1{$exb0};
                                my $match_pair_middle = substr $exb, $left, $read_length-$left-$right;
FILTER_1:                       foreach my $line (keys %hash_read_short_end)
                                {
                                    my $found_seq = '0';
                                    my $match_pair_middle_tmp = $match_pair_middle;

                                    $found_seq = $match_pair_middle_tmp =~ s/.$line/$line/;
                                    if ($found_seq > 0)
                                    {
                                        if (exists($filter_before1_pair{$exb0}))
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $filter_before1_pair{$exb0}." FOUND1\n";
                                                print OUTPUT5 $exb." FOUND1_PAIR\n";
                                            }
                                            if (exists($extensions_yuyu{$exb0}))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 $extensions_yuyu{$exb0}." EXT1\n";
                                                }
                                                push @extensions_before, $extensions_yuyu{$exb0};
                                                push @extensions_before1, $extensions_yuyu{$exb0};
                                            }
                                        }
                                        $count1b++;
                                        last FILTER_1;
                                    }
                                }
                            }
                            foreach my $exb0 (keys %filter_before2)
                            { 
                                my $exb = $filter_before2{$exb0};
                                my $match_pair_middle = substr $exb, $left, $read_length-$left-$right;
FILTER_2:                       foreach my $line (keys %hash_read_short_end)
                                {
                                    my $found_seq = '0';
                                    my $match_pair_middle_tmp = $match_pair_middle;

                                    $found_seq = $match_pair_middle_tmp =~ s/.$line/$line/;
                                    if ($found_seq > 0)
                                    {
                                        if (exists($filter_before2_pair{$exb0}))
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $filter_before2_pair{$exb0}." FOUND2\n";
                                                print OUTPUT5 $exb." FOUND2_PAIR\n";
                                            }
                                            if (exists($extensions_yuyu{$exb0}))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 $extensions_yuyu{$exb0}." EXT2\n";
                                                }
                                                push @extensions_before, $extensions_yuyu{$exb0};
                                                push @extensions_before2, $extensions_yuyu{$exb0};
                                            }
                                        }
                                        $count2b++;
                                        last FILTER_2;
                                    }
                                }
                            }
                            foreach my $exb0 (keys %filter_before3)
                            { 
                                my $exb = $filter_before3{$exb0};
                                my $match_pair_middle = substr $exb, $left, $read_length-$left-$right;
FILTER_3:                       foreach my $line (keys %hash_read_short_end)
                                {
                                    my $found_seq = '0';
                                    my $match_pair_middle_tmp = $match_pair_middle;

                                    $found_seq = $match_pair_middle_tmp =~ s/.$line/$line/;
                                    if ($found_seq > 0)
                                    {
                                        if (exists($filter_before3_pair{$exb0}))
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $filter_before3_pair{$exb0}." FOUND3\n";
                                                print OUTPUT5 $exb." FOUND3_PAIR\n";
                                            }
                                            if (exists($extensions_yuyu{$exb0}))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 $extensions_yuyu{$exb0}." EXT3\n";
                                                }
                                                push @extensions_before, $extensions_yuyu{$exb0};
                                                push @extensions_before3, $extensions_yuyu{$exb0}
                                            }
                                        }
                                        $count3b++;
                                        last FILTER_3;
                                    }
                                }
                            }
                            foreach my $exb0 (keys %filter_before4)
                            { 
                                my $exb = $filter_before4{$exb0};
                                my $match_pair_middle = substr $exb, $left, $read_length-$left-$right;
FILTER_4:                       foreach my $line (keys %hash_read_short_end)
                                {
                                    my $found_seq = '0';
                                    my $match_pair_middle_tmp = $match_pair_middle;

                                    $found_seq = $match_pair_middle_tmp =~ s/.$line/$line/;
                                    if ($found_seq > 0)
                                    {
                                        if (exists($filter_before4_pair{$exb0}))
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $filter_before4_pair{$exb0}." FOUND4\n";
                                                print OUTPUT5 $exb." FOUND4_PAIR\n";
                                            }
                                            if (exists($extensions_yuyu{$exb0}))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 $extensions_yuyu{$exb0}." EXT4\n";
                                                }
                                                push @extensions_before, $extensions_yuyu{$exb0};
                                                push @extensions_before4, $extensions_yuyu{$exb0};                           
                                            }
                                        }
                                        $count4b++;
                                        last FILTER_4;
                                    }
                                }
                            }
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 $count1b." COUNT1B\n";
                                print OUTPUT5 $count2b." COUNT2B\n";
                                print OUTPUT5 $count3b." COUNT3B\n";
                                print OUTPUT5 $count4b." COUNT4B\n";
                            }
                            my $f = '5';
                            if ($repetitive_detect ne "")
                            {
                                $f = '10';
                            }
                            
                            $count1234b{'1'} = $count1b;
                            $count1234b{'2'} = $count2b;
                            $count1234b{'3'} = $count3b;
                            $count1234b{'4'} = $count4b;
                            
                            my $differenceb1 = "";
                            my $differenceb2 = "";
                            my $differenceb3 = "";
                            my $differenceb4 = "";
                            
                            foreach my $count1234b (keys %count1234b)
                            {   
                                if ($count1234b ne '1' && $count1234b{$count1234b} > 0 && $count1b > $f*$count1234b{$count1234b})
                                {
                                    $differenceb1 = "yes";
                                }
                                if ($count1234b ne '2' && $count1234b{$count1234b} > 0 && $count2b > $f*$count1234b{$count1234b})
                                {
                                    $differenceb2 = "yes";
                                }
                                if ($count1234b ne '3' && $count1234b{$count1234b} > 0 && $count3b > $f*$count1234b{$count1234b})
                                {
                                    $differenceb3 = "yes";
                                }
                                if ($count1234b ne '4' && $count1234b{$count1234b} > 0 && $count4b > $f*$count1234b{$count1234b})
                                {
                                    $differenceb4 = "yes";
                                }
                            }
                            
                            if ((($count1b > 2 && ($count2b+$count3b+$count4b) eq '0') || $differenceb1 eq "yes"))
                            {
                                if (@extensions_before > 3000000)
                                {
                                    $ext_before = "yes";
                                    @extensions_before = @extensions_before1;
                                    goto EXT_BEFORE;
                                }
                                delete $seed{$id_split1};
                                delete $seed{$id_split2};
                                delete $seed{$id_split3};
                                $best_extension = $best_extension1;
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE5\n\n";
                                }
                                delete $SNP_active{$id};
                                delete $before{$id};
                                delete $before_shorter_skip{$id};
                                goto AFTER_EXT;
                            }
                            elsif ((($count2b > 2 && ($count1b+$count3b+$count4b) eq '0') || $differenceb2 eq "yes"))
                            {
                                if (@extensions_before > 3000000)
                                {
                                    $ext_before = "yes";
                                    @extensions_before = @extensions_before2;
                                    goto EXT_BEFORE;
                                }
                                delete $seed{$id};
                                delete $seed{$id_split2};
                                delete $seed{$id_split3};
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE6\n\n";
                                }
                                delete $SNP_active{$id_split1};
                                delete $before_shorter_skip{$id_split1};
                                delete $before{$id_split1};
                                $id = $id_split1;
                                if ($correction > 0)
                                {
                                    delete $seed{$id_split1};
                                    $best_extension = $best_extension2;
                                    goto AFTER_EXT;
                                }
                                else
                                {
                                    $read_new = $seed{$id_split1};                       
                                    $read_new1 = $read_new;
                                    goto BACK;
                                }
                            }
                            elsif ((($count3b > 2 && ($count1b+$count2b+$count4b) eq '0') || $differenceb3 eq "yes"))
                            {
                                if (@extensions_before > 300000)
                                {
                                    $ext_before = "yes";
                                    @extensions_before = @extensions_before3;
                                    goto EXT_BEFORE;
                                }
                                delete $seed{$id};
                                delete $seed{$id_split1};
                                delete $seed{$id_split3};
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE7\n\n";
                                }
                                delete $SNP_active{$id_split2};
                                delete $before{$id_split2};
                                delete $before_shorter_skip{$id_split2};
                                $id = $id_split2;
                                if ($correction > 0)
                                {
                                    delete $seed{$id_split2};
                                    $best_extension = $best_extension3;
                                    goto AFTER_EXT;
                                }
                                else
                                {
                                    $read_new = $seed{$id_split2};                       
                                    $read_new1 = $read_new;
                                    goto BACK;
                                }
                            }
                            elsif ((($count4b > 2 && ($count1b+$count2b+$count3b) eq '0') || $differenceb4 eq "yes"))
                            {
                                if (@extensions_before > 300000)
                                {
                                    $ext_before = "yes";
                                    @extensions_before = @extensions_before4;
                                    goto EXT_BEFORE;
                                }
                                delete $seed{$id};
                                delete $seed{$id_split1};
                                delete $seed{$id_split2};
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE8\n\n";
                                }
                                delete $SNP_active{$id_split3};
                                delete $before{$id_split3};
                                delete $before_shorter_skip{$id_split3};
                                $id = $id_split3;
                                if ($correction > 0)
                                {
                                    delete $seed{$id_split3};
                                    $best_extension = $best_extension4;
                                    goto AFTER_EXT;
                                }
                                else
                                {
                                    $read_new = $seed{$id_split3};                       
                                    $read_new1 = $read_new;
                                    goto BACK;
                                }
                            }
                            elsif (($count1b+$count2b+$count3b+$count4b) < 10 && $overhang < $read_length-($overlap*1.5) && $before_shorter_skip ne "yes" && $skip_overhang ne "yes")
                            {
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nMAKE BEFORE SHORTER3\n\n";
                                }
                                $before_shorter = "yes";
                                goto BEFORE;
                            }
                            elsif(($count1b+$count2b+$count3b+$count4b) < 10 && $overhang < $read_length-($overlap*1.5) && $before_shorter_skip ne "yes" && $insert_range_shorter eq "yes")
                            {
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nMAKE BEFORE SHORTER3a\n\n";
                                }
                                $insert_range_shorter = "";
                                goto BEFORE;
                            }
                            elsif(($count1b+$count2b+$count3b+$count4b) > 9 && $before_shorter_skip ne "yes" && $insert_range_shorter eq "")
                            {
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nSMALLER INSERT\n\n";
                                }
                                $insert_range_shorter = "yes";
                                goto BEFORE;
                            }
                            elsif (($count1b+$count2b+$count3b+$count4b) > 9 && $extensions_before ne "yes" && $repetitive_detect2 ne "yes" && $before_shorter_skip ne "yes")
                            {
                                $l = 0;
                                $best_extension = "";
                                $SNP = "";
                                @extensions = @extensions_before;
                                $ext = '0';
                                foreach (@extensions)
                                {
                                    $ext++;
                                }
                                $extensions_before = "yes";
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nEXTENSIONS FROM BEFORE\n\n";
                                }
                                $SNP_active = "yes";
                                delete $before{$id};
                                $before_shorter_skip{$id} = "yes";
                                $split_forward = "";
                                delete $seed{$id};
                                delete $seed{$id_split1};
                                delete $seed{$id_split2};
                                delete $seed{$id_split3};
                                undef @extensions_group1;
                                undef @extensions_group2;
                                undef @extensions_group3;
                                undef @extensions_group4;
                                goto NUCLEO;
                            }          
                        }
EXT_BEFORE:             if ($ext_before eq "yes")
                        {
                            $l = 0;
                            $best_extension = "";
                            $SNP = "";
                            @extensions = @extensions_before;
                            $ext = '0';
                            foreach (@extensions)
                            {
                                $ext++;
                            }
                            $extensions_before = "yes";
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nEXTENSIONS FROM BEFORE\n\n";
                            }
                            $SNP_active = "yes";
                            delete $before{$id};
                            $before_shorter_skip{$id} = "yes";
                            $split_forward = "";
                            delete $seed{$id};
                            delete $seed{$id_split1};
                            delete $seed{$id_split2};
                            delete $seed{$id_split3};
                            undef @extensions_group1;
                            undef @extensions_group2;
                            undef @extensions_group3;
                            undef @extensions_group4;
                            goto NUCLEO;
                        }    
                                if ($repetitive_detect eq "yes" && $jump_rep eq "yes")
                                {
                                    delete $jump_rep{$id};
                                    my $end_repetitive = substr $read, -$insert_size-150;
                                    if ($y > $startprint2)
                                    {
                                        print OUTPUT5 $end_repetitive." END_REP\n";
                                    }
                                    my $second_try = "";
                                    my %rep_pair;
                                    my %rep_pair_exclude;
REP_PAIR0:                          undef %rep_pair;
                                    undef %rep_pair_exclude;
REP_PAIR:                           foreach my $rep_pair (keys %repetitive_pair)
                                    {
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $rep_pair." REP_PAIR_TEST\n";
                                        }
                                        my $part1 = substr $rep_pair, 0, ($read_length/3)*2;
                                        my $part2 = substr $rep_pair, ($read_length/3)*2-5;
                                        my $r = '0';
                                        while ($r < length($part1)-15)
                                        {
                                            my $testit = substr $part1, $r, 15;
                                            my $found8 = $end_repetitive =~ s/$testit/$testit/;
                                            
                                            if ($found8 > 0 || $second_try eq "yes2" || $second_try eq "yes3")
                                            {
                                                my $s = '0';
                                                while ($s < length($part2)-15)
                                                {
                                                    my $testit2 = substr $part2, $s, 15;
                                                    my $found7 = $end_repetitive =~ s/$testit2/$testit2/;
                                                    my $found9 = $part1 =~ s/$testit2/$testit2/;
                                                    if ($found7 > 0 || $found9 > 0)
                                                    {
                                                        next REP_PAIR;
                                                    }
                                                    $s = $s+2;
                                                }
                                                $rep_pair{$rep_pair} = $repetitive_pair{$rep_pair};
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 $rep_pair." REP_PAIR_FOUND\n";
                                                }
                                                next REP_PAIR;
                                            }
                                            $r = $r+5;
                                        }
                                    }
                                    my $hg = '0';
                                    foreach (keys %rep_pair)
                                    {
                                        $hg++;
                                    }
                                    my %temp_rep;
                                    undef %temp_rep;
                                    my $search_rev;
                                    if ($hg eq '0')
                                    {
PAIR_OF_PAIR:                           foreach my $rep_pair (keys %repetitive_pair)
                                        {
                                            my $ft = '0';
                                            while ($ft < length($rep_pair) - $overlap)
                                            {
                                                my $part_rep = substr $rep_pair, $ft, $overlap;
                                                if (exists($hash2b{$part_rep}))
                                                {                       
                                                    my $search = $hash2b{$part_rep};
                                                    $search = substr $search, 1;
                                                    my @search = split /,/,$search;
                                                                                                
                                                    foreach my $search (@search)
                                                    {                             
                                                        my $search_tmp = substr $search, 0, -1;
                                                        my $search_end = substr $search, -1;
                                                        if (exists($hash{$search_tmp}))
                                                        {
                                                            my @search_tmp = split /,/,$hash{$search_tmp};
                                                            my $found;
                                                            if ($search_end eq "1")
                                                            {
                                                                $found = $search_tmp[1];
                                                                $search_rev = $search_tmp."2";
                                                            }
                                                            elsif ($search_end eq "2")
                                                            {
                                                                $found = $search_tmp[0];
                                                                $search_rev = $search_tmp."1";
                                                            }
                                                            if ($encrypt eq "yes")
                                                            {
                                                                $found = decrypt $found;
                                                            }
                                                            $found =~ tr/ACTG/TGAC/;
                                                            my $found_reverse = reverse($found);
                                                            delete $repetitive_pair{$rep_pair};
                                                            $temp_rep{$found_reverse} = $search_rev;
                                                            next PAIR_OF_PAIR;
                                                        }
                                                    }
                                                }
                                                $ft++;
                                            }
                                            delete $repetitive_pair{$rep_pair};
                                        }
                                        %repetitive_pair = %temp_rep;
                                        my $mm = '0';
                                        foreach (keys %repetitive_pair)
                                        {
                                            $mm++;
                                        }
                                        if ($mm > 0)
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 "\n"
                                            }
                                            if ($second_try eq "yes")
                                            {
                                                $second_try = "yes2";
                                            }
                                            elsif ($second_try eq "")
                                            {
                                                $second_try = "yes";
                                            }
                                            elsif ($second_try eq "yes2")
                                            {
                                                 $second_try = "yes3";
                                            }
                                            elsif ($second_try eq "yes3")
                                            {
                                                 goto REP_PAIR1;
                                            }
                                            goto REP_PAIR0;
                                        }
                                    }
                                    my $most_match_max = '0';
                                    my $id_rep;
                                    my $rep_pair2;
REP_PAIR1:                          foreach my $rep_pair (keys %rep_pair)
                                    {
                                        if (exists($rep_pair_exclude{$rep_pair}))
                                        {
                                            next REP_PAIR1;
                                        }
                                        my $part2 = substr $rep_pair, -$overlap-10;
                                        my $s = '0';
                                        my $most_match_total = '0';
                                        while ($s < length($part2)-$overlap)
                                        {
                                            my $part2b = substr $part2, $s, $overlap;
                                            my $part2b_reverse = reverse($part2b);
                                            $part2b_reverse =~ tr/ACTG/TGAC/;
                                            if (exists($hash2b{$part2b}))
                                            {
                                                $most_match_total++;
                                            }
                                            if (exists($hash2c{$part2b}))
                                            {
                                                $most_match_total++;
                                            }
                                            if (exists($hash2b{$part2b_reverse}))
                                            {
                                                $most_match_total++;
                                            }
                                            if (exists($hash2c{$part2b_reverse}))
                                            {
                                                $most_match_total++;
                                            }
                                            $s++;
                                        }
                                        if ($most_match_total > $most_match_max)
                                        {
                                            $most_match_max = $most_match_total;
                                            $id_rep = $rep_pair{$rep_pair};
                                            while (exists($seeds_check{$id_rep}))
                                            {
                                                $id_rep = "0".$id_rep;
                                            }
                                            foreach my $test (keys %contigs)                                                
                                            {
                                                if ($test =~ m/.*\+$id_rep$/)
                                                {
                                                    $id_rep = "0".$id_rep;
                                                }
                                                if ($test =~ m/.*\+0$id_rep$/)
                                                {
                                                    $id_rep = "00".$id_rep;
                                                }
                                            }
                                            $rep_pair2 = $rep_pair;
                                        }      
                                    }
                                    if ($id_rep ne "")
                                    {
                                        $noforward{$id} = "stop";
                                        $noforward = "stop";
                                        my $rep_pair2_tmp = $rep_pair2;
                                        if ($rep_pair2 eq "")
                                        {
                                            $rep_pair_exclude{$rep_pair2_tmp} = undef;
                                            $id_rep = "";
                                            goto REP_PAIR1;
                                        }
                                        $seed{$id_rep} = $rep_pair2;
                                        $seeds_check{$id_rep} = undef;
                                        $nosecond{$id} = undef;
                                        $nosecond = "yes";
                                        $insert_size2{$id_rep} = $insert_size;
                                        $position{$id_rep} = length($rep_pair2);
                                        $old_id{$id_rep} = $id;
                                        $old_id2{$id_rep} = undef;
                                        $old_rep{$id_rep} = undef;
                                        $old_rep_old{$id} = undef;
                                        $seed_old{$id} = $read;
                                        
                                        $read_new = $read;     
                                        $read_new1 = $read_new;
                                        delete $seed{$id_split1};
                                        
                                        if (exists($old_id{$id}))
                                        {
                                            $merge_now = "yes";
                                            goto MERGE;
                                        }                                       
                                        if ($noback eq "stop")
                                        {  
                                            print OUTPUT5 ">".$id."\n";
                                            print OUTPUT5 $read."\n";
                                            delete $seed{$id};
                                            if (!keys %seed)
                                            {
                                                $circle = "contigs";
                                                goto FINISH;
                                            }
                                            else
                                            {
                                                goto ITERATION;
                                            }
                                        }
                                        else
                                        {
                                            print OUTPUT5 ">\n";
                                            $noforward = "stop";
                                            $noforward{$id} = "stop";
                                            $seed_split{$id} = undef;
                                            $best_extension = "";
                                            goto BACK;
                                        }
                                    }

                                    $end_repetitive =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                    my $test_dot = $end_repetitive =~ tr/\./\./;
                                    my %end_repetitive;
                                    undef %end_repetitive;
                                    my %rep_check;
                                    undef %rep_check;
                                    my %extensions_rep;
                                    undef %extensions_rep;
                                    if ($test_dot > 0)
                                    {
                                        %end_repetitive = build_partial3b $end_repetitive;
                                    }
                                    else
                                    {
                                        $end_repetitive{$end_repetitive} = undef;
                                    }
                                    my $read_end_shortie =  substr $read_end, -15;
REP_CHECK0:                         foreach my $exts (keys %extensions_original)
                                    {
                                        if (length($exts) > 25)
                                        {
                                            my $part1 = substr $exts, 0, length($exts)/2;
                                            my $part2 = substr $exts, length($exts)/2;
                                            my $r = '0';
                                            
                                            while ($r < length($part1)-11)
                                            {
                                                my $testit = substr $part1, $r, 11;
                                                foreach my $end_repetitive_tmp (keys %end_repetitive)
                                                {
                                                    my $found8 = $end_repetitive_tmp =~ s/$testit/$testit/;
                                                    if ($found8 > 0)
                                                    {
                                                        my $s = '0';
                                                        
                                                        while ($s < length($part2)-11)
                                                        {
                                                            my $testit2 = substr $part2, $s, 11;
                                                            my $found7 = $end_repetitive_tmp =~ s/$testit2/$testit2/;
                                                            if ($found7 > 0)
                                                            {
                                                                next REP_CHECK0;
                                                            }
                                                            $s++;
                                                        }
                                                        my $search_tmp = substr $extensions_original{$exts}, 0, -1;
                                                        my $search_end = substr $extensions_original{$exts}, -1;
                                                        if (exists($hash{$search_tmp}))
                                                        {
                                                            my @search_tmp = split /,/,$hash{$search_tmp};
                                                            my $found;
                                                            if ($search_end eq "1")
                                                            {
                                                                $found = $search_tmp[0];
                                                            }
                                                            elsif ($search_end eq "2")
                                                            {
                                                                $found = $search_tmp[1];
                                                            }
                                                            if ($encrypt eq "yes")
                                                            {
                                                                $found = decrypt $found;
                                                            }
                                                            $found =~ tr/ACTG/TGAC/;
                                                            my $found2 = reverse($found);

                                                            my $checkedie = $found2 =~ s/$exts/$exts/;
                                                            if ($checkedie > 0)
                                                            {
                                                                $rep_check{$found2} = $extensions_original{$exts};
                                                                if ($y > $startprint2)
                                                                {
                                                                    print OUTPUT5 $found2." REP_CHECK_FULL\n";
                                                                    print OUTPUT5 $extensions_original{$exts}." REP_CHECK_FULL2\n";
                                                                }
                                                                next REP_CHECK0;
                                                            }
                                                        }  
                                                    }
                                                }
                                                $r++;
                                            }
                                        }
                                    }
                                    foreach my $rep_pair (keys %rep_check)
                                    {
                                        $noforward{$id} = "stop";
                                        $noforward = "stop";
                                        my $id_rep = $rep_check{$rep_pair};
                                        $seed{$id_rep} = $rep_pair;
                                        $seeds_check{$id_rep} = undef;
                                        $nosecond{$id} = undef;
                                        $nosecond = "yes";
                                        $insert_size2{$id_rep} = $insert_size;
                                        $position{$id_rep} = length($rep_pair);
                                        $old_id{$id_rep} = $id;
                                        $old_id2{$id_rep} = undef;
                                        $old_rep{$id_rep} = undef;
                                        $old_rep_old{$id} = undef;
                                        $seed_old{$id} = $read;
                                        
                                        $read_new = $read;     
                                        $read_new1 = $read_new;
                                        delete $seed{$id_split1};
                                        
                                        if (exists($old_id{$id}))
                                        {
                                            $merge_now = "yes";
                                            goto MERGE;
                                        }
                                        if ($noback eq "stop")
                                        {
                                            print OUTPUT5 ">".$id."\n";
                                            print OUTPUT5 $read."\n";
                                            
                                            delete $seed{$id};
                                            if (!keys %seed)
                                            {
                                                $circle = "contigs";
                                                goto FINISH;
                                            }
                                            else
                                            {
                                                goto ITERATION;
                                            }
                                        }
                                        else
                                        {
                                            print OUTPUT5 ">\n";
                                            $noforward = "stop";
                                            $noforward{$id} = "stop";
                                            $seed_split{$id} = undef;
                                            $best_extension = "";
                                            goto BACK;
                                        }
                                    }
                                }
                                
                                delete $seed{$id_split1};
                                delete $seed{$id_split2};
                                delete $seed{$id_split3};
                                $before_shorter_skip{$id} = "yes";
                                if ($before_shorter_skip eq "yes")
                                {
                                    $before{$id} = "yes";
                                }
                                $best_extension = "";
                                $read_new = $read;     
                                $read_new1 = $read_new;
                                if ($count_split > 2)
                                {
                                    my $tmp = '0';
                                    if ($first_yuyu ne "yes")
                                    {
                                        $yuyu_option{$id.$firstSNP_max[$tmp]} = "yes";
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $firstSNP_max[$tmp]." FIRST_YUYU\n";
                                        }
                                    }
                                    $tmp++;
                                    if ($second_yuyu ne "yes")
                                    {
                                        $yuyu_option{$id.$firstSNP_max[$tmp]} = "yes";
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $firstSNP_max[$tmp]." SECOND_YUYU\n";
                                        }
                                    }
                                    $tmp++;
                                    if ($third_yuyu ne "yes")
                                    {
                                        $yuyu_option{$id.$firstSNP_max[$tmp]} = "yes";
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $firstSNP_max[$tmp]." THIRD_YUYU\n";
                                        }
                                    }
                                    $tmp++;
                                    if ($fourth_yuyu ne "yes")
                                    {
                                        $yuyu_option{$id.$firstSNP_max[$tmp]} = "yes";
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $firstSNP_max[$tmp]." FOURTH_YUYU\n";
                                        }
                                    }
                                }
                                goto BACK;
                    }
                }
                my $split_no_rep = "";
                if ($jump_rep ne "yes" && $repetitive_detect eq "yes")
                {
                    my $read_end9 = substr $read_short_end2, -$read_length;
                    my $ext1 = substr $best_extension1, 0, 9;
                    my $ext2 = substr $best_extension2, 0, 9;
                    my $check_ext1 = $read_end9 =~ s/(.)$ext1/$1$ext1/;
                    my $check_ext2 = $read_end9 =~ s/(.)$ext2/$1$ext2/;
                    if ($check_ext1 > 0 || $check_ext2 > 0)
                    {
                        $split_no_rep = "yes";
                    }
                }
                if ($deletion ne "yes" && $reference_guided ne "yes" && $type ne "mito_plant" && $last_chance ne "yes" && $SNP_active eq "" && $contig_end ne "yes" && $indel_split_skip ne "yes" && ($delete_first ne "yes" || (length($best_extension1) < 5 && length($best_extension2) < 5)) && $ext_total > 18 && $indel_split < ($read_length-25 -$overlap))
                {
                    $best_extension = "";
                    $indel_split{$id} = $indel_split+10;
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "INCREASE_INDEL_SPLIT\n";
                         print OUTPUT5 $indel_split." INDEL_SPLIT\n";
                    }
                    $read_new = $read;     
                    $read_new1 = $read_new;
                    goto BACK;
                }
                elsif (($SNP_active eq "" || ($extensions_before eq "yes" && $deletion ne "yes" && $reference_guided ne "yes")) && $contig_end ne "yes" && ($delete_first ne "yes" || (length($best_extension1) < 5 && length($best_extension2) < 5)) && $deletion ne "yes" && $reference_guided ne "yes")
                {
                    $best_extension = "";
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "SNP_ACTIVE\n";
                    }
                    $SNP_active{$id} = "yes";
                    delete $indel_split{$id};

                    $read_new = $read;     
                    $read_new1 = $read_new;
                    goto BACK;
                }
                elsif ($deletion eq "yes" || $reference_guided eq "yes")
                {
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "DELETION_DETECT\n";
                    }
                    delete $seed{$id_split1};
                    delete $seed{$id};
                    delete $before{$id};
                    delete $before_shorter_skip{$id};
                    $read_new = $read;     
                    $read_new1 = $read_new;
                }
                elsif ($jump_rep ne "yes" && $repetitive_detect eq "yes" && ($split_no_rep eq "yes" || $repetitive_detect2 eq "yes"))
                {
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "JUMP REP\n";
                    }
                    delete $indel_split{$id};
                    $jump_rep{$id} = "yes";
                    $best_extension = "";
                    $read_new = $read;     
                    $read_new1 = $read_new;
                    goto BACK;
                }
                
                
                elsif($delete_first ne "yes" && $deletion ne "yes" && $reference_guided ne "yes")
                {            
                    if (exists $old_id{$id} && ($noback eq "stop" || $position_back >= ($insert_size*3)))
                    {                                     
                        $merge_read_length = length ($read);
                        $merge_read = "yes"; 
                        $read = $seed_old{$old_id{$id}} ."LLLLLLLLLLLLLLL".$read;
                        $seed{$id} = $read;
                        $seeds_check{$id} = undef;
                        $hasL = "yes";
                        foreach my $tree_tmp (keys %tree)
                        {
                            my $old = $old_id{$id};
                            my $tree2 = $tree{$tree_tmp};
                            my $tree3 = $tree{$tree_tmp};
                            if ($old_id{$id} =~ m/.*_(\d+)$/)
                            {
                                $old = $1;
                            }
                            my $id_tmp = $id;
                            if ($id =~ m/.*_(\d+)$/)
                            {
                                $id_tmp = $1;
                            }
                            my @ids_split = split /\*/, $tree2;
                            foreach my $id_split (@ids_split)
                            {
                                if ($id_split  =~ m/^$old(REP)*$/)
                                {
                                    if ($tree2 =~ m/^(.*\*)*$old(REP)*(\*.*)*$/)
                                    {
                                        if (defined($1))
                                        {
                                            $tree3 = $1.$id_tmp;
                                        }
                                        else
                                        {
                                            $tree3 = $id_tmp;
                                        }
                                        if (defined($2))
                                        {
                                            $tree3 = $tree3."REP";
                                        }
                                        if (defined($3))
                                        {
                                            $tree3 = $tree3.$3;
                                        }
                                    }
                                }
                            }
                            delete $tree{$tree_tmp};
                            $tree{$tree_tmp} = $tree3;
                            foreach my $contigs_end (keys %contigs_end)
                            {
                                if ($contigs_end{$contigs_end} eq $old)
                                {
                                    delete $contigs_end{$contigs_end};
                                    $contigs_end{$contigs_end} = $id_tmp;
                                }
                            }
                        }
                        delete $old_id{$id};
                        $noback{$id} = "stop";
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "Merged contigs with LLLLLLLLLLL!\n";
                        }
                        $contig_gap_min{$id."_".$contig_count} = ($contig_gap_min{$id."_".$contig_count}-$position_back);
                        $contig_gap_max{$id."_".$contig_count} = ($contig_gap_max{$id."_".$contig_count}-$position_back);
                    }
                    
                    $indel_split = '0';
                    delete $indel_split{$id};
                    $best_extension = "";
                    
                    delete $seed{$id_split1};
                    delete $seed{$id};
                    $seed{$id} = $read;
                    $seeds_check{$id} = undef;
                    $read_new = $read;
                    $read_new1 = $read;                 
                    
                    if ($SNR_read eq "yes" && length($best_extension1) < 7 && length($best_extension2) < 7)
                    {
                        $SNR_next_seed = "yes";
                        goto BACK;
                    }
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 ">".$id."\n";
                        print OUTPUT5 $read."\n\n\n";
                    }
                    my $best_ex1 = substr $best_extension1, 0, 7;
                    my $best_ex2 = substr $best_extension2, 0, 7;
                    $best_ex1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    $best_ex2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    my $contigs_endi = substr $read_end, -10;
                    my $read_short_end2_tmp = $read_short_end2;
                    my $check_rep1 = $read_short_end2_tmp =~ s/$contigs_end0$best_ex1/$contigs_end0$best_ex1/g;
                    my $check_rep2 = $read_short_end2_tmp =~ s/$contigs_end0$best_ex2/$contigs_end0$best_ex2/g;
                    
CORRECT:            my $contig_id2_tmp = substr $contig_id2, 0,-1;
                                         
                    my $contigs_end1 = substr $best_extension1, 0, 7;
                    my $contigs_end2 = substr $best_extension2, 0, 7;
                    $contigs_end1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    $contigs_end2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    my $contigs_end0 = substr $read_end, -5;
                    my $tree_empty = "";
                    
                    if (exists($contigs_id{$contig_id2}) || exists($contigs_end{$contigs_end0.$contigs_end2}))
                    {
                        $tree{$id} = $contigs_end{$contigs_end0.$contigs_end2};
                        $tree_empty = "yes";
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "TREE_EMPTY\n";
                        }
                    }
                    elsif ($no_contig_id2 eq "yes")
                    {
                        $nosecond{$id} = undef;
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "NO_CONTIG_ID2\n";
                        }
                    }
                    else
                    {
                        $seed{$contig_id2} = $contig_read2;
                        $seeds_check{$contig_id2} = undef;
                        $insert_size2{$contig_id2} = $insert_size;
                        $position{$contig_id2} = length($contig_read2);
                        $old_id2{$contig_id2} = undef;
                        $noback{$contig_id2} = "stop";
                        
                        $contigs_id{$contig_id2} = undef;
                        $contigs_end{$contigs_end0.$contigs_end2} = $contig_id2;
                        $correct_after_split = "yes";
                        $tree{$id} = $contig_id2;
                        
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 $contig_id2." CONTIG_ID2\n";
                            print OUTPUT5 $best_extension2." best_ext2\n";
                            my $contig_id2hhh = substr $contig_id2, 0 ,-1;
                            if (exists($hash{$contig_id2hhh}))
                            {
                                print OUTPUT5 $hash{$contig_id2hhh}." HASH2\n";
                            }
                        }
                        if (length($read) < 251)
                        {
                            delete $seed{$contig_id2};
                            $seed{$contig_id2} = $read.$best_extension2;
                            $seeds_check{$contig_id2} = undef;
                            delete $tree{$id};
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "SHORT2\n";
                            }
                            foreach my $tree_tmp (keys %tree)
                            {
                                my $tmp = $tree{$tree_tmp};
                                my @ids_split = split /\*/, $tmp;
                                foreach my $id_split (@ids_split)
                                {
                                    if ($id_split eq $id)
                                    {
                                        delete $tree{$tree_tmp};
                                        if ($tmp =~ m/^(.*\*)*$id(\*.*)*$/)
                                        {
                                            if (defined($1))
                                            {
                                                $tmp = $1.$contig_id2
                                            }
                                            else
                                            {
                                                $tmp = $contig_id2
                                            }
                                            if (defined($2))
                                            {
                                                $tmp = $tmp.$2
                                            }
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $tmp." TESTOOO5\n";
                                            }
                                            $tree{$tree_tmp} = $tmp;
                                            foreach my $end_tmp (keys %contigs_end)
                                            {
                                                if ($contigs_end{$end_tmp} eq $id)
                                                {
                                                    delete $contigs_end{$end_tmp};
                                                    $contigs_end{$end_tmp} = $contig_id2;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (exists($contigs_id{$contig_id1}) || exists($contigs_end{$contigs_end0.$contigs_end1}))
                    {
                        if ($tree_empty eq "yes")
                        {
                            $tree{$id} = $contigs_end{$contigs_end0.$contigs_end2}."*".$contigs_end{$contigs_end0.$contigs_end1};
                        }
                        else
                        {
                            $tree{$id} = $contig_id2."*".$contigs_end{$contigs_end0.$contigs_end1};
                        }
                    }
                    elsif ($no_contig_id1 eq "yes")
                    {
                        $nosecond{$id} = undef;
                    }
                    else
                    {
                        
                        $seed{$contig_id1} = $contig_read1;
                        $seeds_check{$contig_id1} = undef;
                        $insert_size2{$contig_id1} = $insert_size;
                        $position{$contig_id1} = length($contig_read1);
                        $old_id2{$contig_id1} = undef;
                        $noback{$contig_id1} = "stop";
                        
                        $contigs_id{$contig_id1} = undef;
                        $contigs_end{$contigs_end0.$contigs_end1} = $contig_id1;
                        $correct_after_split = "yes";
                        
                        if ($tree_empty eq "yes")
                        {
                            $tree{$id} =  $contigs_end{$contigs_end0.$contigs_end2}."*".$contig_id1;
                        }
                        else
                        {
                            $tree{$id} = $contig_id2."*".$contig_id1;
                        }
                                                                    
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 $contig_id1." CONTIG_ID1\n";
                            print OUTPUT5 $best_extension1." best_ext1\n";
                            my $contig_id1hhh = substr $contig_id1, 0 ,-1;
                            if (exists($hash{$contig_id1hhh}))
                            {
                                print OUTPUT5 $hash{$contig_id1hhh}." HASH1\n";
                            }
                        }
                        if (length($read) < 251)
                        {
                            delete $seed{$contig_id1};
                            $seed{$contig_id1} = $read.$best_extension1;
                            $seeds_check{$contig_id1} = undef;
                            delete $tree{$id};
                            foreach my $tree_tmp (keys %tree)
                            {
                                my $tmp = $tree{$tree_tmp};
                                my @ids_split = split /\*/, $tmp;
                                foreach my $id_split (@ids_split)
                                {
                                    if ($tree_empty eq "yes")
                                    {
                                        if ($id_split eq $id)
                                        {
                                            delete $tree{$tree_tmp};
                                            if ($tmp =~ m/^(.*\*)*$id(\*.*)*$/)
                                            {
                                                if (defined($1))
                                                {
                                                    $tmp = $1.$contig_id1
                                                }
                                                else
                                                {
                                                    $tmp = $contig_id1
                                                }
                                                if (defined($2))
                                                {
                                                    $tmp = $tmp.$2
                                                }
                                                $tree{$tree_tmp} = $contigs_end{$contigs_end0.$contigs_end2}."*".$tmp;
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 $tree{$tree_tmp}." TESTOOO4a-----------------\n";
                                                }
                                                foreach my $end_tmp (keys %contigs_end)
                                                {
                                                    if ($contigs_end{$end_tmp} eq $id)
                                                    {
                                                        delete $contigs_end{$end_tmp};
                                                        $contigs_end{$end_tmp} = $contig_id1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if ($id_split eq $contig_id2)
                                        {
                                            delete $tree{$tree_tmp};
                                            $tree{$tree_tmp} = $contig_id1."*".$tmp;
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $tree{$tree_tmp}." TESTOOO5------------\n";
                                            }
                                            foreach my $end_tmp (keys %contigs_end)
                                            {
                                                if ($contigs_end{$end_tmp} eq $contig_id2)
                                                {
                                                    delete $contigs_end{$end_tmp};
                                                    $contigs_end{$end_tmp} = $contig_id2."*".$contig_id1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }                    
                    if ($noback eq "stop")
                    {                    
                        if (exists($old_id{$id}))
                        {
                            my $read_tmp = $seed_old{$old_id{$id}};
                            if (length($read_tmp) > 250)
                            {
                                if ($contig_num eq '1')
                                {
                                    $contigs{$contig_num."+".$id} = $read_tmp;
                                    my $start_point = '500';
                                    my $check_repetitive = '3';
                                                                                                
                                    while ($check_repetitive > 2)
                                    {
                                        my $repetitive = substr $read_tmp, $start_point, 15;
                                        $repetitive =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                        my $read_short_area = substr $read_tmp, $start_point -170, 340;
                                        $check_repetitive = $read_short_area =~ s/$repetitive/$repetitive/g;
                                        if ($check_repetitive > 2)
                                        {
                                            $start_point += 20;
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 "DETECT_REPETITIVE_IN_START_SEQUENCE\n";
                                            }
                                        } 
                                    }
                                    $first_contig_start = substr $read_tmp, $start_point, $overlap;
                                    my $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                    while ($check_start > 0)
                                    {
                                        $start_point += 10;
                                        $first_contig_start = substr $read_tmp, $start_point, $overlap;
                                        $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                    }
                                    if ($y > $startprint2)
                                    {
                                        print OUTPUT5 $first_contig_start." CONTIG_START3\n";
                                    }
                                    foreach my $seedie (keys %seed)
                                    {
                                        my $seedie_part = substr $seed{$seedie}, 0, 1000;
                                        my $check_seedie = $seedie_part =~ s/$first_contig_start/$first_contig_start/;
                                        if ($check_seedie > 0)
                                        {
                                            delete $seed{$seedie};
                                        }
                                    }
                                }
                                else
                                {
                                    $contigs{$contig_num."+".$old_id{$id}} = $read_tmp;                                
                                }
                                $contig_num++;    
                            }
                        }
                        if (length($read) > 250)
                        {   
                            if ($contig_num eq '1')
                            {
                                $contigs{$contig_num."+".$id} = $read;
                                my $start_point = '500';
                                my $check_repetitive = '3';
                                                                                            
                                while ($check_repetitive > 2)
                                {
                                    my $repetitive = substr $read, $start_point, 15;
                                    $repetitive =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                    my $read_short_area = substr $read, $start_point -170, 340;
                                    $check_repetitive = $read_short_area =~ s/$repetitive/$repetitive/g;
                                    if ($check_repetitive > 2)
                                    {
                                        $start_point += 20;
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 "DETECT_REPETITIVE_IN_START_SEQUENCE\n";
                                        }
                                    } 
                                }
                                $first_contig_start = substr $read, $start_point, $overlap;
                                my $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                while ($check_start > 0)
                                {
                                    $start_point += 10;
                                    $first_contig_start = substr $read, $start_point, $overlap;
                                    $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                }
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 $first_contig_start." CONTIG_START4\n";
                                }
                                foreach my $seedie (keys %seed)
                                {
                                    my $seedie_part = substr $seed{$seedie}, 0, 1000;
                                    my $check_seedie = $seedie_part =~ s/$first_contig_start/$first_contig_start/;
                                    if ($check_seedie > 0)
                                    {
                                        delete $seed{$seedie};
                                    }
                                }
                            }
                            else
                            {
                                $contigs{$contig_num."+".$id} = $read;
                            }
                            $contig_num++;
                        }
                        delete $seed{$id};
                        if (!keys %seed)
                        {
                            $circle = "contigs";
                            goto FINISH;
                        }
                        else
                        {
                            goto ITERATION;
                        }
                    }
                    else
                    {
                        $noforward = "stop";
                        $noforward{$id} = "stop";
                        $seed_split{$id} = undef;
                        $best_extension = "";
                        $nosecond{$id} = undef;
                        goto BACK;
                    }
                }
                else
                {
                    $indel_split = '0';
                    delete $indel_split{$id};
                }
            }
            else
            {                             
                if ($best_extension ne "")
                {
                    $indel_split = '0';
                    delete $indel_split{$id};
                }

                             
                my $best_extension_no_dot = $best_extension;
                my @ext = split //, $best_extension;
                                   
                my $u = length($best_extension);
                my $v = '1';
                if (($SNR_test eq "yes2" || $SNR_test eq "yes2_double") && $u > 2)
                {
                    my @ext3 = split //, $best_extension;
                    $v = '1';
                    $u = length($best_extension);
                    while ($ext3[$u-$v-1] eq $ext3[$u-$v] && ($u-$v-1) > 1)
                    {              
                        chop($best_extension);
                        $v++;
                    }
                }
                    if (($SNR_test eq "yes222222" || $SNR_test eq "yes2_double") && $u > 2)
                    {
                        my @ext3 = split //, $best_extension;
                        $v = '1';
                        $u = length($best_extension);
                        while ($ext3[$u-$v-1] eq $ext3[$u-$v] && ($u-$v-1) > 1)
                        {              
                            chop($best_extension);
                            $v++;
                        }
                        my $SNR_max = '0';
                        my $SNR_min = '1000';
                        my $n = '0';
                        my $extensions_after_SNR = substr $best_extension, 0,4;

                        foreach my $SNR_ext (keys %SNR_count)
                        {
                            my $SNR_count = $SNR_count{$SNR_ext};
                            my $checkie = substr $SNR_ext, 0, 4;
                            if ($SNR_count ne "" && $checkie eq $extensions_after_SNR)
                            {
                                if ($SNR_count > $SNR_max)
                                {
                                    $SNR_max = $SNR_count;
                                }
                                if ($SNR_count < $SNR_min)
                                {
                                    $SNR_min = $SNR_count;
                                }
                            }    
                        }
                        if ($SNR_min eq '1000')
                        {
                            $SNR_min = '0';
                        }
                        my $p = '0';
                        my $ut = '0';

                        if ($SNR_test eq "yes2" && $SNR_max > 0)
                        {
                            $best_extension = $SNR_nucleo.$best_extension;
                            $ut = '1';
                            while ($p < $SNR_max - $SNR_min)
                            {
                                $best_extension = "X".$best_extension;
                                $p++;
                            }
                        }
                        elsif ($SNR_test eq "yes2_double"  && $SNR_max > 0)
                        {
                            $best_extension = $SNR_nucleo.$best_extension;
                            $ut = '2';
                            while ($p < ($SNR_max - $SNR_min)/2)
                            {
                                $best_extension = "X2".$best_extension;
                                $p++;
                            }
                        }
                        $p = '0';
                        
                        while ($p+$ut < $SNR_min)
                        {
                            $best_extension = $SNR_nucleo.$best_extension;
                            if ($SNR_test eq "yes2")
                            {
                                $p++; 
                            }
                            if ($SNR_test eq "yes2_double")
                            {
                                $p++;
                                $p++; 
                            }     
                        }
                        delete $SNR{$id};
                    }
                    elsif ($SNR_test eq "yes2222222" || $SNR_test eq "yes2_double")
                    {
                        my $ee2 = '10000'; 
                        foreach (keys %SNR_count)
                        {
                            my $ee = $SNR_count{$_}; 
                            if ($ee < $ee2)
                            {
                                $ee2 = $ee;
                            }
                        }
                        $best_extension = "";
                        my $ee3 = '0';
                        if ($ee2 ne '10000')
                        {
                             while ($ee3 < $ee2)
                            {
                                $best_extension = $SNR_nucleo.$best_extension;
                                $ee3++;
                            }  
                        }
                        else
                        {
                            $best_extension = "";
                        }
                    }
                $SNR_test = "";
                if ($y > $startprint2)
                {
                    print OUTPUT5 $best_extension." BEST_EXTENSION\n\n";
                }
                $best_extension_prev{$id} = $best_extension;    
            }  
AFTER_EXT:

                                            chomp $best_extension;
                                            my $vk2 = '0';
                                            if ($SNR_read eq "")
                                            {
                                                my $best_extension_tmp7 = $best_extension;
                                                $best_extension_tmp7 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                
                                                my @dot2 = split //, $best_extension_tmp7;
                                                my $ut2 = length($best_extension_tmp7);
                                                
                                                while ($dot2[$ut2-1] eq "." || $dot2[$ut2-1] eq "*")
                                                {              
                                                    if ($dot2[$ut2-1] eq "*")
                                                    {
                                                        chop $best_extension;
                                                        chop $best_extension;
                                                        $vk2++;
                                                        $vk2++;
                                                        $ut2--;
                                                        $ut2--;
                                                    }
                                                    else
                                                    {                                                     
                                                        chop $best_extension;
                                                        $vk2++;
                                                        $ut2--;
                                                    }
                                                }
                                            }
                                                my $best_extension_tmpp = substr $best_extension, -15;
                                                my $SNR_check = $best_extension_tmpp =~ qr/AAAAAAA|CCCCCCC|GGGGGGG|TTTTTTT/;
                                                if ($SNR_check > 0)
                                                {
                                                    if ($best_extension =~ m/(.*?(AAAAAA|CCCCCC|GGGGGG|TTTTTT))(.*)/)
                                                    {
                                                        my $best_extension_short = $1;
                                                        if (length($3) < length($best_extension_short)+5)
                                                        {
                                                            $best_extension = $best_extension_short;
                                                            if ($y > $startprint2)
                                                            {
                                                                print OUTPUT5 $best_extension." BEST_EXTENSION_SHORT\n\n";
                                                            }
                                                        }
                                                    }
                                                }
                                            
 
                                            if ($best_extension ne "" && $best_extension ne " ")
                                            {                                                                               
                                                $read_new = $read.$best_extension;                                           
                                                $read_new1 = $read_new;
                                                
                                                $position += length($best_extension);
                                                $position -= $vk2;                                    
                                                                                 
EXTEND_READ:
                                                delete $position{$id};
                                                $position{$id} = $position;

                                                if ($best_extension ne "" && $split_forward eq "")
                                                {
                                                    delete $indel_split_skip{$id};
                                                    delete $before{$id};
                                                    delete $before_shorter_skip{$id};
                                                }
                                                
                                                $best_extension = "";
                                                delete $regex{$id};
                                                if ($split eq "")
                                                {
                                                    delete $yuyu_option{$id.'A'};
                                                    delete $yuyu_option{$id.'C'};
                                                    delete $yuyu_option{$id.'T'};
                                                    delete $yuyu_option{$id.'G'};
                                                }
                                                $use_regex = "";
                                                if ($SNR_read ne "" && $last_chance eq "yes")
                                                {
                                                    delete $last_chance{$id};
                                                    $last_chance{$id} = "yes";
                                                }
                                                else
                                                {
                                                    delete $last_chance{$id};
                                                }
                                                                        
                                                $id_test = $id;
                                                if ($split eq "")
                                                {
                                                    foreach my $add_read2 (keys %extensions)
                                                    {
                                                        my $add_read = substr $extensions{$add_read2}, 0, -1;
                                                        $count_reads{$add_read} = undef;
                                                    }
                                                }
                                            }
                                            elsif ($use_regex ne "yes" && $indel_split > 0)
                                            {                                              
                                                delete $regex{$id};
                                                if ($split_forward eq "")
                                                {
                                                    $regex{$id} = "yes";
                                                }
                                                elsif ($split_forward ne "" && exists($SNP_active{$id}))
                                                {
                                                    $indel_split_skip{$id} = "yes";
                                                    $indel_split = '0';
                                                    delete $indel_split{$id};
                                                }
                                                elsif ($split_forward ne "")
                                                {
                                                    $SNP_active{$id} = "yes";
                                                }
                                                $read_new = $read;
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "OPTION1\n";
                                                }
                                            }
                                            elsif ($use_regex ne "yes" && $repetitive_detect eq "")
                                            {                                              
                                                delete $regex{$id};
                                                $regex{$id} = "yes";
                                                $read_new = $read;
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "OPTION2\n";
                                                }
                                            }
                                            elsif ($use_regex ne "yes" && $last_chance ne "yes" && $indel_split > 0)
                                            {
                                                $read_new = $read;
                                                $indel_split_skip{$id} = "yes";
                                                $indel_split = '0';
                                                delete $indel_split{$id};
                                                $SNP_active{$id} = "yes";
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "OPTION3\n";
                                                }
                                            }
                                            elsif (($use_regex eq "yes" || $repetitive_detect ne "") && $last_chance ne "yes" && $AT_rich_before eq "")
                                            {
                                                $read_new = $read;
                                                delete $last_chance{$id};
                                                $last_chance{$id} = "yes";
                                                delete $regex{$id};
                                                $use_regex = "";
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "OPTION4\n";
                                                }
                                            }
                                            elsif ($last_chance eq "yes" && $use_regex eq "")
                                            {
                                                $read_new = $read;
                                                delete $last_chance{$id};
                                                $last_chance{$id} = "yes";
                                                delete $regex{$id};
                                                $regex{$id} = "yes";
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "OPTION5\n";
                                                }
                                            }
                                            elsif (($last_chance eq "yes" || $AT_rich_before ne "") && $use_regex ne "")
                                            {
                                                $read_new = $read;
                                                delete $last_chance{$id};
                                                $noforward = "stop";
                                                $noforward{$id} = $noforward;
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "OPTION6\n";
                                                }
                                            }
                                            else
                                            {
                                                delete $seed{$id};
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "OPTION7\n";
                                                }
                                            }
                                            if ($split eq "yes2")
                                            {                      
                                                $seed{$id} = $read_new;
                                                $insert_size2{$id} = $insert_size2;
                                                
                                                if ($count_split eq '2')
                                                {
                                                    $split = "yes4";
                                                }
                                                goto SPLIT;
                                            }
                                            elsif ($split eq "yes3")
                                            {                      
                                                $seed{$id} = $read_new;
                                                $insert_size2{$id} = $insert_size2;

                                                if ($count_split eq '3')
                                                {
                                                    $split = "yes4";
                                                }                                               
                                                goto SPLIT;
                                            }
                                            if ($split eq "yes4")
                                            {                      
                                                $seed{$id} = $read_new;
                                                $insert_size2{$id} = $insert_size2;
                                                
                                                goto SPLIT;
                                            }
                                            
BACK:   if (keys %merged_match_back eq 0 && $use_regex_back ne "yes" && $noback eq "" && $circle eq "" )
        {                                              
            delete $regex_back{$id};
            $regex_back{$id} = "yes";
            $best_extension_back_prev{$id} = "";
        }
        elsif (keys %merged_match_back eq 0 && $use_regex_back eq "yes")
        {
            $noback{$id} = "stop";
            $best_extension_back_prev{$id} = "";
        }
        if (keys %merged_match_back > 0 && $noback ne "stop" )
        {
            undef %extensions1;
            undef %extensions2;
            undef %extensions;
            undef %extensions1b;
            undef %extensions2b;
            undef %extensionsb;
            undef @extensions1;
            undef @extensions2;
            undef @extensions;
            undef %extensions1b;
            undef %extensions2b;
            undef %extensionsb;
            undef @matches;
            undef @matches1;
            undef @matches2;
            undef %SNR_length;
            undef %filter_before1;
            undef %filter_before2;
            undef %filter_before3;
            undef %filter_before4;
            undef %filter_before1_pair;
            undef %filter_before2_pair;
            undef %filter_before3_pair;
            undef %filter_before4_pair;
            $extra_regex = "";
            $split_forward = "";
            $split = "";
            $extensions_before = "";
            $ext_before = "";
            $insert_range_shorter = "";
        
REGEX_BACK:
            $read_count = '0';
            $read_ex = '0';
            $read2_ex = '0';

            if ($SNR_read_back ne "")
            {
            }
            my $X4 = $read_short_start =~ tr/\*//;
            
            if ($X4 > 0)
            {
                %read_short_start_tmp = build_partial3c ($read_short_start, "reverse");
            }
            else
            {
               $read_short_start_tmp{$read_short_start} = undef;
            } 
  
            if ($y > $startprint2)
            {
                if ($extra_regex eq "yes" || $use_regex_back eq "yes")
                {
                    print OUTPUT5 "USE_REGEX_BACK_REVERSE\n";
                }
            }

NO_MATCH_BACK:  foreach my $ln (keys %merged_match_back)
                {              
                    $match = $merged_match_back{$ln};
                    $id_match = $ln;
                    chomp $id_match;
                    chomp $match;
    
                    my %match;
                    my %match2;
                    undef %match;
                    undef %match2;
            
                    my $n = '1';
                    while ($n < length($match)-($overlap+1))
                    {
                        my $match_part = substr $match, $n, $overlap;
                        my $match_part2 = substr $match, 0, $n;
                        my $match_part3 = substr $match, $n+$overlap;
                        $match{$match_part} = $match_part2;
                        $match2{$match_part} = $match_part3;
                        $n++;     
                    }
                    if ($last_chance_back eq "yes")
                    {                             
                        my $forward = "";
                        foreach my $line (keys %read_start_b)
                        {
                            my @read_start_b_sub = build_partialb $line;
                                                        
                            my $found_seq = '0';
                            my $match2 = $match;
                                                                
                            foreach my $read_start_b_sub (@read_start_b_sub)
                            {
                                $found_seq = $match2 =~ s/.$read_start_b_sub/+/;
                                if ($found_seq > 0)
                                {
                                    my @ext = split /\+/, $match2;
                                    my $extension5 = $ext[1];
                                    $extension5 =~ tr/ATCG/TAGC/;
                                    $extension = $extension5;                                                  
                                    $read_count++;
                                    goto LAST1_BACK;
                                }
                            }
                        }
                        foreach my $line (keys %read_start)
                        {
                            my @read_start_sub = build_partialb $line;                                
                            my $found_seq = '0';
                            my $match4 = $match;
                                                                
                            foreach my $read_start_sub (@read_start_sub)
                            {
                                $found_seq = $match4 =~ s/(.)$read_start_sub/$1+/g;
                                if ($found_seq > 1)
                                {
                                    my $pos = $merged_match_back_pos{$ln};
                                    my $match4b = substr $match, 0, -$pos;
                                    $match4b =~ s/(.+)$line/$1+/;
                                    my @ext = split /\+/, $match4;
                                    $extension = reverse($ext[0]);                                                
                                    $read_count++;
                                    $forward = "yes";
                                    goto LAST1_BACK;
                                }
                                elsif ($found_seq > 0)
                                {
                                    my @ext = split /\+/, $match4;
                                    $extension = reverse($ext[0]);                                                
                                    $read_count++;
                                    $forward = "yes";
                                    goto LAST1_BACK;
                                }
                            }
                        }
                        next NO_MATCH_BACK;
LAST1_BACK:             my $id_match_end = substr $id_match, -1, 1,"",;
                        my $id_match_id = substr $id_match, 0, -2,;
                                   
                        if (index ($id_match_id, $id) eq "-1")
                        {
                            if ($id_match_end eq "1")
                            {
                                $id_pair_match = $id_match."2";
                            }
                            elsif ($id_match_end eq "2")
                            {
                                $id_pair_match = $id_match."1";
                            }
    
                            my $nuc_exlude = "test";
                            if ($yuyu_option_A_back eq "A" || $yuyu_option_C_back eq "C" || $yuyu_option_T_back eq "T" || $yuyu_option_G_back eq "G")
                            {
                                $nuc_exlude = substr $extension, 0, 1;
                            }
                            if ($extension ne " " && $extension ne "" && $nuc_exlude ne $yuyu_option_A_back && $nuc_exlude ne $yuyu_option_C_back && $nuc_exlude ne $yuyu_option_T_back && $nuc_exlude ne $yuyu_option_G_back)
                            {
                                push @matches2, $id_match.",".$extension.","."".",".$match.","."";
                                $extensions2{$extension} = $id_match;
                                push @extensions2, $extension;
                                if ($forward eq "yes")
                                {
                                    $extensions2b{$id_match} = $extension;
                                }
                                else
                                {
                                    $extensions1b{$id_match} = $extension;
                                }
                            }
                        }
                        next NO_MATCH_BACK;
                    }                   
                    elsif ($extra_regex eq "yes" || $use_regex_back eq "yes")
                    {  
                        if ($repetitive_detect_back eq "")
                        {  
                            foreach my $line (keys %read_start)
                            {
                                if (exists($match{$line}))
                                {
                                    $extension = reverse($match{$line});                                                                    
                                    $read_count++;;
                                    goto FOUND_BACK;
                                }
                            }
                        }
                        foreach my $line (keys %read_start)
                        {
                            my @read_start_sub = build_partialb $line;                                
                            my $found_seq = '0';
                            my $match4 = $match;
                                                                
                            foreach my $read_start_sub (@read_start_sub)
                            {
                                $found_seq = $match4 =~ s/(.)$read_start_sub/$1+/g;
                                if ($found_seq > 1)
                                {
                                    my $pos = $merged_match_back_pos{$ln};
                                    my $match4b = substr $match, 0, -$pos;
                                    $match4b =~ s/(.+)$line/$1+/;
                                    my @ext = split /\+/, $match4;
                                    $extension = reverse($ext[0]);
                                    $read_count++;
                                    goto FOUND_BACK; 
                                }
                                elsif ($found_seq > 0)
                                {
                                    my @ext = split /\+/, $match4;
                                    $extension = reverse($ext[0]);
                                    $read_count++;
                                    goto FOUND_BACK;     
                                }
                            }
                        }
                        next NO_MATCH_BACK;
                    }
                    else
                    {
                        if ($repetitive_detect_back eq "")
                        {
                            foreach my $line (keys %read_start)
                            {
                                if (exists($match{$line}))
                                {
                                    $extension = reverse($match{$line});                                                                    
                                    $read_count++;;
                                    goto FOUND_BACK;
                                }
                            }
                        }
                        foreach my $line (keys %read_start)
                        {                             
                            my $found_seq = '0';
                            my $match4 = $match;
                                                                
                            $found_seq = $match4 =~ s/(.)$line/$1+/g;
                            if ($found_seq > 1)
                            {
                                my $pos = $merged_match_back_pos{$ln};
                                my $match4b = substr $match, 0, -$pos;
                                $match4b =~ s/(.+)$line/$1+/;
                                my @ext = split /\+/, $match4;
                                $extension = reverse($ext[0]);
                                $read_count++;
                                goto FOUND_BACK; ; 
                            }
                            elsif ($found_seq > 0)
                            {
                                my @ext = split /\+/, $match4;
                                $extension = reverse($ext[0]);
                                $read_count++;
                                goto FOUND_BACK; ;
                            }                          
                        }
                        next NO_MATCH_BACK;
                    }
FOUND_BACK:     if ($last_chance_back eq "yes")
                {           
                    next NO_MATCH_BACK;
                }         
                    if ($extension ne "NOOO")
                    {
                            my $id_match_b = $id_match;
                            my $id_match_end = substr $id_match_b, -1, 1,"",;
                        
                            if (exists($hash{$id_match_b}))
                            {
                                my @id_match_b = split /,/, $hash{$id_match_b};
                                
                                if ($id_match_end eq "1")
                                {
                                    $match_pair = $id_match_b[1];
                                }
                                elsif ($id_match_end eq "2")
                                {
                                    $match_pair = $id_match_b[0];
                                }
                                else
                                {
                                    next NO_MATCH_BACK;
                                }                                
                                chomp($match_pair);
                                if ($encrypt eq "yes")
                                {
                                    $match_pair = decrypt $match_pair;
                                }
                                
                                my $is = $overlap;
                                if ($indel_split_back ne '0')
                                {
                                    $is = length($match_pair2)-40;
                                    $is = $overlap+$indel_split_back;
                                }
                                
                                $match_end = substr $match_pair, -($overlap+$right), $overlap;
                                $match_start = substr $match_pair, $left, $overlap;
                                                                  

                                                my $match_pair3 = $match_pair;
                                                $match_pair3 =~ tr/ATCG/TAGC/;
    
                                                $match_pair2 = reverse ($match_pair3);
                                                
                                                my $match_end3 = $match_end;
                                                $match_end3 =~ tr/ATCG/TAGC/;
                                                    
                                                $match_end2 = reverse ($match_end3);
                                                
                                                my $match_start3 = $match_start;
                                                $match_start3 =~ tr/ATCG/TAGC/;
                                                    
                                                $match_start2 = reverse ($match_start3);
                                                
                                                my %match_pair2;
                                                undef %match_pair2;
        
                                                my $v = '1';
                                                while ($v < (length($match_pair2)-($overlap+1)))
                                                {
                                                    my $match_pair2_part = substr $match_pair2, $v, $overlap;
                                                    my $match_pair2_part2 = substr $match_pair2, 0, $v;
                                                    $match_pair2{$match_pair2_part} = $match_pair2_part2;
                                                    $v++;
                                                }
                                                if ($extra_regex eq "yes")
                                                {                                                   
                                                    my $match_pair_middle2 = substr $match_pair2, -26-$is, $is+10;
                                                    $match_pair_middle2 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    my @match_pair_middle_sub = split //, $match_pair_middle2;
                                                    
                                                    my $F = ($insert_size*$insert_range_back)-length($extension) - (((($insert_size*$insert_range_back)-$insert_size)*2)+$is+10+8);
                                                    if ($F < 0)
                                                    {
                                                        $F = '0';
                                                    }

                                                    my $read_short_start_tempie = substr $read, $F, ((($insert_size*$insert_range_back)-$insert_size)*2)+$is+10+8;                                                
                                                    $read_short_start_tempie =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    undef %read_short_start_tmp;
                                                    %read_short_start_tmp = build_partial3c $read_short_start_tempie;
                                                                                         
                                                    foreach my $line (keys %read_short_start_tmp)
                                                    {
                                                        my @line = split //,$line;
                                                        my $gh = '0';
                                                        my $th = '0';
                                                        $line =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                     
CHECK_PAIR_BACK:                                        while ($gh < length($line)-length($match_pair_middle2))
                                                        {
                                                            my $d = '0';
                                                            my $next = "";
                                                                                                
                                                            while ($d < length($match_pair_middle2))
                                                            {
                                                                $th = $d + $gh;
                                                                if ($match_pair_middle_sub[$d] eq $line[$th])
                                                                {
                                                                }
                                                                elsif ($line[$th] eq ".")
                                                                {
                                                                }
                                                                elsif ($match_pair_middle_sub[$d] eq ".")
                                                                {
                                                                }
                                                                elsif ($next eq "")
                                                                {
                                                                    $next = "yes";
                                                                }
                                                                elsif ($next eq "yes")
                                                                {
                                                                    $next = "yes2";
                                                                }
                                                                elsif ($next eq "yes2")
                                                                {
                                                                    $next = "yes3";
                                                                }
                                                                else
                                                                {
                                                                    $gh++;
                                                                    goto CHECK_PAIR_BACK;
                                                                }
                                                                $d++    
                                                            }                                                         
                                                            $extension_match = "";
                                                            goto SKIP_BACK;    
                                                        }
                                                    }   
                                                    $extension_match = "NOOO";
                                                }
                                                else
                                                {
                                                    my $match_pair_middle2 = substr $match_pair2, -26-$is, $is+10;
                        
                                                    my $F = ($insert_size*$insert_range_back)-length($extension) - (((($insert_size*$insert_range_back)-$insert_size)*2)+$is+10+8);
                                                    if ($F < 0)
                                                    {
                                                        $F = '0';
                                                    }
                                                    my $read_short_start_tempie = substr $read, $F, ((($insert_size*$insert_range_back)-$insert_size)*2)+$is+10+8;                                                
                                                    $read_short_start_tempie =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    undef %read_short_start_tmp;
                                                    my $test_dot = $read_short_start_tempie =~ tr/\./\./;
                                                    
                                                    if ($test_dot > 1)
                                                    {
                                                        %read_short_start_tmp = build_partial3c $read_short_start_tempie;
                                                        my $ff = '0';
                                                        my %hash_read_short_start;
                                                        undef %hash_read_short_start;
                                                        foreach my $read_short_start_tempie2 (keys %read_short_start_tmp)
                                                        {
                                                            $ff = '0';
                                                            while ($ff < (length($read_short_start_tempie2)-length($match_pair_middle2)-1))
                                                            {
                                                                my $read_short_start_part = substr $read_short_start_tempie2, $ff, length($match_pair_middle2);
                                                                                               
                                                                $hash_read_short_start{$read_short_start_part} = undef;
                                                                $ff++;
                                                            }
                                                        }
                                                        foreach my $line (keys %hash_read_short_start)
                                                        {
                                                            my $found_seq = '0';
                                                            my $star_first = substr $line, 0, 1;
                                                            if ($star_first eq "*")
                                                            {
                                                                $line = substr $line, 1;
                                                            }  
                        
                                                            $found_seq = $match_pair_middle2 =~ s/$line/$line/;
                                                            if ($found_seq > 0)
                                                            {
                                                                $extension_match = "";
                                                                goto SKIP_BACK;
                                                            }
                                                        }
                                                        $extension_match = "NOOO";
                                                    }
                                                    else
                                                    {
                                                        %read_short_start_tmp = build_partial3b $read_short_start_tempie;
                                                        foreach my $line (keys %read_short_start_tmp)
                                                        {
                                                            my $found_seq = '0';
                        
                                                            $found_seq = $line =~ s/$match_pair_middle2/$match_pair_middle2/;
                                                            if ($found_seq > 0)
                                                            {
                                                                $extension_match = "";
                                                                goto SKIP_BACK;
                                                            }
                                                        }
                                                        $extension_match = "NOOO";
                                                    }                              
                                                }    
SKIP_BACK:
                                                my $nuc_exlude = "test";
                                                if ($yuyu_option_A_back eq "A" || $yuyu_option_C_back eq "C" || $yuyu_option_T_back eq "T" || $yuyu_option_G_back eq "G")
                                                {
                                                    $nuc_exlude = substr $extension, 0, 1;
                                                }
                                                if ($extension_match ne "NOOO" && $extension ne " " && $extension ne "" && $nuc_exlude ne $yuyu_option_A_back && $nuc_exlude ne $yuyu_option_C_back && $nuc_exlude ne $yuyu_option_T_back && $nuc_exlude ne $yuyu_option_G_back)
                                                {
                                                    $read_ex++;
                                                    push @matches1, $id_match.",".$extension.",".$id_pair_match.",".$match.",".$match_pair;

                                                    if ($extension ne " " && $extension ne "")
                                                    {
                                                        $extensions2{$extension} = $id_match;
                                                        $extensions2b{$id_match} = $extension;
                                                        push @extensions2, $extension;
                                                    }               
                                                }                                                                                         
        
                                            }                                
                }                
            }
            %extensions = (%extensions1, %extensions2);
            %extensionsb = (%extensions1b, %extensions2b);
            @extensions = (@extensions1, @extensions2);
            @matches = (@matches1, @matches2);
           
            my $ext = '0';
            my $ext_total_back = '0';
            foreach (@extensions)
            {
                $ext++;
            }
            $ext_total_back = $ext;
            
            if ($y > $startprint2)
            {
                print OUTPUT5 "\n".$read_count ." READ_COUNT_BACK\n";
                print OUTPUT5 $read_ex ." READ_EX_BACK\n";
                print OUTPUT5 $ext ." EXTENSIONS_BACK\n";
            }               

            if ($ext < 8 && $extra_regex eq "" && $last_chance_back ne "yes" && $indel_split_back eq '0' && $ext < $mmbr && $repetitive_detect_back eq "")
            {
                
                undef %extensions1;
                undef %extensions2;
                undef %extensions;
                undef @extensions1;
                undef @extensions2;
                undef @extensions;
                undef @matches;
                undef @matches1;
                undef @matches2;

                $extra_regex = "yes";
                goto REGEX_BACK;
            }
            
            $extra_regex = "";
            undef %SNR_length;
            
            if ($y > $startprint && $print_log eq '2')
            {
                foreach my $matches (@matches)
                {
                    my @matchesb;
                    undef @matchesb;
                    @matchesb = split /,/, $matches;
                    
                    print OUTPUT5 $matchesb[0].",".$matchesb[1]."\n";
                }               
            }
            my $id_original = $id;

            my $l = '0';
            my $best_extension = "";
            my $SNP = "";
            my $A_SNP = '0';
            my $C_SNP = '0';
            my $T_SNP = '0';
            my $G_SNP = '0';
            my $position_SNP = $position_back;
            my $pos_SNP = '0';
            
            my $A_SNP2 = '0';
            my $C_SNP2 = '0';
            my $T_SNP2 = '0';
            my $G_SNP2 = '0';
            my $position_SNP2 = $position_back;
            my $pos_SNP2 = '0';
            
            my $A_SNP3 = '0';
            my $C_SNP3 = '0';
            my $T_SNP3 = '0';
            my $G_SNP3 = '0';
            my $position_SNP3 = $position_back;
            my $pos_SNP3 = '0';
            
            if ($SNR_read_back ne "" && $SNR_read_back2 ne "" && $last_chance_back ne "yes" && $ext < 5 && $indel_split_back eq '0')
            {
                goto AFTER_EXT_BACK;
            }
            undef @extensions_group1;
            undef @extensions_group2;
            undef %extensions_group1;
            undef %extensions_group2;
            undef @extensions_group3;
            undef @extensions_group4;
            undef %extensions_group3;
            undef %extensions_group4;
            
            $no_SNR = "";
            undef %extensions_backup;
            undef @extensions_backup;
            my %extensions_backup = %extensions;
            my @extensions_backup = @extensions;
            
SPLIT_BACK:

            if ($split eq "yes_back")
            {
                %extensions = %extensions_group2;
                @extensions = @extensions_group2;
                $split = "yes2_back";
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            elsif ($split eq "yes2_back")
            {
                @extensions = @extensions_group3;
                %extensions = %extensions_group3;
                $split = "yes3_back";
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            elsif ($split eq "yes3_back")
            {
                @extensions = @extensions_group4;
                %extensions = %extensions_group4;
                $split = "yes4_back";
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            elsif ($split eq "yes4_back")
            {
                @extensions = @extensions_group1;
                %extensions = %extensions_group1;
                $split = "yes5_back";
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            $id = $id_original;
            $position_back = $position_back{$id};                   
            $l  = '0';
            $best_extension = "";          
            $SNP   = "";
            $A_SNP = '0';
            $C_SNP = '0';
            $T_SNP = '0';
            $G_SNP = '0';
            $position_SNP = $position_back;
            $pos_SNP = '0';
            
            $A_SNP2 = '0';
            $C_SNP2 = '0';
            $T_SNP2 = '0';
            $G_SNP2 = '0';
            $position_SNP2 = $position_back;
            $pos_SNP2 = '0';
            
            $A_SNP3 = '0';
            $C_SNP3 = '0';
            $T_SNP3 = '0';
            $G_SNP3 = '0';
            $position_SNP3 = $position_back;
            $pos_SNP3 = '0';
            my %SNR_count_back;
            my %extensions_new;
            my @extensions_new;
            undef %SNR_count_back;
            undef %extensions_new;
            undef @extensions_new;
            my $SNR_test = "";
            my $most_SNR = '0';
            $most_SNR2 = '0';
                      
            if ($SNR_read_back ne "" && $split eq "" && $SNR_read_back2 ne "")
            {     
                $SNR_test = "yes2_back";
                if ($SNR_read_back eq "yes")
                {
                    $SNR_test = "yes2_back";
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $SNR_test." SNR_TEST_BACK\n";
                    }
                    my $G = '0';
                    my $G2 = '0';
                    my $no_SNR1 = "";
                    my $second_round;
SNR1_BACK:          foreach my $extensions (@extensions)
                    { 
                        my @chars = split//, $extensions;
                        my $e = '0';
                        my $check_only_SNR = $extensions =~ s/$SNR_nucleo_back/$SNR_nucleo_back/g;                       
                        if (length($extensions) eq $check_only_SNR)
                        {
                            $G2++;
                        }
                        $G++;
                        if ($second_round eq "")
                        {
                            $no_SNR1 = "yes";
                        }
                        while ($SNR_nucleo_back eq $chars[$e] || $no_SNR1 eq "")
                        {
                            if ($SNR_nucleo_back ne $chars[$e])
                            {
                                $no_SNR1 = "yes";
                            }
                            $e++;
                        }
                        if ($e < length($extensions))
                        {                      
                            $SNR_count_back{$extensions} = $e;
                            $SNR_length{$e} .= exists $SNR_length{$e} ? ",$extensions" : $extensions;
                        }
                        $no_SNR1 = "";
                    }
                    my $SNR_length_count2 = '0';
                    my $SNR_length_reads = "";
                    my $first = "";
                    foreach my $SNR_length (keys %SNR_length)
                    {
                        my $SNR_length_count = $SNR_length{$SNR_length} =~ tr/,/,/;
                        if ($SNR_length_count > $SNR_length_count2)
                        {
                            $SNR_length_count2 = $SNR_length_count;
                            $SNR_length_reads = $SNR_length{$SNR_length};
                            $most_SNR = $SNR_length;
                            $first = $SNR_length_count;
                        }
                    }
                    my @SNR_length = split/,/, $SNR_length_reads;
                    my $repetitive_test = substr $read_short_start2, 0, 10;
                    my $SNR_checkSNR = $repetitive_test =~ s/$SNR_nucleo_back/$SNR_nucleo_back/g;
                    
                    $SNR_length_count2 = '0';
                    foreach my $SNR_length (keys %SNR_length)
                    {
                        my $SNR_length_count = $SNR_length{$SNR_length} =~ tr/,/,/;
                        if ($SNR_length_count > $SNR_length_count2 && $SNR_length_count ne $first)
                        {
                            $SNR_length_count2 = $SNR_length_count;
                            $most_SNR2 = $SNR_length;
                        }
                    }
print OUTPUT5 $first." 1st ".$SNR_length_count2." 2d ".$G." G\n";
                    if ($first < 0.8*$G && $second_round eq "" && $first > 0.35*$G && $SNR_length_count2 > 0.35*($G-$G2))
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "SPLIT_SNR\n";
                            $no_SNR = "yes";
                        }
                        goto NUCLEO0_BACK;
                    }
                    elsif ($first < 0.8*$G && $second_round eq "")
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 $most_SNR." MOST_SNR SECOND ROUND\n";
                        }
                        $second_round = "yes";
                        $no_SNR1 = "";
                        $most_SNR = '0';
                        $G = '0';
                        undef %SNR_count_back;
                        undef %SNR_length;
                        goto SNR1_BACK;
                    }
                    else
                    {
                        foreach my $SNRie (@SNR_length)
                        {
                            if (exists($extensions{$SNRie}))
                            {
                                $extensions_new{$SNRie} = $extensions{$SNRie};
                                push @extensions_new, $SNRie;
                            }
                        }
                        %extensions = %extensions_new;
                        @extensions = @extensions_new;
                    }            
                }
            
                if ($SNR_back{$id} eq "yes2_double_back")
                {
                    $SNR_test = "yes2_double_back";
                    foreach my $extensions (@extensions)
                    { 
                        my @chars = split//, $extensions;
                        my $e = '0';
                        
                        if ($SNR_nucleo_back eq $chars[$e].$chars[$e+1] )
                        {
                            while ($SNR_nucleo_back eq $chars[$e].$chars[$e+1])
                            {
                                my $tempie = reverse $extensions;
                                chop $tempie;
                                chop $tempie;
                                $extensions = reverse $tempie;
                                $e++;
                                $e++;
                            }
                        }
                        else
                        {  
                            while ($SNR_nucleo_back eq $chars[$e+1].$chars[$e])
                            {
                                my $tempie = reverse $extensions;
                                chop $tempie;
                                chop $tempie;
                                $extensions = reverse $tempie;
                                $e++;
                                $e++;
                            }
                        }
                        $extensions_new{$extensions} = $extensions{$extensions};
                        push @extensions_new, $extensions;
                        if ($e < length($extensions))
                        {                      
                            $SNR_count_back{$extensions} = $e;
                        }
                    }
                    %extensions = %extensions_new;
                    @extensions = @extensions_new;
                }
                delete $SNR_back{$id};
            }
NUCLEO0_BACK:            
            if ($SNR_read_back ne "")
            {
                $ext = '0';
                foreach (@extensions)
                {
                    $ext++;
                }
            }
            my $extra_l = '0';

NUCLEO_BACK: while ($l < $read_length - ($overlap+$left-1) + $extra_l)
             {
                my $A = '0';
                my $C = '0';
                my $T = '0';
                my $G = '0';
                
                if ($SNR_read_back ne "" && $l > 0)
                {
                    my $best_extension_tmp8 = $best_extension;
                    $best_extension_tmp8 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    my $last_nuc = substr $best_extension_tmp8, -1;
                    my $arrSize1 = @extensions;

                    if ($last_nuc ne "." && $arrSize1 > 4)
                    {
                        my @extensions_tmp;
                        undef @extensions_tmp;
                        foreach my $extensions (@extensions)
                        {
                            my @chars = split //, $extensions;
                            if ($chars[$l-1] eq $last_nuc || length($extensions) < $l)
                            {
                                push @extensions_tmp, $extensions;
                            }
                        }
                        
                        my $arrSize2 = @extensions_tmp;
                        if ($arrSize1 ne $arrSize2)
                        {
                            undef @extensions;
                            @extensions = @extensions_tmp;
                            my $best_extension_dot = $best_extension_tmp8 =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\./\./;
              
                            if ($best_extension_dot > 0)
                            {
                                $l = 0;
                                $best_extension = "";
                                $SNP = "";
                                goto NUCLEO_BACK;
                            }
                        }
                    }               
                }             
                
                foreach my $extensions (@extensions)
                {                                
                    my @chars = split//, $extensions;
                    
                    if ($chars[$l] eq "A")
                    {
                        $A++;
                    }
                    elsif ($chars[$l] eq "C")
                    {
                        $C++;
                    }
                    elsif ($chars[$l] eq "T")
                    {
                        $T++;
                    }
                    elsif ($chars[$l] eq "G")
                    {
                        $G++;
                    }
                }
                my $c = '2.8';               
                if ($ext > 22 && $SNR_read eq "")
                {
                    $c = '3.7';
                }
                if ($ext > 6 && ($type ne "chloro" || $extensions_before ne "") && $SNR_read_back eq "")
                {
                    $c = '5';
                }
                if ($ext > 22 && $SNR_read eq "" && $type ne "chloro")
                {
                    $c = '6.5';
                }
                if ($ext > 38 && $SNR_read eq "" && $type ne "chloro")
                {
                    $c = '8.4';
                }
                if ($repetitive_detect_back ne "" && $ext < 23 && $SNR_read_back eq "")
                {
                   $c = '7';
                }
                if ($repetitive_detect_back ne "" && $ext >= 23 && $SNR_read_back eq "")
                {
                   $c = '9';
                }
                if ($repetitive_detect_back2 eq "yes" || ($repetitive_detect_back ne "" && $ext > 120))
                {
                   $c = '15';
                }
                if ($type eq "mito_plant")
                {
                    $c += 2;
                }
                my $q = '2';
                my $z = '1';
                my $v = '6';
                my $s = '3';
                if ($split ne "" && $indel_split_back eq 0)
                {
                    $v = '10';
                    $s = '2';
                    $z = '0';
                }
                if ($A > ($C + $T + $G)*$c && (($A > $s && ($ext)/($A + $T + $G + $C) < $q) || ($A > $z && $l < $v && ($C + $T + $G) eq 0 && ($ext)/($A + $T + $G + $C) < $q)))
                {
                    $best_extension = $best_extension."A";
                }
                elsif ($C > ($A + $T + $G)*$c && (($C > $s && ($ext)/($A + $T + $G + $C) < $q) || ($C > $z && $l < $v && ($A + $T + $G) eq 0 && ($ext)/($A + $T + $G + $C) < $q)))
                {
                    $best_extension = $best_extension."C";  
                }
                elsif ($T > ($A + $C + $G)*$c && (($T > $s && ($ext)/($A + $T + $G + $C) < $q) || ($T > $z && $l < $v && ($A + $C + $G) eq 0 && ($ext)/($A + $T + $G + $C) < $q)))
                {
                    $best_extension = $best_extension."T"; 
                }
                elsif ($G > ($C + $T + $A)*$c && (($G > $s && ($ext)/($A + $T + $G + $C) < $q) || ($G > $z && $l < $v && ($C + $T + $A) eq 0 && ($ext)/($A + $T + $G + $C) < $q)))
                {
                    $best_extension = $best_extension."G";  
                }
                elsif (($SNP_active_back eq "yes" || ($SNR_read_back ne "" && $l > 0) || ($extensions_before eq "yes" && $ext_before ne "yes")) && $SNP eq "" && ($A + $T + $G + $C) > 4 && (($l < 15 && $split eq "") || ($l < 11 && $split ne "")) && ($ext)/($A + $T + $G + $C) < $q) 
                {
                    delete $SNP_active_back{$id};
                    $SNP = "yes_back";
                    $A_SNP = $A;
                    $C_SNP = $C;
                    $T_SNP = $T;
                    $G_SNP = $G;
                    $position_SNP += $l;
                    $pos_SNP = $l;
                    
                    my $IUPAC = IUPAC($A,$C,$T,$G);
                    $best_extension = $best_extension.$IUPAC;
                }
                elsif ($SNP eq "yes_back" && ($A + $T + $G + $C) > 4  && $l < 15)
                {
                    $SNP = "yes2_back";
                    $A_SNP2 = $A;
                    $C_SNP2 = $C;
                    $T_SNP2 = $T;
                    $G_SNP2 = $G;
                    $position_SNP2 += $l;
                    $pos_SNP2 = $l;            
                    
                    my $IUPAC = IUPAC($A,$C,$T,$G);
                    $best_extension = $best_extension.$IUPAC;
                }
                elsif ($SNP eq "yes2_back" && ($A + $T + $G + $C) > 4  && $l < 15)
                {
                    $SNP = "yes3_back";
                    $A_SNP3 = $A;
                    $C_SNP3 = $C;
                    $T_SNP3 = $T;
                    $G_SNP3 = $G;
                    $position_SNP3 += $l;
                    $pos_SNP3 = $l;      
                    
                    my $IUPAC = IUPAC($A,$C,$T,$G);
                    $best_extension = $best_extension.$IUPAC;
                }
                elsif ($SNP eq "yes3_back" && ($pos_SNP ne 0 || ($pos_SNP3 > $pos_SNP+12 && $l > 11 ) || ($extensions_before eq "yes" && $l > 9)) && $split eq "")
                {
                    $SNP = "yes4_back";
                    my $g = $l;
                    my $pos_SNP_tmp = $pos_SNP;
                    if ($pos_SNP3 > $pos_SNP+12)
                    {
                        $pos_SNP_tmp = $pos_SNP3;
                    }
                    if ($extensions_before ne "yes")
                    {
                        while ($g > $pos_SNP_tmp)
                        {                                         
                            chop($best_extension);
                            $g--;
                        }
                    }                    
                    last  NUCLEO_BACK;
                }
                
                elsif (($SNP eq "yes3_back" && $pos_SNP eq 0 && $l <= 15) || ($indel_split_skip_back ne "yes" && $l eq 0 && $ext > 4))
                {
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $SNP." SNP\n";
                    }
                    if ($SNR_test ne "" && $no_SNR ne "yes")
                    {
                        $l = '0';
                        $best_extension = "";
                        $no_SNR = "yes";
                        %extensions = %extensions_backup;
                        @extensions = @extensions_backup;
                        
                        goto SPLIT_BACK;
                    }
                    if ($SNP eq "")
                    {
                        $A_SNP = $A;
                        $C_SNP = $C;
                        $T_SNP = $T;
                        $G_SNP = $G;
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 $best_extension." BEST_EXTENSIONll\n";
                            print OUTPUT5 $A_SNP." A\n";
                            print OUTPUT5 $C_SNP." C\n";
                            print OUTPUT5 $T_SNP." T\n";
                            print OUTPUT5 $G_SNP." G\n";
                        }
                    }
                    if ($split ne "")
                    {
                        $read_new = $read_new1;
                        $best_extension = "";
                        delete $seed{$id};
                        $seed{$id} = $read_new;    
                        
                        delete $last_chance_back{$id};
                        $noback = "stop";
                        $noback{$id} = "stop";
                        if ($noforward ne "stop")
                        {
                            goto SEED;
                        }
                        else
                        {
                            delete $seed{$id};
                            $id = $id_original;
                            $split = "";
                            goto AFTER_EXT_BACK;
                        }
                    }
                    $best_extension = "";
                    $split = "yes_back";
                    
                    undef @firstSNP_max;

                    if ($A_SNP >= ($C_SNP+$A_SNP+$T_SNP+$G_SNP)*0.035 && $A_SNP > 0)
                    {
                        if (exists($yuyu_option_back{$id.'A'}))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "YUYU_A_EXISTS\n";
                            }
                        }
                        else
                        {
                            push @firstSNP_max, "A";
                        }  
                    }
                    if ($C_SNP >= ($C_SNP+$A_SNP+$T_SNP+$G_SNP)*0.035 && $C_SNP > 0)
                    {
                        if (exists($yuyu_option_back{$id.'C'}))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "YUYU_C_EXISTS\n";
                            }
                        }
                        else
                        {
                            push @firstSNP_max, "C";
                        }  
                    }
                    if ($T_SNP >= ($C_SNP+$A_SNP+$T_SNP+$G_SNP)*0.035 && $T_SNP > 0)
                    {
                        if (exists($yuyu_option_back{$id.'T'}))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "YUYU_T_EXISTS\n";
                            }
                        }
                        else
                        {
                            push @firstSNP_max, "T";
                        }  
                    }
                    if ($G_SNP >= ($C_SNP+$A_SNP+$T_SNP+$G_SNP)*0.035 && $G_SNP > 0)
                    {
                        if (exists($yuyu_option_back{$id.'G'}))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "YUYU_G_EXISTS\n";
                            }
                        }
                        else
                        {
                            push @firstSNP_max, "G";
                        }  
                    }
                    
                    $count_split = @firstSNP_max;
                    my $count_split_tmp = '0';
                    
                    if ($count_split eq '2')
                    {
                        $delete_third = "yes_back";
                        $delete_second = "yes_back";
                    }
                    if ($count_split eq '3')
                    {
                        $delete_third = "yes_back";
                    }
                    
                    foreach my $firstSNP_max (@firstSNP_max)
                    {                                            
                        foreach my $extensions_tmp (@extensions)
                        {                       
                            my @chars = split//, $extensions_tmp;
                            if ($chars[0] eq $firstSNP_max && $count_split_tmp eq '0')
                            {
                                $extensions_group1{$extensions_tmp} = $extensions{$extensions_tmp};
                                push @extensions_group1, $extensions_tmp;
                            }
                            elsif ($chars[0] eq $firstSNP_max && $count_split_tmp eq '1')
                            {
                                $extensions_group2{$extensions_tmp} = $extensions{$extensions_tmp};
                                push @extensions_group2, $extensions_tmp;
                            }
                            elsif ($chars[0] eq $firstSNP_max && $count_split_tmp eq '2')
                            {
                                $extensions_group3{$extensions_tmp} = $extensions{$extensions_tmp};
                                push @extensions_group3, $extensions_tmp;
                            }
                            elsif ($chars[0] eq $firstSNP_max && $count_split_tmp eq '3')
                            {
                                $extensions_group4{$extensions_tmp} = $extensions{$extensions_tmp};
                                push @extensions_group4, $extensions_tmp;
                            }
                        }
                        $count_split_tmp++;
                    }
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $count_split." COUNT_SPLIT\n";
                    }
                    goto SPLIT_BACK;
                }  
                else
                {  
                    last  NUCLEO_BACK;
                }
                $l++;
            }
            my $last_nucleo = substr $best_extension, -1;
            $last_nucleo =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
            while ($last_nucleo eq '.')
            {
                chop($best_extension);
                $last_nucleo = substr $best_extension, -1;
                $last_nucleo =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
            }
            if ($split eq "yes2_back" || ($split eq "yes3_back" && $count_split > 2) || ($split eq "yes4_back" && $count_split > 3))
            {  
                if ($split eq "yes2_back")
                {
                    $best_extension2 = $best_extension;
                }
                elsif ($split eq "yes3_back")
                {
                    $best_extension3 = $best_extension;
                }
                elsif ($split eq "yes4_back")
                {
                    $best_extension4 = $best_extension;
                }
                my $best_extension_tmp = $best_extension;
                $best_extension_tmp =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;

                if ($y > $startprint2)
                {
                    if ($split eq "yes2_back")
                    {
                        print OUTPUT5 "GROUP2\n";
                        foreach my $extensions_tmp (@extensions_group2)
                        {  
                            print OUTPUT5 $extensions_tmp."\n";
                        }
                        print OUTPUT5 $best_extension_tmp." BEST_EXTENSION_BACK2\n\n";
                    }
                    elsif ($split eq "yes3_back")
                    {
                        print OUTPUT5 "GROUP3\n";
                        foreach my $extensions_tmp (@extensions_group3)
                        {  
                            print OUTPUT5 $extensions_tmp."\n";                        
                        }
                        print OUTPUT5 $best_extension_tmp." BEST_EXTENSION_BACK3\n\n";
                    }
                    elsif ($split eq "yes4_back")
                    {
                         print OUTPUT5 "GROUP4\n";
                        foreach my $extensions_tmp (@extensions_group4)
                        {  
                            print OUTPUT5 $extensions_tmp."\n";                        
                        }
                        print OUTPUT5 $best_extension_tmp." BEST_EXTENSION_BACK4\n\n";
                    }    
                }           
                if ((length($best_extension2) < 3 || (length($best_extension2) < 6 && $ext > 15)) && $repetitive_detect_back eq "" && $before_back eq "yesss")
                {
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "\nDELETE BEST EXTENSION_BACK 2\n\n";
                    }
                    $delete_first = "yes_back";
                    goto SPLIT_BACK;                       
                }
                if ($type eq "chloroplst" && length($best_extension2) > 20)
                {
                    my $best_extension2_reverse2 = $best_extension2;
                    $best_extension2_reverse2 =~ tr/ATCG/TAGC/;
                    my $best_extension2_reverse = reverse($best_extension2_reverse2);
                    $best_extension2_reverse =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    
                    my $read_cp = $read;
                    my $read_start_rev_tmp = reverse($read_start);
                    $read_start_rev_tmp =~ tr/ATCG/TAGC/;
                    $read_start_rev_tmp =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    
                    if (length($read) > $genome_range_low)
                    {
                    }                 
                    my $found_seq_cp2 = $read_cp =~ s/$best_extension2_reverse/$best_extension2_reverse/;
                    my $found_seq_cp4 = $read_cp =~ s/$read_start_rev_tmp/$read_start_rev_tmp/;
                    
                    if ($found_seq_cp2 > 0 && $found_seq_cp4 > 0)
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nDELETE BEST EXTENSION 2 (CP)\n\n";
                        }
                        $delete_first = "yes_back";
                        goto SPLIT_BACK;                               
                    }
                }
                my $end_SNR = substr $read_start, 0,4;
                $end_SNR =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                my $GGGG = $end_SNR =~ tr/G/G/;
                my $TTTT = $end_SNR =~ tr/T/T/;
                my $CCCC = $end_SNR =~ tr/C/C/;
                my $AAAA = $end_SNR =~ tr/A/A/;
                if (($GGGG eq '4' || $TTTT eq '4' || $CCCC eq '4' || $AAAA eq '4') && $before_back eq "yesss")
                {
                    $GGGG = $best_extension2 =~ tr/G/G/;
                    $TTTT = $best_extension2 =~ tr/T/T/;
                    $CCCC = $best_extension2 =~ tr/C/C/;
                    $AAAA = $best_extension2 =~ tr/A/A/;
                    if ($GGGG eq length($best_extension2) || $TTTT eq length($best_extension2) || $CCCC eq length($best_extension2) || $AAAA eq length($best_extension2))
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nDELETE BEST EXTENSION 2 (SNR)\n\n";
                        }
                        $delete_first = "yes_back";
                        goto SPLIT_BACK;    
                    }             
                }
                if (length($best_extension2) > 9 && $before_back eq "yesss")
                {
                    my $end_tmp = substr $read_start, 0, -10;
                    if (length($best_extension2) < 15)
                    {
                        $end_tmp = substr $read_start, 0, -length($best_extension2)+5;
                    }
                    $end_tmp = reverse($best_extension2).$end_tmp;
                    $end_tmp =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    my $s = '0';
                    my $foundit = "";
                    while ($s < length($end_tmp)-$overlap)
                    {
                        my $end_tmp_d = substr $end_tmp, -($s+$overlap), $overlap;
                        if ($containX_short_start2 > 0)
                        {
                            my $star = $end_tmp_d =~ tr/\*//;
                    
                            $end_tmp_d = substr $end_tmp, -($s+$overlap+($star*2)), $overlap+($star*2);
                            my $star2 = $end_tmp_d =~ tr/\*//;                                                
                            while ($star2 > $star)
                            {
                                $end_tmp_d = substr $end_tmp, -($s+$overlap+($star*2)), $overlap+($star*2);
                                $star = $star2;
                                $star2 = $end_tmp_d =~ tr/\*//;
                            }
                        }    
                        
                        my %end_tmp_d = build_partial3b $end_tmp_d, "reverse";
                        foreach my $end_tmp_d (keys %end_tmp_d)
                        {
                            if (exists($hash2b{$end_tmp_d}))
                            {
                                $foundit = "yes_back";
                            }
                            elsif (exists($hash2c{$end_tmp_d}))
                            {
                                $foundit = "yes_back";
                            }
                        }
                        $s++;
                    }
                    if ($foundit ne "yes_back")
                    {
                        print OUTPUT5 "\nDELETE BEST EXTENSION 2 (no reverse match)\n\n";
                        $delete_first = "yes_back";
                        goto SPLIT_BACK;    
                    }  
                }              
                if ($split eq "yes2_back")
                {
                    delete $seed{$id};
                    $id      = "b".$position_back."_".$id;
                    $id_split1 = $id;
                }
                elsif ($split eq "yes3_back")
                {
                    $id      = "bb".$position_back."_".$id;
                    $id_split2 = $id;
                }
                elsif ($split eq "yes4_back")
                {
                    $id      = "bbb".$position_back."_".$id;
                    $id_split3 = $id;
                }
                
                $position_back{$id} = $position_back;
                $position{$id} = $position;
                $contig_count{$id} = $contig_count{$id_original};
                $position_adjust{$id} = $position_adjust{$id_original};
        
                my $count_contig_tmp = $contig_count;
                while ($count_contig_tmp > 0)
                {
                    $contig_gap_min{$id."_".$count_contig_tmp} = $contig_gap_min{$id_original."_".$count_contig_tmp};
                    $contig_gap_max{$id."_".$count_contig_tmp} = $contig_gap_max{$id_original."_".$count_contig_tmp};
                    $count_contig_tmp--;    
                }
                if (exists($old_id{$id_original}))
                {
                    $old_id{$id} = $old_id{$id_original};
                }
                if (exists($old_rep{$id_original}))
                {
                    $old_rep{$id} = $old_rep{$id_original};
                }
                if (exists($old_rep_old{$id_original}))
                {
                    $old_rep_old{$id} = $old_rep_old{$id_original};
                }
                if (exists($noforward{$id_original}))
                {
                    $noforward{$id} = $noforward;
                }
                if (exists($seed_split{$id_original}))
                {
                    $seed_split{$id} = undef;
                }
                if (exists($nosecond{$id_original}))
                {
                    $nosecond{$id} = undef;
                }
                if (exists($before_shorter_skip{$id_original}))
                {
                    $before_shorter_skip{$id} = undef;
                }
                if (exists($jump_rep{$id_original}))
                {
                    $jump_rep{$id} = undef;
                }
                if (exists($rep_return{$id_original}))
                {
                    $rep_return{$id} = $rep_return{$id_original};
                }
                if (exists($rep_return_back{$id_original}))
                {
                    $rep_return_back{$id} = $rep_return_back{$id_original};
                }
            }
            elsif ($split eq "yes5_back")
            {
                $best_extension1 = $best_extension;
                
                if ($y > $startprint2)
                {
                    print OUTPUT5 "GROUP1\n";
                    foreach my $extensions_tmp (@extensions_group1)
                    {  
                        print OUTPUT5 $extensions_tmp."\n";
                    }
                    print OUTPUT5 $best_extension1." BEST_EXTENSION_BACK1\n\n";
                }
                if ((length($best_extension1) < 3 || (length($best_extension1) < 6 && $ext > 15)) && $before_back eq "yessss" && $repetitive_detect_back eq "")
                {
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "\nDELETE BEST EXTENSION 1\n\n";
                    }
                    if ($delete_first eq "yes_back" && $last_chance_back eq "yes")
                    {
                        delete $seed{$id};                                         
                        $delete_first = "yes2_back";
                        $noback = "stop";
                        $noback{$id} = "stop";
                        goto FINISH;                
                    }
                    elsif ($delete_first eq "yes_back")
                    {
                        $best_extension = "";
                        goto AFTER_EXT_BACK;
                    }
                    else
                    {
                        goto SEED;  
                    }                      
                }
                if ($type eq "chlorop" && length($best_extension1) > 20)
                {
                    my $best_extension1_reverse2 = $best_extension1;
                    $best_extension1_reverse2 =~ tr/ATCG/TAGC/;
                    my $best_extension1_reverse = reverse($best_extension1_reverse2);
                    $best_extension1_reverse =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    
                    my $read_cp = $read;
                    my $read_start_rev_tmp = reverse($read_start);
                    $read_start_rev_tmp =~ tr/ATCG/TAGC/;
                    if (length($read) > $genome_range_low)
                    {
                        $read_cp = substr $read, 0, -5000;
                        $read_cp =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    }
                    
                    my $found_seq_cp2 = $read_cp =~ s/$best_extension1_reverse/$best_extension1_reverse/;
                    my $found_seq_cp4 = $read_cp =~ s/$read_start_rev_tmp/$read_start_rev_tmp/;
                    
                    if ($found_seq_cp2 > 0 && $found_seq_cp4 > 0)
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nDELETE BEST EXTENSION 1 (CP)\n\n";
                        }
                        if ($delete_first eq "yes_back")
                        {
                            delete $seed{$id};                                         
                            $delete_first = "yes2_back";
                            $noback = "stop";
                            $noback{$id} = "stop";
                            goto FINISH;                
                        }
                        goto SEED;                             
                    }
                }
                my $end_SNR = substr $read_start,0, 4;
                $end_SNR =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                my $GGGG = $end_SNR =~ tr/G/G/;
                my $TTTT = $end_SNR =~ tr/T/T/;
                my $CCCC = $end_SNR =~ tr/C/C/;
                my $AAAA = $end_SNR =~ tr/A/A/;
                if (($GGGG eq '4' || $TTTT eq '4' || $CCCC eq '4' || $AAAA eq '4') && $before_back eq "yesss")
                {
                    $GGGG = $best_extension1 =~ tr/G/G/;
                    $TTTT = $best_extension1 =~ tr/T/T/;
                    $CCCC = $best_extension1 =~ tr/C/C/;
                    $AAAA = $best_extension1 =~ tr/A/A/;
                    if ($GGGG eq length($best_extension1) || $TTTT eq length($best_extension1) || $CCCC eq length($best_extension1) || $AAAA eq length($best_extension1))
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nDELETE BEST EXTENSION 1 (SNR)\n\n";
                        }
                        if ($delete_first eq "yes_back")
                        {
                            delete $seed{$id};                                         
                            $delete_first = "yes2_back";
                            $noback = "stop";
                            $noback{$id} = "stop";
                            goto FINISH;                
                        }
                        goto SEED;    
                    }             
                }
                if (length($best_extension1) > 9 && $before_back eq "yessss")
                {
                    my $end_tmp = substr $read_start,0, -10;
                    if (length($best_extension1) < 15)
                    {
                        $end_tmp = substr $read_start, 0, -length($best_extension1)+5;
                    }
                    $end_tmp = reverse($best_extension1).$end_tmp;
                    $end_tmp =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                    my $s = '0';
                    my $foundit = "";
                    while ($s < length($end_tmp)-$overlap)
                    {
                        my $end_tmp_d = substr $end_tmp, -($s+$overlap), $overlap;
                        
                        if ($containX_short_start2 > 0)
                        {
                            my $star = $end_tmp_d =~ tr/\*//;
                    
                            $end_tmp_d = substr $end_tmp, -($s+$overlap+($star*2)), $overlap+($star*2);
                            my $star2 = $end_tmp_d =~ tr/\*//;                                                
                            while ($star2 > $star)
                            {
                                $end_tmp_d = substr $end_tmp, -($s+$overlap+($star*2)), $overlap+($star*2);
                                $star = $star2;
                                $star2 = $end_tmp_d =~ tr/\*//;
                            }
                        }    
                        my %end_tmp_d = build_partial3b $end_tmp_d, "reverse";
                        foreach my $end_tmp_d (keys %end_tmp_d)
                        {
                            if (exists($hash2b{$end_tmp_d}))
                            {
                                $foundit = "yes_back";
                            }
                            elsif (exists($hash2c{$end_tmp_d}))
                            {
                                $foundit = "yes_back";
                            }
                        }
                        $s++;
                    }
                    if ($foundit ne "yes_back")
                    {
                        print OUTPUT5 "\nDELETE BEST EXTENSION 1 (no reverse match)\n\n";
                        if ($delete_first eq "yes_back")
                        {
                            delete $seed{$id};                                         
                            $delete_first = "yes2_back";
                            $noback = "stop";
                            $noback{$id} = "stop";
                            goto FINISH;                
                        }                     
                        goto SEED;    
                    }  
                }
                if ($SNP_active_back eq "yes" && ($before_back eq "yes" || $extensions_before eq "yes") && $count_split eq '2' && $SNR_read_back eq "yes")
                {
                    my @chars;
                    my @chars2;
                    undef @chars;
                    undef @chars2;
                    my $p = '0';
                    my $amatch = '0';
                    my $nomatch = '0';
                    my $best_extension_short = "";
                    my $best_extension_long = "";
                    
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "CHECK_DELETION_BACK\n";
                    } 
                    
                    if (length($best_extension1) >= length($best_extension2))
                    {
                        @chars = split //, $best_extension2;
                        @chars2 = split //, $best_extension1;
                        $best_extension_short = $best_extension2;
                        $best_extension_long = $best_extension1;
                    }
                    elsif (length($best_extension2) > length($best_extension1))
                    {
                        @chars = split //, $best_extension1;
                        @chars2 = split //, $best_extension2;
                        $best_extension_short = $best_extension1;
                        $best_extension_long = $best_extension2;
                    }
                    while ($p < length($best_extension_short) && $p < 4)
                    {
                        my $i = '1';
                        $amatch = '0';
                        $nomatch = '0';
                        while ($i < @chars-$p)
                        {
                            if ($chars[$i+$p] eq $chars2[$i-1])
                            {
                                $amatch++;   
                            }
                            elsif ($chars[$i+$p] eq '.')
                            {                            
                            }
                            else
                            {                            
                                $nomatch++;
                            }
                            $i++;
                        }
                        $p++;
                        if ($amatch > ($nomatch*8) && (($amatch > 2 && $p < 3) || $amatch > 7))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "DELETION\n";
                            } 
                            $deletion_back = "yes";
                            delete $seed{$id};
                            $read_new = $read_new1;
                                                  
                            my $indel = substr $best_extension_short,0, $p;
                            $indel =~ s/A/*A/g;
                            $indel =~ s/T/*T/g;
                            $indel =~ s/G/*G/g;
                            $indel =~ s/C/*C/g;
                            
                            my $after_indel;
                            if ($nomatch > 0)
                            {
                                my $after_indel1 = substr $best_extension_short, $p;
                                $after_indel1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                my @chars3 = split //, $after_indel1;
                                my $f = '0';
                                
INDEL0_BACK:                    while ($f < @chars3)
                                {
                                    if ($chars3[$f] eq $chars2[$f])
                                    {
                                        $after_indel .= $chars3[$f];
                                    }
                                    elsif ($f < 20)
                                    {
                                        $after_indel .= ".";
                                    }
                                    else
                                    {
                                        last INDEL0_BACK;
                                    }
                                    $f++;
                                }
                            }
                            else
                            {
                                my $after_indel1 = substr $best_extension_short, $p;
                                $after_indel1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                my @chars3 = split //, $after_indel1;
                                my $f = '0';
                                
INDEL1_BACK:                    while ($f < @chars3)
                                {
                                    if ($chars3[$f] eq $chars2[$f])
                                    {
                                        $after_indel .= $chars3[$f];
                                    }
                                    else
                                    {
                                        last INDEL1_BACK;
                                    }
                                    $f++;
                                }
                            }
                            $best_extension = $indel.$after_indel;
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 reverse($best_extension)." BEST_EXTENSION_BACK(1)\n";
                            }
                            
                            goto INDEL_BACK;
                        }  
                    }
                    $p = '0';
                    while ($p < length($best_extension_short) && $p < 4)
                    {
                        my $i = '1';
                        $amatch = '0';
                        $nomatch = '0';
                        while ($i <= @chars && $i <= @chars2-$p)
                        {
                            if ($chars[$i-1] eq $chars2[$i+$p])
                            {
                                $amatch++;
                            }
                            elsif ($chars[$i-1] eq '.')
                            {                            
                            }
                            else
                            {
                                $nomatch++;
                            }
                            $i++;
                        }
                        $p++;
                        if ($amatch > ($nomatch*8) && (($amatch > 2 && $p < 3) || $amatch > 7))
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "DELETION\n";
                            } 
                            $deletion_back = "yes";
                            delete $seed{$id};
                            $read_new = $read_new1;
                           
                            my $indel = substr $best_extension_long,0, $p;
                            $indel =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                            my $indel2 = $indel;
                            $indel =~ s/A/*A/g;
                            $indel =~ s/T/*T/g;
                            $indel =~ s/G/*G/g;
                            $indel =~ s/C/*C/g;
                            
                            my $after_indel;
                            if ($nomatch > 0)
                            {
                                my $after_indel1 = substr $best_extension_long, $p;
                                $after_indel1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                my @chars3 = split //, $after_indel1;
                                my $f = '0';
                                
INDEL2_BACK:                    while ($f < @chars3)
                                {
                                    if ($chars3[$f] eq $chars[$f])
                                    {
                                        $after_indel .= $chars3[$f];
                                    }
                                    elsif ($f < 20)
                                    {
                                        $after_indel .= ".";
                                    }
                                    else
                                    {
                                        last INDEL2_BACK;
                                    }
                                    $f++;
                                }
                            }
                            else
                            {
                                my $after_indel1 = substr $best_extension_long, $p;
                                $after_indel1 =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                my @chars3 = split //, $after_indel1;
                                my $f = '0';
                                
INDEL3_BACK:                    while ($f < @chars3)
                                {
                                    if ($chars3[$f] eq $chars[$f])
                                    {
                                        $after_indel .= $chars3[$f];
                                    }
                                    else
                                    {
                                        last INDEL3_BACK;
                                    }
                                    $f++;
                                }
                            }
                            $best_extension = $indel.$after_indel;
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 reverse($best_extension)." BEST_EXTENSION_BACK(2)\n";
                            }       
                            goto INDEL_BACK;
                        }
                    }
                }
                if (length($best_extension1) > 4 && length($best_extension2) > 4 && $cp_input ne "" && $type eq "mito_plant")
                {
                   my $p = -$overlap;
                    while ($p > (-$overlap*2))
                    {
                        my $ref_part2 = substr $read_short_start2, $p, $overlap;
                        my %ref_part = build_partial3b $ref_part2, "";
                        foreach my $ref_part (keys %ref_part)
                        {
                            if (exists($hashref{$ref_part}))
                            {
                                my $ref_loc = -(($p + $overlap)/10) + 6;
                                $ref_loc =~ s/\..*$//;
                                
                                my $ref_id3 = $hashref{$ref_part};       print OUTPUT5 $ref_id3." REFID!!".$p."\n\n";
                                my $ref_id2 = substr $ref_id3, 1;
                                my @ref_id3 = split /,/,$ref_id2;
                           
                                foreach my $ref_id (@ref_id3)
                                {
                                    my $prev_loc1 = $ref_id -6;
                                    my $prev_loc2 = $ref_id -12;
                                    my $ref_part_prev = substr $read_short_start2, $p-($overlap*2), ($overlap*2)+15;
                                    
                                    if (exists($hashref2{$prev_loc1}))
                                    {
                                        my $ref_check_prev1 = $hashref2{$prev_loc1};
                                        
                                        if (exists($hashref2{$prev_loc2}))
                                        {
                                            my $ref_check_prev2 = $hashref2{$prev_loc2};
                                            my @ref_check_prev1 = build_partialb $ref_check_prev1;
                                            
                                            foreach my $ref_check_prev1 (@ref_check_prev1)
                                            {
                                                my $found_ref1 = $ref_part_prev =~ s/$ref_check_prev1/$ref_check_prev1/;
                                                if ($found_ref1 > 0)
                                                {
                                                    my @ref_check_prev2 = build_partialb $ref_check_prev2;
                                                    
                                                    foreach my $ref_check_prev2 (@ref_check_prev2)
                                                    {
                                                        my $found_ref2 = $ref_part_prev =~ s/$ref_check_prev2/$ref_check_prev2/;
                       
                                                        if ($found_ref2 > 0)
                                                        {
                                                            if (exists($hashref2{$ref_id+$ref_loc}))
                                                            {
                                                                my $ref_check = $hashref2{$ref_id+$ref_loc};      
                                                                my $best_extension1_part = substr $best_extension1,0 ,30;
                                                                my $best_extension2_part = substr $best_extension2,0 ,30;
                                                                my @ref1 =  build_partial $best_extension1_part;
                                                                
                                                                my $ref_part3 = substr $read_short_start2, $p;
                                                                my $star_ref = $ref_part3 =~ tr/\*//;
                                                                my $loc_start = $overlap*3-($ref_loc*10)-($p + ($overlap*3))-($star_ref*2);
                                                                
                                                                foreach my $best_extension1_part2 (@ref1)
                                                                {                                  
                                                                    if ($ref_check =~ m/.{$loc_start}$best_extension1_part2.*/)
                                                                    {    
                                                                        print OUTPUT5 "REFERNCE_GUIDED\n";
                                                                        print OUTPUT5 $best_extension1." BEST_EXTENSION1\n\n";
                                                                        delete $seed{$id};
                                                                        goto INDEL_BACK;
                                                                    }
                                                                }
                                                                my @ref2 =  build_partial $best_extension2_part;
                                                                foreach my $best_extension2_part2 (@ref2)
                                                                {
                                                                    if ($ref_check =~ m/.{$loc_start}$best_extension2_part2.*/)
                                                                    {
                                                                        print OUTPUT5 "REFERNCE_GUIDED\n";
                                                                        print OUTPUT5 $best_extension2." BEST_EXTENSION2\n\n";
                                                                        delete $seed{$id};
                                                                        $best_extension = $best_extension2;
                                                                        goto INDEL_BACK;
                                                                    }
                                                                }
                                                            }  
                                                        }
                                                    }
                                                }
                                            }
                                        }    
                                    }
                                }                        
                            }
                        }
                        $p--;
                    }
                }
INDEL_BACKa:  
                if ($last_chance_back ne "yes" && (($SNP_active_back eq "" && $before_back eq "yes") || $extensions_before eq "yes") && $indel_split_skip_back ne "yes" && $ext_total_back > 15 && (($delete_first ne "yes_back" || $count_split ne "2") || (length($best_extension1) < 5 && length($best_extension2) < 5 && length($best_extension3) < 5 && length($best_extension4) < 5)) && $indel_split_back < ($read_length-25 -$overlap))
                {
                    delete $seed{$id};
                    delete $seed{$id_split1};
                    delete $seed{$id_split2};
                    delete $seed{$id_split3};
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "\nDELETE SPLIT_BACK\n\n";
                    }
                }
                elsif ((($before_back eq "yes" && $SNP_active_back eq "") || $extensions_before eq "yes") && (($delete_first ne "yes_back" && $count_split eq "2") || (length($best_extension1) < 5 && length($best_extension2) < 5 && length($best_extension3) < 5 && length($best_extension4) < 5)))
                {
                    delete $seed{$id};
                    delete $seed{$id_split1};
                    delete $seed{$id_split2};
                    delete $seed{$id_split3};
                    $SNP_active_back{$id_original} = "yes";
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "\nDELETE SPLIT_BACK2\n\n";
                    }
                }
                elsif (($before_back eq "yes" || $extensions_before eq "yes") && (($delete_first ne "yes_back" || $count_split ne "2") || (length($best_extension1) < 5 && length($best_extension2) < 5 && length($best_extension3) < 5 && length($best_extension4) < 5)))
                {
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $delete_first." DELETE_FIRST\n";
                        print OUTPUT5 "STOP BACK\n";
                    }
                    $read_new = $read_new1;
                    $best_extension = "";
                    delete $seed{$id};
                    delete $seed{$id_split1};
                    delete $seed{$id_split2};
                    delete $seed{$id_split3};
                    $seed{$id} = $read_new;    
                    
                    delete $last_chance_back{$id};
                    $noback = "stop";
                    $noback{$id} = "stop";
                    if ($noforward ne "stop")
                    {
                        goto SEED;
                    }
                    else
                    {
                        delete $seed{$id};
                        $id             = $id_original;
                        $split          = "";
                        goto AFTER_EXT_BACK;
                    }
                }
INDEL_BACK:     
                $id                   = $id_original;
                $split                = "";
                $split_forward        = "yes";
                

                if ($jump_rep_back eq "yes" || ($extensions_before ne "yes" && $indel_split_skip_back ne "yes" && ($delete_first eq "" || $delete_second eq "" || $delete_third eq "") && $before_back ne "yes"))
                {                   
                    my $before_shorter = "";
                    my $overhangb = 1+($read_length/40);
                    my $count_one = '0';
                    my $count_two = '0';
                    while (($count_one < 2 || $count_two < 2) && $overhangb < 15)
                    {
                        $overhangb++;
                        foreach my $one (@extensions_group1)
                        {
                            if (length($one) < $overhangb && length($one) > 0)
                            {
                                $count_one++;
                            }
                        }
                        foreach my $two (@extensions_group2)
                        {
                            if (length($two) < $overhangb && length($two) > 0)
                            {
                                $count_two++;
                            }
                        }
                    }
                    my $overhang = sprintf("%.0f", $overhangb);
                    delete $jump_rep_back{$id};
BEFORE_BACK:
                    my $s = '0';
                    
                    if ($before_shorter eq "yes")
                    {
                        $overhang += 5;
                    }
                    $before_shorter = "";
                   
                    print OUTPUT5 $overhang." OVERHANG_BACK\n";
                    my @extensions_yuyu;
                    undef @extensions_yuyu;
                    my %extensions_yuyu;
                    my %before1F;
                    my %before2F;
                    my %before3F;
                    my %before4F;
                    my %before1B;
                    my %before2B;
                    my %before3B;
                    my %before4B;
                    undef %repetitive_pair;
                    
                    my $skip_overhang = "";
                    if ($AT_rich_before_back eq "yes")
                    {
                        my $before_split_tmp = substr $read_short_start2, 0, $read_length-$overhang;
                        my $A_rich_test = $before_split_tmp =~ tr/AX\.//;
                        my $T_rich_test = $before_split_tmp =~ tr/TX\.//;
                        my $G_rich_test = $before_split_tmp =~ tr/GX\.//;
                        my $C_rich_test = $before_split_tmp =~ tr/CX\.//;
                        if ($A_rich_test > ($read_length-$overhang-3) || $T_rich_test > ($read_length-$overhang-3) || $G_rich_test > ($read_length-$overhang-3) || $C_rich_test > ($read_length-$overhang-3))
                        {
                            $skip_overhang = "yes";
                            $overhang -= 5;
                        }
                    }
                    $s = '0';
                    my $before_split = substr $read_short_start2, $read_length-$left-1-$overhang-$overlap, $overhang+$overlap;
                    my $star0 = '0';
                    my $star0b = '0';
                    my $star0b2 = '0';
                    if ($containX_short_start2 > 0)
                    {
                        my $before_split2 = substr $read_short_start2, 0, $read_length-$left-1-$overlap-$overhang;
                        my $before_split3 = substr $read_short_start2, $read_length-$left-1-$overlap-$overhang, $overhang+$overlap;
                        $star0 = $before_split2 =~ tr/\*//;
                        $star0b = $before_split3 =~ tr/\*//;
                        if ($star0 > 0 || $star0b > 0)
                        {
                            $before_split = substr $read_short_start2, $read_length-$left-1-$overhang-$overlap+($star0*2), $overhang+$overlap+($star0)+($star0b*2);
                            my $star0c = $before_split =~ tr/\*//;                                                
                            while ($star0c > $star0b)
                            {
                                $before_split = substr $read_short_start2, $read_length-$left-1-$overhang-$overlap+($star0*2), $overhang+$overlap+($star0)+($star0b*2);
                                $star0b = $star0c;
                                $star0c = $before_split =~ tr/\*//;
                            }
                        }
                    }
                      
                  print OUTPUT5 $before_split." BEFORE_SPLIT\n"; 
                    while ($s <= length($before_split)-$overlap-($star0b*2))
                    {
                        my $line_tmp_reverse = substr $before_split, $s, $overlap;
                        
                        if ($star0 > 0 || $star0b > 0)
                        {
                            my $star = $line_tmp_reverse =~ tr/\*//;
                            $line_tmp_reverse = substr $before_split, $s, $overlap+($star*2);
                            my $star2 = $line_tmp_reverse =~ tr/\*//;                                                
                            while ($star2 > $star)
                            {
                                $line_tmp_reverse = substr $before_split, $s, $overlap+($star*2)+(($star2-$star)*2);
                                $star = $star2;
                                $star2 = $line_tmp_reverse =~ tr/\*//;
                            }
                            my $before_split2_tmp = substr $read_short_start2, 0, $read_length-$left-1-$overhang-$overlap+($star0*2)+(length($before_split)-$s-$overlap);
                            $star0b2 = $before_split2_tmp =~ tr/\*//;      
                        }
                        
                        my %line_tmp = build_partial3b $line_tmp_reverse, "reverse";
                        foreach my $line_reverse (keys %line_tmp)
                        {
                            $line_reverse =~ tr/ACTG/TGAC/;
                            my $line = reverse($line_reverse);
                            if (exists($hash2b{$line}))
                            {                       
                                my $search = $hash2b{$line};
                                $search = substr $search, 1;
                                my @search = split /,/,$search;
                                
                                foreach my $search (@search)
                                {                             
                                    my $search_tmp = substr $search, 0, -1;
                                    my $search_end = substr $search, -1;
                                    if (exists($hash{$search_tmp}))
                                    {
                                        my @search_tmp = split /,/,$hash{$search_tmp};
                                        my $found;
                                        if ($search_end eq "1")
                                        {
                                            $found = $search_tmp[0];
                                        }
                                        elsif ($search_end eq "2")
                                        {
                                            $found = $search_tmp[1];
                                        }
                                        if ($encrypt eq "yes")
                                        {
                                            $found = decrypt $found;
                                        }
                                        my $found_new;
                                        my $first_10 = substr $found, $read_length-1-$overhang+$s-10+$star0, 10;
                                        my $first_10b = substr $found, $read_length-1-$overhang+$s-10, 10;
                                        my $first_nuc;
                                        $first_10 =~ tr/ACTG/TGAC/;
                                        $first_10b =~ tr/ACTG/TGAC/;
                                        my $first_10r = reverse($first_10);
                                        my $first_10rb = reverse($first_10b);
                                        my $first_10_read_start = substr $read_start, 0, 10;
                                        my $check_first10 = $first_10r =~ s/$first_10_read_start/$first_10_read_start/;
                                        my $check_first10b = $first_10rb =~ s/$first_10_read_start/$first_10_read_start/;
      
                                        if ($check_first10 > 0)
                                        {
                                            $found_new = substr $found, 0, $read_length-1-$overhang+$s+$star0,"";
                                            $first_nuc = substr $found, 0, 1;
                                        }
                                        elsif ($check_first10b > 0)
                                        {
                                            $found_new = substr $found, 0, $read_length-1-$overhang+$s,"";
                                            $first_nuc = substr $found, 0, 1;
                                        }
                                        if ($check_first10 > 0 || $check_first10b > 0)
                                        {
                                            $first_nuc =~ tr/ACTG/TGAC/;
                                            $found_new =~ tr/ACTG/TGAC/;
                                            $found =~ tr/ACTG/TGAC/;
                                            my $first_nuc1 = substr $best_extension1, 0, 1;
                                            my $first_nuc2 = substr $best_extension2, 0, 1;
                                            my $first_nuc3 = substr $best_extension3, 0, 1;
                                            my $first_nuc4 = substr $best_extension4, 0, 1;
                                            
                                            my $found_reverse = reverse($found_new);
                                            my $extension_yuyu = $found;
                                            $extensions_yuyu{$search} = $extension_yuyu;
    
                                            if ($first_nuc eq $first_nuc1 && ($check_first10 > 0 || $check_first10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                $before1B{$search} = $found_reverse;
                                            }
                                            elsif ($first_nuc eq $first_nuc2 && ($check_first10 > 0 || $check_first10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                $before2B{$search} = $found_reverse;
                                            }
                                            elsif ($first_nuc eq $first_nuc3 && ($check_first10 > 0 || $check_first10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                $before3B{$search} = $found_reverse;
                                            }
                                            elsif ($first_nuc eq $first_nuc3 && ($check_first10 > 0 || $check_first10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                $before4B{$search} = $found_reverse;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        $s++;
                    }
                    $s = '0';
                    $before_split = substr $read_short_start2, $read_length-$right-1-$overhang-$overlap, $overhang+$overlap;
                    my $star3 = '0';
                    if ($containX_short_start2 > 0)
                    {
                        my $before_split2 = substr $read_short_start2, 0, $read_length-$right-1;
                        $star3 = $before_split2 =~ tr/\*//;
                        if ($star0 > 0)
                        {
                            $before_split = substr $read_short_start2, $read_length-$right-1-$overhang-$overlap+($star3*2), $overhang+$overlap+($star3*2);
                        }
                    }

                    while ($s <= length($before_split)-$overlap)
                    {
                        my $line_tmp = substr $before_split, $s, $overlap;
                        if ($star3 > 0)
                        {
                            my $star = $line_tmp =~ tr/\*//;
                            $line_tmp = substr $before_split, $s, $overlap+($star*2);
                            my $star2 = $line_tmp =~ tr/\*//;                                                
                            while ($star2 > $star)
                            {
                                $line_tmp = substr $before_split, $s, $overlap+($star*2)+(($star2-$star)*2);
                                $star = $star2;
                                $star2 = $line_tmp =~ tr/\*//;
                            }   
                        }
                        my %line_tmp = build_partial3b $line_tmp, "";
                        foreach my $line (keys %line_tmp)
                        {
                            if (exists($hash2c{$line}))
                            {                       
                                my $search = $hash2c{$line};
                                $search = substr $search, 1;
                                my @search = split /,/,$search;
                                                                            
                                foreach my $search (@search)
                                {                             
                                    my $search_tmp = substr $search, 0, -1;
                                    my $search_end = substr $search, -1;
                                    if (exists($hash{$search_tmp}))
                                    {
                                        my @search_tmp = split /,/,$hash{$search_tmp};
                                        my $found;
                                        if ($search_end eq "1")
                                        {
                                            $found = $search_tmp[0];
                                        }
                                        elsif ($search_end eq "2")
                                        {
                                            $found = $search_tmp[1];
                                        }
                                        if ($encrypt eq "yes")
                                        {
                                            $found = decrypt $found;
                                        }
                                        my $found_new;
                                        my $first_nuc;
                                        my $first_10 = substr $found, $overhang-$s+1, 10;
                                        my $first_10b = substr $found,  $overhang-$s+1-$star3, 10;
                                        
                                        my $first_10_read_start = substr $read_start, 0, 10;
                                        my $check_first10 = $first_10 =~ s/$first_10_read_start/$first_10_read_start/;
                                        my $check_first10b = $first_10b =~ s/$first_10_read_start/$first_10_read_start/;
                                        if ($check_first10 > 0)
                                        {
                                            $found_new = substr $found,  $overhang-$s+1, 100000, "";
                                            $first_nuc = substr $found, -1, 1;
                                        }
                                        elsif ($check_first10b > 0)
                                        {
                                            $found_new = substr $found, $overhang-$s+1-$star3, 100000, "";
                                            $first_nuc = substr $found,  -1, 1;
                                        }
                                        if ($check_first10 > 0 || $check_first10b > 0)
                                        {
                                            my $first_nuc1 = substr $best_extension1, 0, 1;
                                            my $first_nuc2 = substr $best_extension2, 0, 1;
                                            my $first_nuc3 = substr $best_extension3, 0, 1;
                                            my $first_nuc4 = substr $best_extension4, 0, 1;
                                            my $extension_yuyu = reverse($found);
                                            $extensions_yuyu{$search} = $extension_yuyu;
                               
                                            if ($first_nuc eq $first_nuc1 && ($check_first10 > 0 || $check_first10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                $filter_before1_pair{$search} = $found_new;
                                                $before1F{$search} = $found_new;
                                            }
                                            elsif ($first_nuc eq $first_nuc2 && ($check_first10 > 0 || $check_first10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                $filter_before2_pair{$search} = $found_new;
                                                $before2F{$search} = $found_new;
                                            }
                                            elsif ($first_nuc eq $first_nuc3 && ($check_first10 > 0 || $check_first10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                $filter_before3_pair{$search} = $found_new;
                                                $before3F{$search} = $found_new;
                                            }
                                            elsif ($first_nuc eq $first_nuc4 && ($check_first10 > 0 || $check_first10b > 0))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                }
                                                $filter_before4_pair{$search} = $found_new;
                                                $before4F{$search} = $found_new;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        $s++;
                    }
                    my %before_all1 = (%before1B, %before1F);
                    my %before_all2 = (%before2B, %before2F);
                    my %before_all3 = (%before3B, %before3F);
                    my %before_all4 = (%before4B, %before4F);


                    my $start_short_tmp = substr $read_short_start2, 0, $read_length+20;
                    my $star_first = substr $start_short_tmp, 0, 1;
                    if ($star_first eq "*")
                    {
                        $start_short_tmp = substr $read_short_start2, 1, $read_length+20;
                    }  
BEFORE_EXTRA_BACK:       
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $start_short_tmp." START_SHORT_TMP\n\n";
                    }
                    my $first_yuyu = "";
                    my $second_yuyu = "";
                    my $third_yuyu = "";
                    my $fourth_yuyu = "";
                    my $first_yuyu2 = '0';
                    my $second_yuyu2 = '0';
                    my $third_yuyu2 = '0';
                    my $fourth_yuyu2 = '0';
                    my @filter_dot_before1;
                    my @filter_dot_before2;
                    my @filter_dot_before3;
                    my @filter_dot_before4;
                    undef @filter_dot_before1;
                    undef @filter_dot_before2;
                    undef @filter_dot_before3;
                    undef @filter_dot_before4;
                    my @extensions_before;
                    my @extensions_before1;
                    my @extensions_before2;
                    my @extensions_before3;
                    my @extensions_before4;
                    undef @extensions_before;
                    undef @extensions_before1;
                    undef @extensions_before2;
                    undef @extensions_before3;
                    undef @extensions_before4;
                    
                    my $start_short_tmp_part = substr $start_short_tmp, 0, $read_length-$overhang-5-1; 
                    
                    foreach my $search (keys %before_all2)
                    {
                        my $overhang_tmp = $overhang;
                        my $yuyu0 = $before_all2{$search};
                        my $yuyu2 = $yuyu0;
                        if (length($yuyu2) >= $read_length-$overhang_tmp-5-1)
                        {
                            $yuyu2 =~ tr/N/\./;
                            my $yuyu2_tmp = $yuyu2;
                            my $check_yuyuy2 = $yuyu2_tmp =~ s/$start_short_tmp_part/$start_short_tmp_part/;

                            if ($check_yuyuy2 > 0)
                            {
                                $second_yuyu = "yes";
                                $second_yuyu2++;
                                my $start_short_tmp_part = substr $start_short_tmp, 0, length($yuyu2);
                                my $start_short_tmp_part2;
                                my $star = $start_short_tmp_part =~ tr/\*//;
                                my $check_yuyuy2b2 = '0';
                                if ($star > 0)
                                {
                                    $start_short_tmp_part = substr $start_short_tmp, 0, length($yuyu2)+($star*2);
                                    $start_short_tmp_part2 = substr $start_short_tmp, 0, length($yuyu2)+($star);
                                    my $yuyu2_tmp2 = $yuyu2;
                                    $check_yuyuy2b2 = $yuyu2_tmp2 =~ s/$start_short_tmp_part2/$start_short_tmp_part2/;
                                }      
                                my $yuyu2_tmp3 = $yuyu2;
                                my $check_yuyuy2b = $yuyu2_tmp3 =~ s/$start_short_tmp_part/$start_short_tmp_part/;
                                if ($check_yuyuy2b > 0 || $check_yuyuy2b2 > 0)
                                {
                                    push @filter_dot_before2, $yuyu0;
                                    if (exists($extensions_yuyu{$search}))
                                    {     
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $yuyu0." FOUND2   ";
                                            print OUTPUT5 $extensions_yuyu{$search}." EXTe2\n";
                                        }
                                        push @extensions_before, $extensions_yuyu{$search};
                                        push @extensions_before2, $extensions_yuyu{$search};
                                    }
                                    $filter_before2{$search} = undef;
                                }
                            }
                        }
                    }
                    foreach my $search (keys %before_all1)
                    {
                        my $overhang_tmp = $overhang;   
                        my $yuyu0 = $before_all1{$search};
                        my $yuyu1 = $yuyu0;
                        if (length($yuyu1) >= $read_length-$overhang_tmp-5-1)
                        {
                            $yuyu1 =~ tr/N/\./;
                            my $yuyu1_tmp = $yuyu1;
                            my $check_yuyuy = $yuyu1_tmp =~ s/$start_short_tmp_part/$start_short_tmp_part/;
                            
                            if ($check_yuyuy > 0)
                            {
                                $first_yuyu = "yes";
                                $first_yuyu2++;
                                my $start_short_tmp_part = substr $start_short_tmp, 0, length($yuyu1)-1;
                                my $star = $start_short_tmp_part =~ tr/\*//;
                                my $check_yuyuy1b2 = '0';
                                my $start_short_tmp_part2;
                                if ($star > 0)
                                {
                                    $start_short_tmp_part = substr $start_short_tmp, 0, length($yuyu1)-1+($star*2);
                                    $start_short_tmp_part2 = substr $start_short_tmp, 0, length($yuyu1)-1+($star);
                                    my $yuyu1_tmp2 = $yuyu1;
                                    $check_yuyuy1b2 = $yuyu1_tmp2 =~ s/$start_short_tmp_part2/$start_short_tmp_part2/;
                                }
                                my $yuyu1_tmp3 = $yuyu1;
                                my $check_yuyuy1b = $yuyu1_tmp3 =~ s/$start_short_tmp_part/$start_short_tmp_part/;
                                if ($check_yuyuy1b > 0 || $check_yuyuy1b2 > 0)
                                {
                                    push @filter_dot_before1, $yuyu0;
                                    if (exists($extensions_yuyu{$search}))
                                    {
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $yuyu0." FOUND1   ";
                                            print OUTPUT5 $extensions_yuyu{$search}." EXTe1\n";
                                        }
                                        push @extensions_before, $extensions_yuyu{$search};
                                        push @extensions_before1, $extensions_yuyu{$search};
                                    }
                                    $filter_before1{$search} = undef;
                                }
                            }
                        }
                    }
                    foreach my $search (keys %before_all3)
                    {
                        my $overhang_tmp = $overhang;
                        my $yuyu0 = $before_all3{$search};
                        my $yuyu3 = $yuyu0;
                        if (length($yuyu3) >= $read_length-$overhang_tmp-5-1)
                        {
                            $yuyu3 =~ tr/N/\./;
                            my $yuyu3_tmp = $yuyu3;
                            my $check_yuyuy = $yuyu3_tmp =~ s/$start_short_tmp_part/$start_short_tmp_part/;
                            
                            if ($check_yuyuy > 0)
                            {
                                $third_yuyu = "yes";
                                $third_yuyu2++;
                                my $start_short_tmp_part = substr $start_short_tmp, 0, length($yuyu3)-1;
                                my $star = $start_short_tmp_part =~ tr/\*//;
                                my $check_yuyuy3b2 = '0';
                                my $start_short_tmp_part2;
                                if ($star > 0)
                                {
                                    $start_short_tmp_part = substr $start_short_tmp, 0, length($yuyu3)-1+($star*2);
                                    $start_short_tmp_part2 = substr $start_short_tmp, 0, length($yuyu3)-1+($star);
                                    my $yuyu3_tmp2 = $yuyu3;
                                    $check_yuyuy3b2 = $yuyu3_tmp2 =~ s/$start_short_tmp_part2/$start_short_tmp_part2/;
                                }
                                my $yuyu3_tmp3 = $yuyu3;
                                my $check_yuyuy3b = $yuyu3_tmp3 =~ s/$start_short_tmp_part/$start_short_tmp_part/;
                                if ($check_yuyuy3b > 0 || $check_yuyuy3b2 > 0)
                                {
                                    push @filter_dot_before3, $yuyu0;
                                    if (exists($extensions_yuyu{$search}))
                                    {
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $yuyu0." FOUND3   ";
                                            print OUTPUT5 $extensions_yuyu{$search}." EXTe3\n";
                                        }
                                        push @extensions_before, $extensions_yuyu{$search};
                                        push @extensions_before3, $extensions_yuyu{$search};
                                    }
                                    $filter_before3{$search} = undef;
                                }
                            }
                        }
                    }
                    foreach my $search (keys %before_all4)
                    {
                        my $yuyu0 = $before_all4{$search};
                        my $yuyu4 = $yuyu0;
                        if (length($yuyu4) >= $read_length-$overhang-5-1)
                        {
                            $yuyu4 =~ tr/N/\./;
                            my $yuyu4_tmp = $yuyu4;
                            my $check_yuyuy = $yuyu4_tmp =~ s/$start_short_tmp_part/$start_short_tmp_part/;
                            
                            if ($check_yuyuy > 0)
                            {
                                $fourth_yuyu = "yes";
                                $fourth_yuyu2++;
                                my $start_short_tmp_part = substr $start_short_tmp, 0, length($yuyu4)-1;
                                my $star = $start_short_tmp_part =~ tr/\*//;
                                my $check_yuyuy4b2 = '0';
                                my $start_short_tmp_part2;
                                if ($star > 0)
                                {
                                    $start_short_tmp_part = substr $start_short_tmp, 0, length($yuyu4)-1+($star*2);
                                    $start_short_tmp_part2 = substr $start_short_tmp, 0, length($yuyu4)-1+($star);
                                    my $yuyu4_tmp2 = $yuyu4;
                                    $check_yuyuy4b2 = $yuyu4_tmp2 =~ s/$start_short_tmp_part2/$start_short_tmp_part2/;
                                }
                                my $yuyu4_tmp3 = $yuyu4;
                                my $check_yuyuy4b = $yuyu4_tmp3 =~ s/$start_short_tmp_part/$start_short_tmp_part/;
                                if ($check_yuyuy4b > 0 || $check_yuyuy4b2 > 0)
                                {
                                    push @filter_dot_before4, $yuyu0;
                                    if (exists($extensions_yuyu{$search}))
                                    {
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $yuyu0." FOUND4   ";
                                            print OUTPUT5 $extensions_yuyu{$search}." EXTe4\n";
                                        }
                                        push @extensions_before, $extensions_yuyu{$search};
                                        push @extensions_before4, $extensions_yuyu{$search};
                                    }
                                    $filter_before4{$search} = undef;
                                }
                            }
                        }
                    }

                    if ($y > $startprint2)
                    {
                        print OUTPUT5 $first_yuyu." FIRST_YUYU\n";
                        print OUTPUT5 $second_yuyu." SECOND_YUYU\n";
                        print OUTPUT5 $third_yuyu." THIRD_YUYU\n";
                        print OUTPUT5 $fourth_yuyu." FOURTH_YUYU\n";
                    }
                    
                    my $correction = '0';
                    my $count_all = '0';
                    $count_all = (@filter_dot_before1) + (@filter_dot_before2) + (@filter_dot_before3) + (@filter_dot_before4);
                    my $start_short_tmp2 = substr $read_short_start2, 0, $read_length-$overhang-5-1;
                    my $check_dot = $start_short_tmp2 =~ tr/\./\./;
                    if ($check_dot eq "")
                    {
                        $check_dot = '0'
                    }
                    my $check_star = $start_short_tmp2 =~ tr/\*/\*/;
                    my $check_star2 = "";
                    my $deletion = "";
print OUTPUT5 $check_star." CHECK_STAR\n";
print OUTPUT5 $count_all." COUNT_ALL\n";
                    if ($count_all > 2 && ($check_dot > 0 || $check_star > 0))
                    {
                        my @split_dot = split /\./, $start_short_tmp2;
                        my $size = '0';
                        my $size2 = '1';
                        my $length_total = '0';
                        my @filter_dot_before_all;
                        undef @filter_dot_before_all;
                        @filter_dot_before_all = (@filter_dot_before1, @filter_dot_before2, @filter_dot_before3, @filter_dot_before4);
                        
                        if ($check_star > 0)
                        {
                            my @deletion = split //, $start_short_tmp2;
                            my $found_star = "";
                            foreach my $nuc (@deletion)
                            {
                                if ($found_star eq "yes")
                                {
                                    $found_star = "";
                                    $deletion .= $nuc;
                                }
                                if ($nuc eq "*")
                                {
                                    $found_star = "yes";
                                }
                            }
                            print OUTPUT5 $deletion." NUC0\n";
                        }
                        foreach my $dot_split (@split_dot)
                        {
                           $size++; 
                        }
                        foreach my $dot_split (@split_dot)
                        {
                            if ($size2 ne $size)
                            {
                                $length_total += length($dot_split);
                                my $A = '0';
                                my $C = '0';
                                my $T = '0';
                                my $G = '0';
                                my $check_star2_part = "";
                                if ($check_star > 0)
                                {
                                    $check_star2_part = substr $start_short_tmp2, 0, $length_total;
                                    $check_star2 = $check_star2_part =~ tr/\*/\*/;
                                    print OUTPUT5 $check_star2." CHECK_STAR2\n";
                                }
                                foreach my $dot_before2 (@filter_dot_before_all)
                                {  
                                    my $dot = substr $dot_before2, $length_total, 1;
                                    if ($check_star2 > 0)
                                    {
                                        print OUTPUT5 $deletion." DEL_DOT\n";
                                        my $check_star2_part_tmp1 = $check_star2_part;
                                        $check_star2_part_tmp1 =~ tr/\*//d;
                                        my $dot_before2_tmp = $dot_before2;
                                        my $check_star1 = $dot_before2_tmp =~ s/$check_star2_part_tmp1/$check_star2_part_tmp1/;
                                        if ($check_star1 > 0)
                                        {
                                            $dot = substr $dot_before2, $length_total-$check_star2, 1;
                                        }
                                        else
                                        {
                                            $dot = substr $dot_before2, $length_total-($check_star2*2), 1;
                                        }
                                    }
                                    if ($dot eq "A")
                                    {
                                        $A++;
                                    }
                                    elsif ($dot eq "C")
                                    {
                                        $C++;
                                    }
                                    elsif ($dot eq "T")
                                    {
                                        $T++;
                                    }
                                    elsif ($dot eq "G")
                                    {
                                        $G++;
                                    } 
                                }
                                if (($A > 2 && ($C+$T+$G) eq 0) || ($A > ($C+$T+$G)*10 && $repetitive_detect eq ""))
                                {
                                    substr $read_new1, $length_total, 1, "A";
                                    print OUTPUT5 "DOT AAAAAAAAAAAA\n";
                                    $correction++;
                                }
                                if (($C > 2 && ($A+$T+$G) eq 0) || ($C > ($A+$T+$G)*10 && $repetitive_detect eq ""))
                                {
                                    substr $read_new1, $length_total, 1, "C";
                                    print OUTPUT5 "DOT CCCCCCCCCCCC\n";
                                    $correction++;
                                }
                                if (($T > 2 && ($C+$A+$G) eq 0) || ($T > ($C+$A+$G)*10 && $repetitive_detect eq ""))
                                {
                                    substr $read_new1, $length_total, 1, "T";
                                    print OUTPUT5 "DOT TTTTTTTTTTTT\n";
                                    $correction++;
                                }
                                if (($G > 2 && ($C+$T+$A) eq 0) || ($G > ($C+$T+$A)*10 && $repetitive_detect eq ""))
                                {
                                    substr $read_new1, $length_total, 1, "G";
                                    print OUTPUT5 "DOT GGGGGGGGGGGGG\n";
                                    $correction++;
                                }
                                $length_total += 1;
                            }
                            $size2++;
                        }
                    }
       
                    if ($first_yuyu eq "yes" && $second_yuyu ne "yes" && $third_yuyu ne "yes" && $fourth_yuyu ne "yes" && $first_yuyu2 > 2 && $correction eq $check_dot)
                    {
                        if (@extensions_before > 3)
                        {
                            $ext_before = "yes";
                            @extensions_before = @extensions_before1;
                            goto EXT_BEFORE_BACK;
                        }
                        delete $seed{$id_split1};
                        delete $seed{$id_split2};
                        delete $seed{$id_split3};
                        $best_extension = $best_extension1;
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nONLY FIRST\n\n";
                        }
                        delete $SNP_active_back{$id};
                        delete $before_back{$id};
                        delete $before_shorter_skip_back{$id};
                        goto AFTER_EXT_BACK;
                    }
                    elsif ($second_yuyu eq "yes" && $first_yuyu ne "yes" && $third_yuyu ne "yes" && $fourth_yuyu ne "yes" && $second_yuyu2 > 2 && $correction eq $check_dot)
                    {
                        if (@extensions_before > 3)
                        {
                            $ext_before = "yes";
                            @extensions_before = @extensions_before2;
                            goto EXT_BEFORE_BACK;
                        }
                        delete $seed{$id};
                        delete $seed{$id_split2};
                        delete $seed{$id_split3};
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nONLY SECOND\n\n";
                            print OUTPUT5 $best_extension2." BEST_E2\n\n";
                        }
                        delete $SNP_active_back{$id_split1};
                        delete $before_back{$id_split1};
                        delete $before_shorter_skip_back{$id_split1};
                        $id = $id_split1;
                        if ($correction > 0)
                        {
                            delete $seed{$id_split1};
                            $best_extension = $best_extension2;
                            goto AFTER_EXT_BACK;
                        }
                        else
                        {
                            $read_new = $seed{$id_split1};                       
                            goto FINISH;
                        }
                    }
                    elsif ($third_yuyu eq "yes" && $first_yuyu ne "yes" && $second_yuyu ne "yes" && $fourth_yuyu ne "yes" && $third_yuyu2 > 2 && $correction eq $check_dot) 
                    {
                        if (@extensions_before > 3)
                        {
                            $ext_before = "yes";
                            @extensions_before = @extensions_before3;
                            goto EXT_BEFORE_BACK;
                        }
                        delete $seed{$id};
                        delete $seed{$id_split3};
                        delete $seed{$id_split1};
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nONLY THIRD\n\n";
                        }
                        delete $SNP_active_back{$id_split2};
                        delete $before_back{$id_split2};
                        delete $before_shorter_skip_back{$id_split2};
                        $id = $id_split2;
                        if ($correction > 0)
                        {
                            delete $seed{$id_split2};
                            $best_extension = $best_extension3;
                            goto AFTER_EXT_BACK;
                        }
                        else
                        {
                            $read_new = $seed{$id_split2};                       
                            goto FINISH;
                        }
                    }
                    elsif ($fourth_yuyu eq "yes" && $first_yuyu ne "yes" && $third_yuyu ne "yes" && $second_yuyu ne "yes" && $fourth_yuyu2 > 2 && $correction eq $check_dot)
                    {
                        if (@extensions_before > 3)
                        {
                            $ext_before = "yes";
                            @extensions_before = @extensions_before4;
                            goto EXT_BEFORE_BACK;
                        }
                        delete $seed{$id};
                        delete $seed{$id_split2};
                        delete $seed{$id_split1};
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nONLY FOURTH\n\n";
                        }
                        delete $SNP_active_back{$id_split3};
                        delete $before_back{$id_split3};
                        delete $before_shorter_skip_back{$id_split3};
                        $id = $id_split3;
                        if ($correction > 0)
                        {
                            $best_extension = $best_extension4;
                            delete $seed{$id_split3};
                            goto AFTER_EXT_BACK;
                        }
                        else
                        {
                            $read_new = $seed{$id_split3};                       
                            goto FINISH;
                        }
                    }
                    elsif($count_all < 3 && $overhang < $read_length-($overlap+15) && $skip_overhang ne "yes")
                    {
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 "\nMAKE BEFORE SHORTER\n\n";
                        }
                        $before_shorter = "yes";
                        goto BEFORE_BACK;
                    }
                    else
                    {
                        my %count1234;
                        undef %count1234;
                        my %count1234b;
                        undef %count1234b;
                        my $count1 = '0';
                        my $count2 = '0';
                        my $count3 = '0';
                        my $count4 = '0';
                        my $count1_pair = '0';
                        my $count2_pair = '0';
                        my $count3_pair = '0';
                        my $count4_pair = '0';
                        foreach my $exb (keys %filter_before1)
                        {
                            my $search_tmp = substr $exb, 0, -1;
                            my $search_end = substr $exb, -1;
                            if (exists($hash{$search_tmp}))
                            {
                                my @search_tmp = split /,/,$hash{$search_tmp};
                                my $found;
                                my $found2;
                                if ($search_end eq "1")
                                {
                                    $found = $search_tmp[1];
                                    $found2 = $search_tmp[0];
                                }
                                elsif ($search_end eq "2")
                                {
                                    $found = $search_tmp[0];
                                    $found2 = $search_tmp[1];
                                }
                                if ($encrypt eq "yes")
                                {
                                    $found = decrypt $found;
                                }
                                $count1++;
                                my $found_reverse = reverse($found);
                                $found_reverse =~ tr/ACTG/TGAC/;
                                
                                if (exists($filter_before1{$exb}))
                                {
                                    if (exists($filter_before1_pair{$exb}))
                                    {
                                        $filter_before1{$exb} = $found_reverse;
                                        $count1_pair++;
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $found_reverse." 1\n";
                                        }
                                    }
                                    else
                                    {
                                        delete $filter_before1{$exb};
                                    }
                                }
                            }
                        }
                        foreach my $exb (keys %filter_before2)
                        {
                            my $search_tmp = substr $exb, 0, -1;
                            my $search_end = substr $exb, -1;
                            if (exists($hash{$search_tmp}))
                            {
                                my @search_tmp = split /,/,$hash{$search_tmp};
                                my $found;
                                my $found2;
                                if ($search_end eq "1")
                                {
                                    $found = $search_tmp[1];
                                    $found2 = $search_tmp[0];
                                }
                                elsif ($search_end eq "2")
                                {
                                    $found = $search_tmp[0];
                                    $found2 = $search_tmp[1];
                                }
                                if ($encrypt eq "yes")
                                {
                                    $found = decrypt $found;
                                }
                                my $found_reverse = reverse($found);
                                $found_reverse =~ tr/ACTG/TGAC/;
                                $count2++;
                                if (exists($filter_before2{$exb}))
                                {
                                    if (exists($filter_before2_pair{$exb}))
                                    {
                                        $filter_before2{$exb} = $found_reverse;
                                        $count2_pair++;
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $found_reverse." 2\n";
                                        }
                                    }
                                    else
                                    {
                                        delete $filter_before2{$exb};
                                    }
                                }
                            }
                        }
                        foreach my $exb (keys %filter_before3)
                        {
                            my $search_tmp = substr $exb, 0, -1;
                            my $search_end = substr $exb, -1;
                            if (exists($hash{$search_tmp}))
                            {
                                my @search_tmp = split /,/,$hash{$search_tmp};
                                my $found;
                                if ($search_end eq "1")
                                {
                                    $found = $search_tmp[1];
                                }
                                elsif ($search_end eq "2")
                                {
                                    $found = $search_tmp[0];
                                }
                                if ($encrypt eq "yes")
                                {
                                    $found = decrypt $found;
                                }
                                my $found_reverse = reverse($found);
                                $found_reverse =~ tr/ACTG/TGAC/;                      
                                $count3++;
                                if (exists($filter_before3{$exb}))
                                {
                                    if (exists($filter_before3_pair{$exb}))
                                    {
                                        $filter_before3{$exb} = $found_reverse;
                                        $count3_pair++;
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $found_reverse." 3\n";
                                        }
                                    }
                                    else
                                    {
                                        delete $filter_before3{$exb};
                                    }
                                }
                            }
                        }
                        foreach my $exb (keys %filter_before4)
                        {
                            my $search_tmp = substr $exb, 0, -1;
                            my $search_end = substr $exb, -1;
                            if (exists($hash{$search_tmp}))
                            {
                                my @search_tmp = split /,/,$hash{$search_tmp};
                                my $found;
                                if ($search_end eq "1")
                                {
                                    $found = $search_tmp[1];
                                }
                                elsif ($search_end eq "2")
                                {
                                    $found = $search_tmp[0];
                                }
                                if ($encrypt eq "yes")
                                {
                                    $found = decrypt $found;
                                }
                                my $found_reverse = reverse($found);
                                $found_reverse =~ tr/ACTG/TGAC/;
                                $count4++;
                                if (exists($filter_before4{$exb}))
                                {
                                    if (exists($filter_before4_pair{$exb}))
                                    {
                                        $filter_before4{$exb} = $found_reverse;
                                        $count4_pair++;
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $found_reverse." 4\n";
                                        }
                                    }
                                    else
                                    {
                                        delete $filter_before4{$exb};
                                    }
                                }
                            }
                        }
                        if ($y > $startprint2)
                        {
                            print OUTPUT5 $count1." COUNT1\n";
                            print OUTPUT5 $count2." COUNT2\n";
                            print OUTPUT5 $count3." COUNT3\n";
                            print OUTPUT5 $count4." COUNT4\n";
                        }
                        $count1234{'1'} = $count1_pair;
                        $count1234{'2'} = $count2_pair;
                        $count1234{'3'} = $count3_pair;
                        $count1234{'4'} = $count4_pair;
                        
                        my $morethan3 = '0';
                        my $difference = "";
                        
                        foreach my $count1234 (keys %count1234)
                        {
                            if ($count1234{$count1234} > 3)
                            {
                                $morethan3++;
                            }
                        }
                        my $best_extension1_tmp = $best_extension1;
                        my $best_extension2_tmp = $best_extension2;
                        my $SNR_check1 = $best_extension1_tmp =~ s/AAAAAAAA|CCCCCCCC|GGGGGGGG|TTTTTTTT//;
                        my $SNR_check2 = $best_extension2_tmp =~ s/AAAAAAAA|CCCCCCCC|GGGGGGGG|TTTTTTTT//;
                        my $GGGG2 = $best_extension2 =~ tr/G/G/;
                        my $TTTT2 = $best_extension2 =~ tr/T/T/;
                        my $CCCC2 = $best_extension2 =~ tr/C/C/;
                        my $AAAA2 = $best_extension2 =~ tr/A/A/;
                        if ($GGGG2 eq length($best_extension2) || $TTTT2 eq length($best_extension2) || $CCCC2 eq length($best_extension2) || $AAAA2 eq length($best_extension2))
                        {
                            $SNR_check2 = '1';
                        }
                        my $GGGG1 = $best_extension1 =~ tr/G/G/;
                        my $TTTT1 = $best_extension1 =~ tr/T/T/;
                        my $CCCC1 = $best_extension1 =~ tr/C/C/;
                        my $AAAA1 = $best_extension1 =~ tr/A/A/;
                        if ($GGGG1 eq length($best_extension1) || $TTTT1 eq length($best_extension1) || $CCCC1 eq length($best_extension1) || $AAAA1 eq length($best_extension1))
                        {
                            $SNR_check1 = '1';
                        }
                        if ((($count1 > 4 && ($count2+$count3+$count4) eq 0) || ($count1 > ($count2+$count3+$count4)*10 && ($count2+$count3+$count4) ne 0)) && $repetitive_detect_back2 ne "yes" && $correction eq $check_dot)
                        {
                            delete $seed{$id_split1};
                            delete $seed{$id_split2};
                            delete $seed{$id_split3};
                            $best_extension = $best_extension1;
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE9\n\n";
                            }
                            delete $SNP_active_back{$id};
                            delete $before_back{$id};
                            delete $before_shorter_skip_back{$id};
                            goto AFTER_EXT_BACK;
                        }
                        elsif ((($count2 > 4 && ($count1+$count3+$count4) eq 0) || ($count2 > ($count1+$count3+$count4)*10 && ($count1+$count3+$count4) ne 0)) && $repetitive_detect_back2 ne "yes" && $correction eq $check_dot)
                        {
                            delete $seed{$id};
                            delete $seed{$id_split1};
                            delete $seed{$id_split2};
                            delete $seed{$id_split3};
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE10\n\n";
                            }
                            delete $SNP_active_back{$id_split1};
                            delete $before_back{$id_split1};
                            delete $before_shorter_skip_back{$id_split1};
                            $id = $id_split1;

                            $best_extension = $best_extension2;
                            goto AFTER_EXT_BACK;
                        }
                        elsif ((($count3 > 4 && ($count1+$count2+$count4) eq 0) || ($count3 > ($count1+$count2+$count4)*10 && ($count1+$count2+$count4) ne 0)) && $repetitive_detect_back2 ne "yes" && $correction eq $check_dot)
                        {
                            delete $seed{$id};
                            delete $seed{$id_split3};
                            delete $seed{$id_split1};
                            delete $seed{$id_split2};
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE11\n\n";
                            }
                            delete $SNP_active_back{$id_split2};
                            delete $before_back{$id_split2};
                            delete $before_shorter_skip_back{$id_split2};
                            $id = $id_split2;
                            $best_extension = $best_extension3;
                            goto AFTER_EXT_BACK;
                        }
                        elsif ((($count4 > 4 && ($count1+$count3+$count2) eq 0) || ($count4 > ($count1+$count3+$count2)*10 && ($count1+$count3+$count2) ne 0)) && $repetitive_detect_back2 ne "yes" && $correction eq $check_dot)
                        {
                            if (@extensions_before > 3)
                            {
                                $ext_before = "yes";
                                goto EXT_BEFORE_BACK;
                            }
                            delete $seed{$id};
                            delete $seed{$id_split2};
                            delete $seed{$id_split1};
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE12\n\n";
                            }
                            delete $SNP_active_back{$id_split3};
                            delete $before_back{$id_split3};
                            delete $before_shorter_skip_back{$id_split3};
                            $id = $id_split3;

                            $best_extension = $best_extension4;
                            goto AFTER_EXT_BACK;
                        }
                        elsif(($morethan3 < 2 && $overhang < $read_length-($overlap*1.5)) || ($morethan3 > 1 && length($read) <= $insert_size && $overhang < 30) && $skip_overhang ne "yes")
                        {
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nMAKE BEFORE SHORTER2\n\n";
                            }
                            $before_shorter = "yes";
                            goto BEFORE_BACK;
                        }
                        
                        elsif ($morethan3 > 1 && length($read) > $insert_size && $last_chance_back ne "yes" && (($SNR_check1 eq "" && $SNR_check2 eq "") || $overhang < 30))
                        {
                            undef @extensions_before;
                            undef @extensions_before1;
                            undef @extensions_before2;
                            undef @extensions_before3;
                            undef @extensions_before4;
                            my $count1b = '0';
                            my $count2b = '0';
                            my $count3b = '0';
                            my $count4b = '0';
                                                                                                   
                            my $size = keys %read_short_start_tmp;
                            my $insert_range_before = 1.65;
                            if ($insert_range_shorter eq "yes")
                            {
                                $insert_range_before = 1.45;
                            }

                            my $F = ($insert_size-$read_length) - ((($insert_size*$insert_range_before)-$insert_size)/2) - $overhang;
                            if ($F < 0)
                            {
                                $F = '0';
                            }
                            my $read_short_start_tempie = substr $read, $F, ($insert_size*$insert_range_before)+$overhang;
                                                 
                            $read_short_start_tempie =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                               
                            undef %read_short_start_tmp;

                            %read_short_start_tmp = build_partial3c $read_short_start_tempie;
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 $read_short_start_tempie." READ_SHORT\n";
                            }
                        
                            my $ff = '0';
                            my %hash_read_short_start;
                            undef %hash_read_short_start;
                            foreach my $read_short_end_tempie (keys %read_short_end_tmp)
                            {
                                $ff = '0';
                                while ($ff < (length($read_short_start_tempie)-$read_length+$left+$right+2))
                                {
                                    my $read_short_start_part = substr $read_short_start_tempie, $ff, $read_length-$left-$right-1;
                                    my $star_start = substr $read_short_start_part, 0, 1;
                                    if ($star_start eq "*")
                                    {
                                        $read_short_start_part = substr $read_short_start_tempie, $ff+1, $read_length-$left-$right-1;
                                    }
                                                                   
                                    $hash_read_short_start{$read_short_start_part} = undef;
                                    $ff++;
                                }
                            }

                            foreach my $exb0 (keys %filter_before1)
                            { 
                                my $exb = $filter_before1{$exb0};
                                my $match_pair_middle = substr $exb, $right, $read_length-$left-$right;
FILTER_1_BACK:                  foreach my $line (keys %hash_read_short_start)
                                {
                                    my $found_seq = '0';
                                    my $match_pair_middle_tmp = $match_pair_middle;

                                    $found_seq = $match_pair_middle_tmp =~ s/.$line/$line/;
                                    if ($found_seq > 0)
                                    {
                                        if (exists($filter_before1_pair{$exb0}))
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $filter_before1_pair{$exb0}." FOUND1\n";
                                            }
                                        }
                                        if (exists($extensions_yuyu{$exb0}))
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $extensions_yuyu{$exb0}." EXT1\n";
                                            }
                                            push @extensions_before, $extensions_yuyu{$exb0};
                                            push @extensions_before1, $extensions_yuyu{$exb0};
                                        }
                                        $count1b++;
                                        last FILTER_1_BACK;
                                    }
                                }
                            }
                            foreach my $exb0 (keys %filter_before2)
                            { 
                                my $exb = $filter_before2{$exb0};
                                my $match_pair_middle = substr $exb, $right, $read_length-$left-$right;
FILTER_2_BACK:                  foreach my $line (keys %hash_read_short_start)
                                {
                                    my $found_seq = '0';
                                    my $match_pair_middle_tmp = $match_pair_middle;

                                    $found_seq = $match_pair_middle_tmp =~ s/.$line/$line/;
                                    if ($found_seq > 0)
                                    {
                                        if (exists($filter_before2_pair{$exb0}))
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $filter_before2_pair{$exb0}." FOUND2\n";
                                            }
                                            if (exists($extensions_yuyu{$exb0}))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 $extensions_yuyu{$exb0}." EXT2\n";
                                                }
                                                push @extensions_before, $extensions_yuyu{$exb0};
                                                push @extensions_before2, $extensions_yuyu{$exb0};
                                            }
                                        }
                                        $count2b++;
                                        last FILTER_2_BACK;
                                    }
                                }
                            }
                            foreach my $exb0 (keys %filter_before3)
                            { 
                                my $exb = $filter_before3{$exb0};
                                my $match_pair_middle = substr $exb, $right, $read_length-$left-$right;
FILTER_3_BACK:                  foreach my $line (keys %hash_read_short_start)
                                {
                                    my $found_seq = '0';
                                    my $match_pair_middle_tmp = $match_pair_middle;

                                    $found_seq = $match_pair_middle_tmp =~ s/.$line/$line/;
                                    if ($found_seq > 0)
                                    {
                                        if (exists($filter_before3_pair{$exb0}))
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $filter_before3_pair{$exb0}." FOUND3\n";
                                            }
                                            if (exists($extensions_yuyu{$exb0}))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 $extensions_yuyu{$exb0}." EXT3\n";
                                                }
                                                push @extensions_before, $extensions_yuyu{$exb0};
                                                push @extensions_before3, $extensions_yuyu{$exb0};
                                            }
                                        }
                                        $count3b++;
                                        last FILTER_3_BACK;
                                    }
                                }
                            }
                            foreach my $exb0 (keys %filter_before4)
                            { 
                                my $exb = $filter_before4{$exb0};
                                my $match_pair_middle = substr $exb, $right, $read_length-$left-$right;
FILTER_4_BACK:                  foreach my $line (keys %hash_read_short_start)
                                {
                                    my $found_seq = '0';
                                    my $match_pair_middle_tmp = $match_pair_middle;

                                    $found_seq = $match_pair_middle_tmp =~ s/.$line/$line/;
                                    if ($found_seq > 0)
                                    {
                                        if (exists($filter_before4_pair{$exb0}))
                                        {
                                            if ($y > $startprint2)
                                            {
                                                print OUTPUT5 $filter_before4_pair{$exb0}." FOUND4\n";
                                            }
                                            if (exists($extensions_yuyu{$exb0}))
                                            {
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 $extensions_yuyu{$exb0}." EXT4\n";
                                                }
                                                push @extensions_before, $extensions_yuyu{$exb0};
                                                push @extensions_before4, $extensions_yuyu{$exb0};
                                            }
                                        }
                                        $count4b++;
                                        last FILTER_4_BACK;
                                    }
                                }
                            }
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 $count1b." COUNT1B\n";
                                print OUTPUT5 $count2b." COUNT2B\n";
                                print OUTPUT5 $count3b." COUNT3B\n";
                                print OUTPUT5 $count4b." COUNT4B\n";
                            }
                            my $f = '7';
                            if ($repetitive_detect_back ne "")
                            {
                                $f = '10';
                            }
                            
                            $count1234b{'1'} = $count1b;
                            $count1234b{'2'} = $count2b;
                            $count1234b{'3'} = $count3b;
                            $count1234b{'4'} = $count4b;
                            
                            my $differenceb1 = "";
                            my $differenceb2 = "";
                            my $differenceb3 = "";
                            my $differenceb4 = "";
                            
                            foreach my $count1234b (keys %count1234b)
                            {   
                                if ($count1234b ne '1' && $count1234b{$count1234b} > 0 && $count1b > $f*$count1234b{$count1234b})
                                {
                                    $differenceb1 = "yes";
                                }
                                if ($count1234b ne '2' && $count1234b{$count1234b} > 0 && $count2b > $f*$count1234b{$count1234b})
                                {
                                    $differenceb2 = "yes";
                                }
                                if ($count1234b ne '3' && $count1234b{$count1234b} > 0 && $count3b > $f*$count1234b{$count1234b})
                                {
                                    $differenceb3 = "yes";
                                }
                                if ($count1234b ne '4' && $count1234b{$count1234b} > 0 && $count4b > $f*$count1234b{$count1234b})
                                {
                                    $differenceb4 = "yes";
                                }
                            }
                            
                            if (($count1b > 2 && ($count2b+$count3b+$count4b) eq '0') || $differenceb1 eq "yes")
                            {
                                if (@extensions_before > 15 && $repetitive_detect_back2 eq "")
                                {
                                    $ext_before = "yes";
                                    @extensions_before = @extensions_before1;
                                    goto EXT_BEFORE_BACK;
                                }
                                delete $seed{$id_split1};
                                delete $seed{$id_split2};
                                delete $seed{$id_split3};
                                $best_extension = $best_extension1;
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE5\n\n";
                                }
                                delete $SNP_active_back{$id};
                                delete $before_back{$id};
                                delete $before_shorter_skip_back{$id};
                                goto AFTER_EXT_BACK;
                            }
                            elsif (($count2b > 2 && ($count1b+$count3b+$count4b) eq '0') || $differenceb2 eq "yes")
                            {
                                if (@extensions_before > 15 && $repetitive_detect_back2 eq "")
                                {
                                    $ext_before = "yes";
                                    @extensions_before = @extensions_before2;
                                    goto EXT_BEFORE_BACK;
                                }
                                delete $seed{$id};
                                delete $seed{$id_split2};
                                delete $seed{$id_split3};
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE6\n\n";
                                }
                                delete $SNP_active_back{$id_split1};
                                delete $before_back{$id_split1};
                                delete $before_shorter_skip_back{$id_split1};
                                $id = $id_split1;
                                if ($correction > 0)
                                {
                                    delete $seed{$id_split1};
                                    $best_extension = $best_extension2;
                                    goto AFTER_EXT_BACK;
                                }
                                else
                                {
                                    $read_new = $seed{$id_split1};                       
                                    goto FINISH;
                                }
                            }
                            elsif (($count3b > 2 && ($count1b+$count2b+$count4b) eq '0') || $differenceb3 eq "yes")
                            {
                                if (@extensions_before > 15 && $repetitive_detect_back2 eq "")
                                {
                                    $ext_before = "yes";
                                    @extensions_before = @extensions_before3;
                                    goto EXT_BEFORE_BACK;
                                }
                                delete $seed{$id};
                                delete $seed{$id_split1};
                                delete $seed{$id_split3};
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE7\n\n";
                                }
                                delete $SNP_active_back{$id_split2};
                                delete $before_back{$id_split2};
                                delete $before_shorter_skip_back{$id_split2};
                                $id = $id_split2;
                                if ($correction > 0)
                                {
                                    delete $seed{$id_split2};
                                    $best_extension = $best_extension3;
                                    goto AFTER_EXT_BACK;
                                }
                                else
                                {
                                    $read_new = $seed{$id_split2};                       
                                    goto FINISH;
                                }
                            }
                            elsif (($count4b > 2 && ($count1b+$count2b+$count3b) eq '0') || $differenceb4 eq "yes")
                            {
                                if (@extensions_before > 15 && $repetitive_detect_back2 eq "")
                                {
                                    $ext_before = "yes";
                                    @extensions_before = @extensions_before4;
                                    goto EXT_BEFORE_BACK;
                                }
                                delete $seed{$id};
                                delete $seed{$id_split1};
                                delete $seed{$id_split2};
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nREVERSE_DELETE SPLIT_BEFORE8\n\n";
                                }
                                delete $SNP_active_back{$id_split3};
                                delete $before_back{$id_split3};
                                delete $before_shorter_skip_back{$id_split3};
                                $id = $id_split3;
                                if ($correction > 0)
                                {
                                    delete $seed{$id_split3};
                                    $best_extension = $best_extension4;
                                    goto AFTER_EXT_BACK;
                                }
                                else
                                {
                                    $read_new = $seed{$id_split3};                       
                                    goto FINISH;
                                }
                            }
                            elsif(($count1b+$count2b+$count3b+$count4b) < 10 && $overhang < $read_length-($overlap*1.5) && $before_shorter_skip_back ne "yes" && $skip_overhang ne "yes")
                            {
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nMAKE BEFORE SHORTER3\n\n";
                                }
                                $before_shorter = "yes";
                                goto BEFORE_BACK;
                            }
                            elsif(($count1b+$count2b+$count3b+$count4b) <= 10 && $overhang < $read_length-($overlap*1.5) && $before_shorter_skip_back ne "yes" && $insert_range_before eq "yes")
                            {
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nMAKE BEFORE SHORTER3\n\n";
                                }
                                $insert_range_shorter = "";
                                goto BEFORE_BACK;
                            }
                            elsif(($count1b+$count2b+$count3b+$count4b) > 10 && $before_shorter_skip_back ne "yes" && $insert_range_shorter eq "")
                            {
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nSMALLER INSERT_BACK\n\n";
                                }
                                $insert_range_shorter = "yes";
                                goto BEFORE_BACK;
                            }
                            elsif (($count1b+$count2b+$count3b+$count4b) > 9 && $extensions_before ne "yes" && $repetitive_detect_back2 ne "yes" && $before_shorter_skip_back ne "yes")
                            {
                                $l = 0;
                                $best_extension = "";
                                $SNP = "";
                                @extensions = @extensions_before;
                                $ext = '0';
                                foreach (@extensions)
                                {
                                    $ext++;
                                }
                                $extensions_before = "yes";
                                if ($y > $startprint2)
                                {
                                    print OUTPUT5 "\nEXTENSIONS FROM BEFORE BACK\n\n";
                                }
                                $SNP_active_back = "yes";
                                delete $before_back{$id};
                                $before_shorter_skip_back{$id} = "yes";
                                delete $seed{$id};
                                delete $seed{$id_split1};
                                delete $seed{$id_split2};
                                delete $seed{$id_split3};
                                undef @extensions_group1;
                                undef @extensions_group2;
                                undef @extensions_group3;
                                undef @extensions_group4;
                                goto NUCLEO_BACK;
                            }         
                        }                                
EXT_BEFORE_BACK:        if ($ext_before eq "yes")
                        {
                            $l = 0;
                            $best_extension = "";
                            $SNP = "";
                            @extensions = @extensions_before;
                            $ext = '0';
                            foreach (@extensions)
                            {
                                $ext++;
                            }
                            $extensions_before = "yes";
                            if ($y > $startprint2)
                            {
                                print OUTPUT5 "\nEXTENSIONS FROM BEFORE BACK2\n\n";
                            }
                            $SNP_active_back = "yes";
                            delete $before_back{$id};
                            $before_shorter_skip_back{$id} = "yes";
                            delete $seed{$id};
                            delete $seed{$id_split1};
                            delete $seed{$id_split2};
                            delete $seed{$id_split3};
                            undef @extensions_group1;
                            undef @extensions_group2;
                            undef @extensions_group3;
                            undef @extensions_group4;
                            goto NUCLEO_BACK;
                        }    
                                delete $seed{$id_split1};
                                delete $seed{$id_split2};
                                delete $seed{$id_split3};
                                $before_back{$id} = "yes";
                                $best_extension = "";
                                $read_new = $read_new1;
                                $before_shorter_skip_back{$id} = "yes";
                                if ($before_shorter_skip_back eq "yes")
                                {
                                    $before_back{$id} = "yes";
                                }
                                
                                if ($count_split > 2)
                                {
                                    my $tmp = '0';
                                    if ($first_yuyu ne "yes")
                                    {
                                        $yuyu_option_back{$id.$firstSNP_max[$tmp]} = "yes";
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $firstSNP_max[$tmp]." FIRST_YUYU\n";
                                        }
                                    }
                                    $tmp++;
                                    if ($second_yuyu ne "yes")
                                    {
                                        $yuyu_option_back{$id.$firstSNP_max[$tmp]} = "yes";
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $firstSNP_max[$tmp]." SECOND_YUYU\n";
                                        }
                                    }
                                    $tmp++;
                                    if ($third_yuyu ne "yes")
                                    {
                                        $yuyu_option_back{$id.$firstSNP_max[$tmp]} = "yes";
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $firstSNP_max[$tmp]." THIRD_YUYU\n";
                                        }
                                    }
                                    $tmp++;
                                    if ($fourth_yuyu ne "yes")
                                    {
                                        $yuyu_option_back{$id.$firstSNP_max[$tmp]} = "yes";
                                        if ($y > $startprint2)
                                        {
                                            print OUTPUT5 $firstSNP_max[$tmp]." FOURTH_YUYU\n";
                                        }
                                    }
                                }
                                goto FINISH;
                    }
                }
                print OUTPUT5 $delete_first." DLETE_FIRST\n";
                if (($extensions_before eq "yes" || ($indel_split_skip_back ne "yes"  && ($delete_first ne "yes_back" || (length($best_extension1) < 5 && length($best_extension2) < 5)) && $ext_total_back > 18 && $indel_split_back < ($read_length-25 -$overlap) && $SNR_read_back eq "")) && $deletion_back ne "yes")
                {
                    $best_extension = "";
                    $indel_split_back{$id} = $indel_split_back+10;
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "INCREASE_INDEL_SPLIT_BACK\n";
                    }
                    goto FINISH;
                }
                elsif ($SNP_active_back eq "" && ($delete_first ne "yes_back" || (length($best_extension1) < 5 && length($best_extension2) < 5)) && $deletion_back ne "yes")
                {
                    $best_extension = "";
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "SNP_ACTIVE_BACK\n";
                    }
                    delete $indel_split_back{$id};
                    $SNP_active_back{$id} = "yes";
                    goto FINISH;
                }
                elsif ($deletion_back eq "yes")
                {
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "DELETION_DETECT\n";
                    }
                    delete $seed{$id_split1};
                    delete $seed{$id};
                    delete $before_back{$id};
                    delete $before_shorter_skip_back{$id};
                    $read_new = $read_new1;     
                }
                else
                {
                    if ($y > $startprint2)
                    {
                        print OUTPUT5 "ELSE\n";
                    }
                    delete $seed{$id_split1};
                    delete $seed{$id};
                    $indel_split_back = '0';
                    delete $indel_split_back{$id};
                    $noback{$id} = "yes";
                    $read_new = $read_new1;
                }
            }
            else
            {                             
                if ($best_extension ne "")
                {
                    $indel_split_back = '0';
                    delete $indel_split_back{$id};
                }    

                my @ext = split //, $best_extension;
                my $u = length($best_extension);
                my $v = '1';
 
                    if (($SNR_test eq "yes2_back22222" || $SNR_test eq "yes2_double_back") && $u > 2)
                    {
                        my @ext3 = split //, $best_extension;
                        $v = '1';
                        $u = length($best_extension);
                        while ($ext3[$u-$v-1] eq $ext3[$u-$v] && ($u-$v-1) > 1)
                        {              
                            chop($best_extension);
                            $v++;
                        }
                        my $SNR_max = '0';
                        my $SNR_min = '1000';
                        my $n = '0';
                        my $extensions_after_SNR = substr $best_extension, 0,4;
                        $extensions_after_SNR =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;

                        foreach my $SNR_ext (keys %SNR_count_back)
                        {
                            my $SNR_count_back = $SNR_count_back{$SNR_ext};
                            my $checkie = substr $SNR_ext, 0, 4;
                            if ($SNR_count_back ne "" && $checkie eq $extensions_after_SNR)
                            {
                                if ($SNR_count_back > $SNR_max)
                                {
                                    $SNR_max = $SNR_count_back;
                                }
                                if ($SNR_count_back < $SNR_min)
                                {
                                    $SNR_min = $SNR_count_back;
                                    print $SNR_ext." MIN\n";
                                }
                            }    
                        }
                        if ($SNR_min eq '1000')
                        {
                            $SNR_min = '0';
                        }
                        
                        my $p = '0';
                        my $ut = '0';
                        if ($SNR_test eq "yes2_back" && $SNR_max > 0)
                        {
                            $best_extension = $SNR_nucleo.$best_extension;
                            $ut = '1';
                            while ($p < $SNR_max - $SNR_min)
                            {
                                $best_extension = "X".$best_extension;
                                $p++;
                            }
                        }
                        elsif ($SNR_test eq "yes2_double_back" && $SNR_max > 0)
                        {
                            $best_extension = $SNR_nucleo.$best_extension;
                            $ut = '2';
                            while ($p < ($SNR_max - $SNR_min)/2)
                            {
                                $best_extension = "X2".$best_extension;
                                $p++;
                            }
                        }
                        $p = '0';
                        
                        while ($p+$ut < $SNR_min)
                        {
                            $best_extension = $SNR_nucleo.$best_extension;
                            if ($SNR_test eq "yes2_back")
                            {
                                $p++; 
                            }
                            if ($SNR_test eq "yes2_double_back")
                            {
                                $p++;
                                $p++; 
                            }   
                        }
                        delete $SNR_back{$id};
                    }
                    elsif ($SNR_test eq "yes2_back2222" || $SNR_test eq "yes2_double_back")
                    {
                        my $ee2 = '10000'; 
                        foreach (keys %SNR_count_back)
                        {
                            my $ee = $SNR_count_back{$_}; 
                            if ($ee < $ee2)
                            {
                                $ee2 = $ee;
                            }
                        }
                        $best_extension = "";
                        my $ee3 = '0';
                        if ($ee2 ne '10000')
                        {
                            while ($ee3 < $ee2)
                            {
                                $best_extension = $SNR_nucleo.$best_extension;
                                $ee3++;
                            }   
                        }
                        else
                        {
                            $best_extension = "";
                        }
                    }
                
                $SNR_test = "";
                if ($y > $startprint2)
                {
                    print OUTPUT5 reverse($best_extension)." BEST_EXTENSION_BACK\n\n";
                }
                $best_extension_back_prev{$id} = $best_extension;
                
               
            }  
AFTER_EXT_BACK:
            undef @matches;
            undef @matches1;
            undef @matches2;
                                            chomp $best_extension;
                                            
                                            my $vk2 = '0';
                                            
                                            if ($SNR_read_back eq "")
                                            {
                                                my @dot2 = split //, $best_extension;
                                                my $ut2 = length($best_extension);
                                                
                                                while ($dot2[$ut2-1] eq "." || $dot2[$ut2-1] eq "*")
                                                {              
                                                    if ($dot2[$ut2-1] eq "*")
                                                    {
                                                        chop $best_extension;
                                                        chop $best_extension;
                                                        $vk2++;
                                                        $vk2++;
                                                        $ut2--;
                                                        $ut2--;
                                                    }
                                                    else
                                                    {                                                     
                                                        chop $best_extension;
                                                        $vk2++;
                                                        $ut2--;
                                                    }
                                                }
                                                
                                                my $SNR_check2 = $best_extension =~ qr/AAAAAAA|CCCCCCC|GGGGGGG|TTTTTTT/;
                                                if ($SNR_check2 > 0)
                                                {
                                                    if ($best_extension =~ m/(.*?(AAAAAA|CCCCCC|GGGGGG|TTTTTT)).*/)
                                                    {
                                                        $best_extension = $1;
                                                        if ($y > $startprint2)
                                                        {
                                                            print OUTPUT5 reverse($best_extension)." BEST_EXTENSION_SHORT_BACK\n\n";
                                                        }
                                                    }
                                                }
                                            }
 
                                            if ($noback ne "stop" && ($best_extension ne "" && $best_extension ne " "))
                                            {                                                                               
                                                $read_new = reverse($best_extension).$read_new1;                                     
                                                $position_back += length($best_extension);                                        
                                                $position_back -= $vk2;
                                                                               
                                                my $position_back_tmp = $position_back{$id};
                                                delete $position_back{$id};
                                                $position_back{$id} = $position_back;
                                                $position_back = $position_back_tmp;
                                                if ($best_extension ne "" && $split_forward eq "")
                                                {
                                                   delete $indel_split_skip_back{$id}; 
                                                }
                                                
                                                $best_extension = "";
                                                delete $regex_back{$id};
                                                $use_regex_back = "";
                                                delete $last_chance_back{$id};
                                                $last_chance_back = "";
                                                if ($split_forward eq "")
                                                {
                                                    delete $before_back{$id};
                                                    delete $before_shorter_skip_back{$id};
                                                }   
                                                
                                                $id_test = $id;
                                                
                                                if ($split eq "")
                                                {
                                                    delete $yuyu_option_back{$id.'A'};
                                                    delete $yuyu_option_back{$id.'C'};
                                                    delete $yuyu_option_back{$id.'T'};
                                                    delete $yuyu_option_back{$id.'G'};
                                                }
                                                if ($split eq "")
                                                {
                                                    foreach my $add_read2 (keys %extensions)
                                                    {
                                                        my $add_read = substr $extensions{$add_read2}, 0, -1;
                                                        $count_reads{$add_read} = undef;
                                                    }
                                                }
                                            }
                                            elsif ($use_regex_back ne "yes" && $noback ne "stop" && $indel_split_back > 0)
                                            {                                              
                                                delete $regex_back{$id};
                                                if ($split_forward eq "")
                                                {
                                                    $indel_split_skip_back{$id} = "yes";
                                                    $indel_split_back = '0';
                                                    delete $indel_split_back{$id};
                                                }
                                                elsif ($split_forward ne "" && exists($SNP_active_back{$id}))
                                                {
                                                    $indel_split_skip_back{$id} = "yes";
                                                    $indel_split_back = '0';
                                                    delete $indel_split_back{$id};
                                                }
                                                elsif ($split_forward ne "")
                                                {
                                                    $SNP_active_back{$id} = "yes";
                                                }
                                                $read_new = $read_new1;
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "2B\n";
                                                }
                                            }
                                            elsif ($use_regex_back ne "yes" && $noback ne "stop" && $repetitive_detect_back eq "")
                                            {                                              
                                                delete $regex_back{$id};
                                                $regex_back{$id} = "yes";
                                                $read_new = $read_new1;
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "3B\n";
                                                }
                                            }
                                            elsif ($use_regex_back ne "yes" && $last_chance_back ne "yes" && $indel_split_back > 0 && $noback ne "stop")
                                            {
                                                $read_new = $read_new1;
                                                $indel_split_skip_back{$id} = "yes";
                                                $indel_split_back = '0';
                                                delete $indel_split_back{$id};
                                                $SNP_active_back{$id} = "yes";
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "4B\n";
                                                }
                                            }
                                            elsif (($use_regex_back eq "yes" || $repetitive_detect_back ne "") && $last_chance_back ne "yes" && $noback ne "stop" && $AT_rich_before_back eq "")
                                            {
                                                $read_new = $read_new1;
                                                delete $last_chance_back{$id};
                                                $last_chance_back{$id} = "yes";
                                                delete $regex_back{$id};
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "5B\n";
                                                }
                                            }
                                            elsif ($last_chance_back eq "yesssssssssssssssssssssssssssssssssss" && $SNR_read_back ne "" && $noback ne "stop")
                                            {
                                                $read_new = $read_new1;
                                                delete $last_chance_back{$id};
                                                $last_chance_back{$id} = "yes";
                                                $last_chance = "";
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "6B\n";
                                                }
                                            }
                                            elsif ($last_chance_back eq "yes" && $noback ne "stop" && $use_regex_back ne "yes")
                                            {
                                                $read_new = $read_new1;
                                                delete $last_chance_back{$id};
                                                $last_chance_back{$id} = "yes";
                                                delete $regex_back{$id};
                                                $regex_back{$id} = "yes";
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "7B\n";
                                                }
                                            }
                                            elsif (($last_chance_back eq "yes" || $AT_rich_before_back ne "") && $noforward ne "stop" && $noback ne "stop" && $use_regex_back eq "yes")
                                            {
                                                $read_new = $read_new1;
                                                delete $last_chance_back{$id};
                                                $noback = "stop";
                                                $noback{$id} = "stop";
                                                delete $regex_back{$id};
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "8B\n";
                                                }
                                            }
                                            elsif (exists($old_id{$id}) && $noforward eq "stop")
                                            {
                                                $merge_now = "yes2";
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "8B_MERGE\n";
                                                }
                                                goto MERGE;
                                            }
                                            elsif (exists($seed_split{$id}))
                                            {
                                                delete $seed{$id};
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "9B\n";
                                                }
                                                my $old_rep_old_check = "";
                                                if (exists ($old_rep_old{$id}))
                                                {
                                                    $old_rep_old_check = "yes";
                                                }
                                                if (length($read) > 250 && $old_rep_old_check eq "")
                                                {
                                                    if ($id =~ m/(.*_|)(\d+)$/)
                                                    {
                                                        my $check_id = $2;
                                                        if ($contig_num eq '1')
                                                        {       
                                                            $contigs{$contig_num."+".$check_id} = $read;
                                                            my $start_point = '500';
                                                            my $check_repetitive = '3';
                                                                                                
                                                            while ($check_repetitive > 2)
                                                            {
                                                                my $check_for_rep = substr $read, 0, $start_point+200;
                                                                my $repetitive = substr $check_for_rep, $start_point, 15;
                                                                $check_repetitive = $check_for_rep =~ s/$repetitive/$repetitive/g;
                                                                {
                                                                    $start_point += 20;
                                                                    if ($y > $startprint2)
                                                                    {
                                                                        print OUTPUT5 "DETECT_REPETITIVE_IN_START_SEQUENCE\n";
                                                                    }
                                                                } 
                                                            }
                                                            $first_contig_start = substr $read, $start_point, $overlap;
                                                            my $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                                            while ($check_start > 0)
                                                            {
                                                                $start_point += 10;
                                                                $first_contig_start = substr $read, $start_point, $overlap;
                                                                $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                                            }
                                                            if ($y > $startprint2)
                                                            {
                                                                print OUTPUT5 $first_contig_start." CONTIG_START5\n";
                                                            }
                                                            foreach my $seedie (keys %seed)
                                                            {
                                                                my $seedie_part = substr $seed{$seedie}, 0, 1000;
                                                                my $check_seedie = $seedie_part =~ s/$first_contig_start/$first_contig_start/;
                                                                if ($check_seedie > 0)
                                                                {
                                                                    delete $seed{$seedie};
                                                                }
                                                            }
                                                        }
                                                        else
                                                        {
                                                           $contigs{$contig_num."+".$check_id} = $read; 
                                                        }
                                                        $contig_num++;
                                                        $seed_old{$id} = $read;
                                                    }
                                                }
SKIP1:                                          if ($y > $startprint2)
                                                {                                                                                 
                                                    print OUTPUT5 ">".$id."\n";
                                                    print OUTPUT5 $read."\n";
                                                }
                                                goto SEED;
                                            }
                                            else
                                            {
                                                $noback = "stop";
                                                $noback{$id} = "stop";
                                                print OUTPUT5 "10B\n";
                                            }
                                            if ($split eq "yes2_back")
                                            {                      
                                                $seed{$id} = $read_new;
                                                $insert_size2{$id} = $insert_size2;
                                                $read_new = $read_new1;
                                                
                                                if ($count_split eq '2')
                                                {
                                                    $split = "yes4_back";
                                                }
                                                goto SPLIT_BACK;
                                            }
                                            elsif ($split eq "yes3_back")
                                            {                      
                                                $seed{$id} = $read_new;
                                                $insert_size2{$id} = $insert_size2;
                                                $read_new = $read_new1;

                                                if ($count_split eq '3')
                                                {
                                                    $split = "yes4_back";
                                                }                                               
                                                goto SPLIT_BACK;
                                            }
                                            if ($split eq "yes4_back")
                                            {                      
                                                $seed{$id} = $read_new;
                                                $insert_size2{$id} = $insert_size2;
                                                $read_new = $read_new1;
                                                
                                                goto SPLIT_BACK;
                                            }
        }
FINISH:
                                            my $count_seed = '0';
                                            foreach my $count_seed2 (keys %seed)
                                            {
                                                $count_seed++;
                                            }
                                           if ($y > $startprint2)
                                            {
                                                if ($noback ne "")
                                                {
                                                    print OUTPUT5 $noback." NOBACK\n";
                                                }
                                                if ($noforward ne "")
                                                {
                                                    print OUTPUT5 $noforward." NOFORWARD\n";
                                                }
                                                if ($last_chance ne "")
                                                {
                                                    print OUTPUT5 $last_chance." LASTCHANCE\n";
                                                }
                                                if ($use_regex_back ne "")
                                                {
                                                    print OUTPUT5 $use_regex_back." REGEX_BACK\n";
                                                }
                                                if ($last_chance_back ne "")
                                                {
                                                    print OUTPUT5 $last_chance_back." LASTCHANCE_BACK\n";
                                                }
                                                print OUTPUT5 $count_seed." COUNTSEED\n";
                                            }
                                            if ($delete_first eq "yes2" && $indel_split > 0 && $SNR_next_seed ne "yes")
                                            {
                                                $indel_split = '0';
                                                $indel_split_skip = "yes";
                                                delete $indel_split{$id};
                                                $noforward{$id} = "";
                                                $seed{$id} = $read_new;
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "FINISH1\n";
                                                }
                                            } 
                                            elsif ($merge ne "yes" && $circle eq "" && ($noback ne "stop" || ($noforward ne "stop" && $delete_first ne "yes2")) && $AT_rich ne "yes" && $bad_read ne "yes" && $SNR_next_seed ne "yes")
                                            {
                                                delete $seed{$id};                                         
                                                $seed{$id} = $read_new;
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "FINISH2\n";
                                                }
                                            }
                                            elsif ($merge ne "yes" && $circle eq "" && $AT_rich ne "yes" && $bad_read ne "yes" && $delete_first eq "yes2" && $last_chance ne "yes" && $SNR_next_seed ne "yes")
                                            {
                                                delete $seed{$id};                                         
                                                $seed{$id} = $read_new;
                                                delete $last_chance{$id};
                                                $last_chance{$id} = "yes";
                                                if ($y > $startprint2)
                                                {
                                                    print OUTPUT5 "FINISH3\n";
                                                }
                                            }

                                            
                                            elsif ($no_next_seed ne "yes" && ($SNR_next_seed eq "yes" || ($nosecond eq "" && $CP_check ne "yes" && length($read) > $read_length+150 && $circle eq "" && (($last_chance eq "yes" || $noforward eq "stop") && ($last_chance_back eq "yes" || $noback eq "stop")) || ($AT_rich eq "yes" && $count_seed ne "0") || ($bad_read eq "yes" && $count_seed ne "0"))))
                                            {
                                                if ($id_original eq "")
                                                {
                                                    $id_original = $id;
                                                }
                                                if ($bad_read eq "yes")
                                                {
                                                    $read = $seed_old{$id_old};
                                                }
                                                else
                                                {
                                                    if ($y > $startprint2)
                                                    {
                                                        print OUTPUT5 ">".$id."\n";
                                                        print OUTPUT5 $read."\n";
                                                    }
                                                    my $lastbit_contig = substr $read, -20;
                                                    $lastbit_contig =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    
                                                    if (length($read) > $read_length+70 && ($lastbit_contig_prev !~ m/.*.$lastbit_contig.*/ || length($read) > $read_length+300))
                                                    {
                                                        print OUTPUT6 ">".$id."\n";
                                                        print OUTPUT6 $read."\n";
                                                        $lastbit_contig_prev = substr $read, -100;
                                                        $lastbit_contig_prev =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                                                    }                                                 

                                                    if (exists $old_id{$id})
                                                    {
                                                        if (exists($contigs_id{$old_id{$id}}) && length($seed_old{$old_id{$id}}) > $read_length+30)
                                                        {
                                                        }
                                                        else
                                                        {  
                                                            print OUTPUT5 $seed_old{$old_id{$id}}." SEED_OLD_CONTIG\n";
                                                            if ($contig_num eq '1')
                                                            {
                                                                $contigs{$contig_num."+".$id} = $seed_old{$old_id{$id}};
                                                                my $start_point = '500';
                                                                my $check_repetitive = '3';
                                                                                            
                                                                while ($check_repetitive > 2)
                                                                {
                                                                    my $check_for_rep = substr $seed_old{$old_id{$id}}, 0, $start_point+200;
                                                                    my $repetitive = substr $check_for_rep, $start_point, 15;
                                                                    $check_repetitive = $check_for_rep =~ s/$repetitive/$repetitive/g;
                                                                    if ($check_repetitive > 2)
                                                                    {
                                                                        $start_point += 20;
                                                                        print OUTPUT5 "DETECT_REPETITIVE_IN_START_SEQUENCE\n";
                                                                    } 
                                                                }
                                                                $first_contig_start = substr $seed_old{$old_id{$id}}, $start_point, $overlap;
                                                                my $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                                                while ($check_start > 0)
                                                                {
                                                                    $start_point += 10;
                                                                    $first_contig_start = substr $seed_old{$old_id{$id}}, $start_point, $overlap;
                                                                    $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                                                }
                                                                print OUTPUT5 $first_contig_start." CONTIG_START6\n";
                                                                
                                                                foreach my $seedie (keys %seed)
                                                                {
                                                                    my $seedie_part = substr $seed{$seedie}, 0, 1000;
                                                                    my $check_seedie = $seedie_part =~ s/$first_contig_start/$first_contig_start/;
                                                                    if ($check_seedie > 0)
                                                                    {
                                                                        delete $seed{$seedie};
                                                                    }
                                                                }
                                                            }
                                                            else
                                                            {
                                                                $contigs{$contig_num."+".$id} = $seed_old{$old_id{$id}};
                                                            }
                                                            $contig_num++;
                                                        }
                                                        delete $old_id{$id};
                                                    }
                                                    $seed_old{$id} = $read;
                                                    $id_old = $id;
                                                }                                     
                                                
                                                my $SNR_skip = "yes";
                                                my $xy = -($overlap+3);
                                                my $tt = '-250';
                                                my $second_seed = "";                                           
                                                my $u = '10';
                                                while ($SNR_skip eq "yes")
                                                { 
                                                    my $new_seed_part_tmp = substr $read, -$u, 10;
                                                    my $SNR_checkA = $new_seed_part_tmp =~ tr/A/A/;
                                                    my $SNR_checkC = $new_seed_part_tmp =~ tr/C/C/;
                                                    my $SNR_checkT = $new_seed_part_tmp =~ tr/T/T/;
                                                    my $SNR_checkG = $new_seed_part_tmp =~ tr/G/G/;
                                                    if ($SNR_checkA > 6 || $SNR_checkC > 6 || $SNR_checkG > 6 || $SNR_checkT > 6)
                                                    {
                                                        $SNR_skip = "yes";
                                                        $u += 10;
                                                    }
                                                    else
                                                    {
                                                        $SNR_skip = "";
                                                    }
                                                }
                                                
NEXT_SEED:                                      while ($xy > $tt)
                                                {
                                                    my $new_seed_part = substr $read, $xy-$u+10, $overlap;
                                                    
                                                    if (exists($hash2c{$new_seed_part}))
                                                    {
                                                        my $id_tmp = $hash2c{$new_seed_part};
                                                        my $id_tmp2 = substr $id_tmp, 1;
                                                        my @id_tmp = split /,/,$id_tmp2;
                                                            
                                                        foreach my $id_tmp3 (@id_tmp)
                                                        {
                                                            chomp ($id_tmp3);
                                                            my $id_match_end = substr $id_tmp3, -1, 1,"",;
                                                            my $seed_tmp;

                                                            if (exists($hash{$id_tmp3}) && $id_tmp3)
                                                            {
                                                                my @id_tmp3 = split /,/, $hash{$id_tmp3};
                                                                if ($id_match_end eq "1")
                                                                {
                                                                    $seed_tmp = $id_tmp3[1];
                                                                }
                                                                elsif ($id_match_end eq "2")
                                                                {
                                                                    $seed_tmp  = $id_tmp3[0];
                                                                }
                                                                if ($encrypt eq "yes")
                                                                {
                                                                    $seed_tmp = decrypt $seed_tmp;
                                                                }
                                                                
                                                                $seed_tmp =~ tr/ACTG/TGAC/;
                                                                my $seed_tmp2 = reverse($seed_tmp);                             
                                                                my $yx = '0';
                                                                while ($yx < length($seed_tmp2)-($overlap+1))
                                                                {
                                                                    my $seed_tmp2_part = substr $seed_tmp2, $yx, $overlap;
                                                                    if (exists($hash2b{$seed_tmp2_part}))
                                                                    {
                                                                        my $id2 = $hash2b{$seed_tmp2_part};
                                                                        my $id5 = substr $id2, 1;
                                                                        my @id_tmp2 = split /,/,$id5;
                                                                        $id = $id_tmp2[0];
                                                                        my $id_b = $id;
                                                                        my $id_tmp2_end = substr $id_b, -1, 1,"",;
                                                                        
                                                                        if ($y > $startprint2)
                                                                        {
                                                                            print OUTPUT5 $id_b." ID_B!!!!\n";
                                                                        }
                                                                    
                                                                        if ($second_seed eq "ddd")
                                                                        {
                                                                            my @id_b = split /,/, $hash{$id_b};
                                                                            if ($id_tmp2_end eq "1")
                                                                            {
                                                                                $read = $id_b[0];
                                                                            }
                                                                            elsif ($id_tmp2_end eq "2")
                                                                            {
                                                                                $read  = $id_b[1];
                                                                            }
                                                                            $second_seed = "yes";
                                                                            if ($encrypt eq "yes")
                                                                            {
                                                                                $read = decrypt $read;
                                                                            }
                                                                            
                                                                            $tt = -length($read);
                                                                            goto NEXT_SEED;
                                                                        }
                                                                        if (exists($id_bad{$id_b}))
                                                                        {
                                                                            goto SAME_ID;
                                                                        }  
                                                                        elsif (exists($hash{$id_b}))
                                                                        {
                                                                            my @id_b = split /,/, $hash{$id_b};
                                                                            if ($id_tmp2_end eq "1")
                                                                            {
                                                                                $seed = $id_b[0];
                                                                            }
                                                                            elsif ($id_tmp2_end eq "2")
                                                                            {
                                                                                $seed  = $id_b[1];
                                                                            }
                                                                            if ($encrypt eq "yes")
                                                                            {
                                                                                $seed = decrypt $seed;
                                                                            }

                                                                            my $part2 = substr $seed, -$overlap-10;
                                                                            my $s = '0';
                                                                            my $most_match_total = '0';
                                                                            while ($s < length($part2)-$overlap)
                                                                            {
                                                                                my $part2b = substr $part2, $s, $overlap;
                                                                                my $part2b_reverse = reverse($part2b);
                                                                                $part2b_reverse =~ tr/ACTG/TGAC/;
                                                                                if (exists($hash2b{$part2b}))
                                                                                {
                                                                                    $most_match_total++;
                                                                                }
                                                                                if (exists($hash2c{$part2b}))
                                                                                {
                                                                                    $most_match_total++;
                                                                                }
                                                                                if (exists($hash2b{$part2b_reverse}))
                                                                                {
                                                                                    $most_match_total++;
                                                                                }
                                                                                if (exists($hash2c{$part2b_reverse}))
                                                                                {
                                                                                    $most_match_total++;
                                                                                }
                                                                                $s++;
                                                                            }
                                                                            if ($most_match_total > 5)
                                                                            {
                                                                            }
                                                                            else
                                                                            {
                                                                               $xy--;
                                                                               goto NEXT_SEED;
                                                                            }
                                                                        }
                                                                        else
                                                                        {
                                                                           goto SAME_ID;
                                                                        }
                                                                        if ($y > $startprint2)
                                                                        {
                                                                            print OUTPUT5 $id." SECOND!!!!\n";
                                                                            print OUTPUT5 $seed." SEED!!!!\n";
                                                                            print OUTPUT5 $id_original." ID_O!!!!\n";
                                                                        }
                                                                        delete $seed{$id_original};
                                                                        delete $seed{$id};
                                                                        delete $seed{$id_split1};
                                                                        delete $seed{$id_split2};
                                                                        delete $seed{$id_split3};
                                                                        $seed = correct ($seed);
   
                                                                        $id_bad{$id_b} = undef;
                                                                                                                                              
                                                                        $seed{$id} = $seed;
                                                                        $seeds_check{$id} = undef;
                                                                        $position{$id} = length($seed);
                                                                        $old_id{$id} = $id_old;
                                                                        $old_id2{$id} = undef;
                  
                                                                        delete $position_adjust{$id};
                                                                        $insert_size2{$id} = $insert_size;
                                                                        
                                                                        my $count_contig_tmp = $contig_count;
                                                                        while ($count_contig_tmp > 0)
                                                                        {
                                                                            $contig_gap_min{$id."_".$count_contig_tmp} = $contig_gap_min{$id_old."_".$count_contig_tmp};
                                                                            $contig_gap_max{$id."_".$count_contig_tmp} = $contig_gap_max{$id_old."_".$count_contig_tmp};
                                                                            $count_contig_tmp--;    
                                                                        } 
                                                                        
                                                                        $contig_count++;
                                                                        $contig_count{$id} = $contig_count;
                                                                        my $gap_min = ($insert_size*0.62)-(2*$read_length)+$xy + $overlap-16 + $yx;
                                                                        my $gap_max= ($insert_size*1.38)-(2*$read_length)+$xy + $overlap-16 + $yx; 
                                                                        $contig_gap_min{$id."_".$contig_count} = $gap_min;
                                                                        $contig_gap_max{$id."_".$contig_count} = $gap_max;
                                                                        goto ITERATION;
SAME_ID:       
                                                                    }
                                                                    $yx++;
                                                                }
                                                            }
                                                        }
                                                    }
                                                    $xy--;
                                                }
                                                $no_next_seed = "yes";
                                                goto FINISH;
                                            }
                                            elsif($circle eq "" && (($last_chance eq "yes" || $noforward eq "stop") && ($last_chance_back eq "yes" || $noback eq "stop")) || ($AT_rich eq "yes" && $count_seed ne "0") || ($bad_read eq "yes" && $count_seed ne "0"))
                                            {
                                                delete $seed{$id};                                         
                                                delete $last_chance{$id};
                                                if ($y > $startprint2)
                                                {                                                                                 
                                                    print OUTPUT5 "DELETE READS AND SEED\n";
                                                    print OUTPUT5 ">".$id."\n";
                                                    print OUTPUT5 $read."\n";
                                                }
                                                if (exists($old_id{$id}))
                                                {
                                                    my $read_tmp = $seed_old{$old_id{$id}};
                                                    if (length($read_tmp) > 250)
                                                    {
                                                        $contigs{$contig_num."+".$old_id{$id}} = $read_tmp;
                                                        $contig_num++;
                                                    }
                                                }                                               
                                                else
                                                {
                                                    if (length($read) > 250)
                                                    {
                                                        if ($id =~ m/(.*_|)(\d+)$/)
                                                        {
                                                            my $check_id = $2;
                                                            if ($contig_num eq '1')
                                                            {       
                                                                $contigs{$contig_num."+".$check_id} = $read;
                                                                my $start_point = '500';
                                                                my $check_repetitive = '3';
                                                                                                    
                                                                while ($check_repetitive > 2)
                                                                {
                                                                    my $check_for_rep = substr $read, 0, $start_point+200;
                                                                    my $repetitive = substr $check_for_rep, $start_point, 15;
                                                                    $check_repetitive = $check_for_rep =~ s/$repetitive/$repetitive/g;
                                                                    {
                                                                        $start_point += 20;
                                                                        print OUTPUT5 "DETECT_REPETITIVE_IN_START_SEQUENCE\n";
                                                                    } 
                                                                }
                                                                $first_contig_start = substr $read, $start_point, $overlap;
                                                                my $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                                                while ($check_start > 0)
                                                                {
                                                                    $start_point += 10;
                                                                    $first_contig_start = substr $read, $start_point, $overlap;
                                                                    $check_start = $first_contig_start =~ tr/N|K|R|Y|S|W|M|B|D|H|V|\.//;
                                                                }
                                                                print OUTPUT5 $first_contig_start." CONTIG_START7\n";
                                                                foreach my $seedie (keys %seed)
                                                                {
                                                                    my $seedie_part = substr $seed{$seedie}, 0, 1000;
                                                                    my $check_seedie = $seedie_part =~ s/$first_contig_start/$first_contig_start/;
                                                                    if ($check_seedie > 0)
                                                                    {
                                                                        delete $seed{$seedie};
                                                                    }
                                                                }
                                                            }
                                                            else
                                                            {
                                                               $contigs{$contig_num."+".$check_id} = $read; 
                                                            }
                                                            $contig_num++;
                                                        }
                                                    }
                                                }
                                            }
                                            elsif($circle eq "yes")
                                            {
                                                my $output_file  = "Circularized_assembly_".$option."_".$project.".fasta";
                                                open(OUTPUT, ">" .$output_file) or die "Can't open file $output_file, $!\n";
                                                
                                                delete $seed{$id};
                                                $read =~ tr/\./N/;
                                                $read =~ tr/X//d;
                                                my @contigs = split /L+/, $read;
                                                my $l = '0';
                                                my $largest_contig = '0';
                                                my $miminum_contig = '100000000000000000000000000000';
                                                
                                                print "\b" x length($progress_before);
                                                print ' ' x length($progress_before);
                                                print "\b" x length($progress_before);
                                                print "\n-----------------Assembly ".$option." finished successfully: The genome has been circularized-----------------\n\n";
                                                print OUTPUT4 "\n-----------------Assembly ".$option." finished successfully: The genome has been circularized-----------------\n\n";
                                                print OUTPUT5 "\n-----------------Assembly ".$option." finished successfully: The genome has been circularized-----------------\n\n";
                                                $assembly_length = '1';            
                                                foreach (@contigs)
                                                {
                                                    $l++;
                                                    my $fin = $_;
                                                    my $fin2 = $fin;
                                                    $fin =~ s/(.{1,150})/$1\n/gs;
                                                    
                                                    if ($l > 1)
                                                    {
                                                        my $gap_min = sprintf("%.0f", $contig_gap_min{$id."_".($l-1)});
                                                        my $gap_max = sprintf("%.0f", $contig_gap_max{$id."_".($l-1)});
                                                        print "Estimated gap      : ".$gap_min." bp to ".$gap_max." bp";
                                                        print OUTPUT4 "Estimated gap      : ".$gap_min." bp to ".$gap_max." bp";
                                                        if ($contig_gap_min{$id."_".($l-1)} < 0)
                                                        {
                                                            print " (Check manually if the two contigs overlap to merge them together!)\n";
                                                            print OUTPUT4 " (Check manually if the two contigs overlap to merge them together!)\n";
                                                        }
                                                        else
                                                        {
                                                            print "\n";
                                                            print OUTPUT4 "\n";
                                                        }
                                                    }                                          
                                                    print OUTPUT ">Contig".$l."\n";
                                                    print OUTPUT $fin;
                                                    if (length($fin2) > $largest_contig)
                                                    {
                                                        $largest_contig = length($fin2);
                                                    }
                                                    if (length($fin2) < $miminum_contig)
                                                    {
                                                        $miminum_contig = length($fin2);
                                                    }
                                                    print "Contig ".$l."                  : ".length($fin2)." bp\n";
                                                    print OUTPUT4 "Contig ".$l."                  : ".length($fin2)." bp\n";
                                                    $assembly_length += length($fin2);
                                                }
                                                $assembly_success = "yes";
                                                if ($y > $startprint2)
                                                {                                                                                 
                                                    print OUTPUT5 ">".$id."\n";
                                                    print OUTPUT5 $read."\n\n\n";
                                                }
                                                
                                                print "\nTotal contigs              : ".$l."\n";
                                                print "Largest contig             : ".$largest_contig." bp\n";
                                                print "Smallest contig            : ".$miminum_contig." bp\n";
                                                print "Average insert size        : ".$insert_size." bp\n\n";
                                                print OUTPUT4 "\nTotal contigs              : ".$l."\n";
                                                print OUTPUT4 "Largest contig             : ".$largest_contig." bp\n";
                                                print OUTPUT4 "Smallest contig            : ".$miminum_contig." bp\n";
                                                print OUTPUT4 "Average insert size        : ".$insert_size." bp\n\n";
                                                $option++;
                                                close OUTPUT;
                                                $finish = "yes";
                                            }
                                            elsif ($circle ne "contigs")
                                            {
                                                my $output_file2  = "Uncircularized_assemblies_".$option."_".$project.".fasta";
                                                open(OUTPUT2, ">" .$output_file2) or die "Can't open file $output_file2, $!\n";
                                                
                                                if (length($read) > 250)
                                                {
                                                    $contigs{$contig_num."+".$id} = $read;
                                                    $contig_num++;
                                                }
                                                
                                                delete $seed{$id};
          
                                                $read =~ tr/\./N/;
                                                $read =~ tr/X//d;
                                                my @contigs = split /L+/, $read;
                                                my $l = '0';
                                                my $largest_contig = '0';
                                                my $miminum_contig = '100000000000000000000000000000';
                                                
                                                print "\b" x length($progress_before);
                                                print ' ' x length($progress_before);
                                                print "\b" x length($progress_before);
                                                print "\n\nAssembly ".$option." finished incomplete: The genome has not been circularized\n\n";
                                                print OUTPUT4 "\n\nAssembly ".$option." finished incomplete: The genome has not been circularized\n\n";
                                                print OUTPUT5 "\n\nAssembly ".$option." finished incomplete: The genome has not been circularized3\n\n";

                                                if ($y > $startprint2)
                                                {                                                                                 
                                                    print OUTPUT5 ">".$id."\n";
                                                    print OUTPUT5 $read."\n\n\n";
                                                }

                                                foreach my $contig_tmp (keys %contigs)
                                                {
                                                    if ($contig_tmp =~ m/(\d+)\+*\d*/)
                                                    {
                                                        my $contig_tmp3 = $1;
                                                        if ($contig_tmp3 < 10)
                                                        {
                                                            $contigs{"0".$contig_tmp} = $contigs{$contig_tmp};
                                                            delete $contigs{$contig_tmp};
                                                        }
                                                    }
                                                }
                                                foreach my $contig_tmp (sort keys %contigs)
                                                {
                                                    my $contig_tmp2;
                                                    if ($contig_tmp =~ m/(\d+)\+*\d*/)
                                                    {
                                                        $contig_tmp2 = $1;
                                                    }
                                                    $read = $contigs{$contig_tmp};
                                                    $read =~ tr/\./N/;
                                                    $read =~ tr/X//d;
                                                    my @contigs = split /L+/, $read;
                                                    my $jj = '0';
                                                    foreach (@contigs)
                                                    {
                                                        $l++;
                                                        $jj++;
                                                        my $fin = $_;
                                                        my $fin2 = $fin;
                                                        $fin =~ s/(.{1,150})/$1\n/gs;
                                                        
                                                        if ($jj > 1)
                                                        {
                                                            my $gap_min = sprintf("%.0f", $contig_gap_min{$id."_".($jj-1)});
                                                            my $gap_max = sprintf("%.0f", $contig_gap_max{$id."_".($jj-1)});
                                                            print "Estimated gap      : ".$gap_min." bp to ".$gap_max." bp";
                                                            print OUTPUT4 "Estimated gap      : ".$gap_min." bp to ".$gap_max." bp";
                                                            if ($contig_gap_min{$id."_".($jj-1)} < 0)
                                                            {
                                                                print " (Check manually if the two contigs overlap to merge them together!)\n";
                                                                print OUTPUT4 " (Check manually if the two contigs overlap to merge them together!)\n";
                                                            }
                                                            else
                                                            {
                                                                print "\n";
                                                                print OUTPUT4 "\n";
                                                            }
                                                        }                                          
                                                        print OUTPUT2 ">Contig".$contig_tmp."\n";
                                                        print OUTPUT2 $fin;
                                                        if (length($fin2) > $largest_contig)
                                                        {
                                                            $largest_contig = length($fin2);
                                                        }
                                                        if (length($fin2) < $miminum_contig)
                                                        {
                                                            $miminum_contig = length($fin2);
                                                        }
                                                        print "Contig ".$contig_tmp2."           : ".length($fin2)." bp\n";
                                                        print OUTPUT4 "Contig ".$contig_tmp2."           : ".length($fin2)." bp\n";
                                                    }
                                                    if ($y > $startprint2)
                                                    {                                                                                 
                                                        print OUTPUT5 ">Contig".$contig_tmp."\n";
                                                        print OUTPUT5 $read."\n";
                                                    }
                                                }
                                                print "\nTotal contigs          : ".$l."\n";
                                                print "Largest contig         : ".$largest_contig." bp\n";
                                                print "Smallest contig        : ".$miminum_contig." bp\n";
                                                print "Average insert size    : ".$insert_size." bp\n\n";
                                                print "----------------------------------------------------------------------------------------------------\n\n";
                                                print OUTPUT4 "\nTotal contigs          : ".$l."\n";
                                                print OUTPUT4 "Largest contig         : ".$largest_contig." bp\n";
                                                print OUTPUT4 "Smallest contig        : ".$miminum_contig." bp\n";
                                                print OUTPUT4 "Average insert size    : ".$insert_size." bp\n\n";
                                                print OUTPUT4 "----------------------------------------------------------------------------------------------------\n\n";
                                                $option++;
                                                close OUTPUT2;
                                                $finish = "yes";
                                            }
    }
}
$y++;
    
    if ($y eq $iterations)
    {
        foreach my $seedprint (keys %seed)
        {
            print OUTPUT5 "\n".$seedprint."\n";
            print OUTPUT5 $seed{$seedprint}."\n\n";
            
            my $id_match_end = substr $seedprint, -1, 1,"",;
            if ($id_match_end eq "1")
            {
                $id_pair_match = $seedprint."2";
            }
            elsif ($id_match_end eq "2")
            {
                $id_pair_match = $seedprint."1";
            }
        }
    }
    if (@insert_size > 500 && @insert_size < 8000 && $insert_size_correct eq "yes")
    {          
        my $insert_total = '0';
        my $k = '1';
        foreach my $insert1 (@insert_size)
        {
            $insert_total += $insert1;
            $k++;
        }
        
        my $insert_size_temp = $insert_total/$k;
        $insert_size = int($insert_size_temp + $insert_size_temp/abs($insert_size_temp*2));
        if ($y > $startprint2)
        {
            print OUTPUT5 $insert_size." Insert Size\n";
        }
    }
    elsif (@insert_size >= 8000)
    {
        $insert_size_correct = "";
        $insert_range = $insert_range_b;
        $insert_range_back = $insert_range_b;
    }
    elsif ($insert_size_correct ne "yes" && length($read) > $insert_size)
    {
        $insert_range = $insert_range_b;
        $insert_range_back = $insert_range_b;
    }
}
FINISH2:

    if ($finish ne "yes")
    {
        foreach my $seed_id_tmp (keys %seed)
        {
            if (length($seed{$seed_id_tmp}) > 250)
            {
                $contigs{$contig_num."+".$seed_id_tmp} = $seed{$seed_id_tmp};
                $contig_num++; 
            }
            delete $seed{$seed_id_tmp};
        }
                                                
        foreach my $contig_tmp (keys %contigs)
        {
            if ($contig_tmp =~ m/(\d+)\+*\d*/)
            {
                my $contig_tmp3 = $1;
                if ($contig_tmp3 < 10)
                {
                    $contigs{"0".$contig_tmp} = $contigs{$contig_tmp};
                    delete $contigs{$contig_tmp};
                }
            }
        }
    }
    my $count_contigs = keys %contigs;

    my $tree_succes = "";
    if($circle ne "yes" && $count_contigs > 1)
    {
        my $h = '0';
        my $terminate = '0';
        my %node;
        my %row_nodes;
        my %contig_num;
        my %contigs2;
        undef %contigs2;
        my %repetitive;
        undef %repetitive;
        foreach my $contig_tmp (keys %contigs)
        {
            if ($contig_tmp =~ m/(\d+)\+*(.*)/)
            {
                my $contig_num = $1;
                my $contig_code = $2;
                if ($contig_code =~ m/.*_(\d+)$/)
                {
                    $contig_code = $1;
                }
                $contig_num{$contig_code} = $contig_num;
                $contigs2{$contig_num} = $contigs{$contig_tmp};
            }
        }
        print OUTPUT7 "LINKS BETWEEN CONTIGS\n";
        print OUTPUT7 "---------------------\n\n";
        foreach my $tree (keys %tree)
        {
            my $tree2 = $tree;
            my $tree3 = $tree{$tree2};
            if ($tree2 =~ m/.*_(\d+(REP)*)$/)
            {
                $tree = $1;
            }
            if ($tree3 =~ m/.*_(\d+(REP))$/)
            {
                $tree3 = $1;
            }
            delete $tree{$tree2};
            $tree{$tree} = $tree3;
            my $tree_tmp = $tree;
            my $tree_tmp2 = $tree{$tree};

            my @ids_split = split /\*/, $tree_tmp2;
            foreach my $id_split (@ids_split)
            {
                foreach my $contig_num (keys %contig_num)
                {
                    if ($id_split  =~ m/^$contig_num(REP)*$/)
                    {
                        if ($tree_tmp2 =~ m/^(.*\*)*$contig_num(REP)*(\*.*)*$/)
                        {
                            if (defined($1))
                            {
                                $tree_tmp2 = $1.$contig_num{$contig_num};
                            }
                            else
                            {
                                $tree_tmp2 = $contig_num{$contig_num};
                            }
                            if (defined($2))
                            {
                                $tree_tmp2 = $tree_tmp2."REP";
                            }
                            if (defined($3))
                            {
                                $tree_tmp2 = $tree_tmp2.$3;
                            }
                        }
                    }
                    if ($contig_num eq $tree_tmp)
                    {
                        $tree_tmp = $contig_num{$contig_num};
                    }
                }
            }
            my $rep_test = substr $tree_tmp2, -3;
            $tree_tmp2 =~ s/\*/ OR /g;
            if ($rep_test eq "REP")
            {
                my $tree_tmp2a = substr $tree_tmp2,0, -3;
                print OUTPUT7 $tree_tmp."----> REPETITIVE REGION----> ".$tree_tmp2a."\n";
            }
            else
            {
                print OUTPUT7 $tree_tmp."----> ".$tree_tmp2."\n";
            }
        }
        if (exists($tree{"START"}))
        {
            $row{$h} = "01";
            $node{$h} = $tree{"START"};
        }
        if ($hasL eq "yes")
        {
            print OUTPUT7 "\n(Contigs broken up by long homopolymer stretches are linked together as one contig with 15 N's, they can be merged manually in some cases)\n";
            print OUTPUT7 "\n(If the region before and after the N zone overlaps, you can merge them by deleting the regions next to N zone from both sides (those are the least reliable))\n\n";

        }                                     
        my %row_circle;
        undef %row_circle;
TERMINATE:                                      while (keys %node)
        {
            foreach my $h1 (keys %node)
            {
                undef %row_nodes;
                %row_nodes = map { $_ => undef } split(/\+/, $row{$h1});
                my @row_nodes = map { $_ => undef } split(/\+/, $row{$h1});
                if(exists($tree{$node{$h1}}))
                {
                    if ($tree{$node{$h1}} eq $node{$h1})
                    {
                        delete $node{$h1};
                        delete $row{$h1};
                    }
                    elsif ($tree{$node{$h1}} eq "END_REVERSE")
                    {
                        delete $node{$h1};
                        delete $row{$h1};
                    }
                    elsif ($tree{$node{$h1}} eq "END")
                    {
                        delete $node{$h1};
                        $row_circle{$h1} = $row{$h1};
                    }
                    elsif ($tree{$node{$h1}} =~ m/(.*)\*(.*)/)
                    {                                                          
                        my @id_node = split /\*/, $tree{$node{$h1}};
                        foreach my $id_node (@id_node)
                        {
                            my $count1 = '0';
                            foreach my $row_n (@row_nodes)
                            {
                                if ($id_node eq $row_n)
                                {
                                    $count1++;
                                }                                                                   
                            }
                            if ($id_node eq $node{$h1})
                            {
                            }
                            elsif ($count1 > 0)
                            {
                            }
                            else
                            {
                                $h++;
                                $node{$h} = $id_node;
                                $row{$h} = $row{$h1}."+".$id_node;
                            }
                        }
                        delete $node{$h1};
                        delete $row{$h1};
                    }
                    elsif ($tree{$node{$h1}} =~ m/(.*)REP/)
                    {
                        $h++;       
                        $row{$h} = $row{$h1}."+".$1."R";
                        $node{$h} = $1;
                        delete $node{$h1};
                        delete $row{$h1};
                        $repetitive{$1} = undef;
                    }
                    else
                    {
                        $h++;
                        if ($row{$h1} =~ m/^(.*)\+\d+$/)
                        {
                            $row{$h1} = $1;
                        }
                        if ($row{$h1} =~ m/^01$/)
                        {
                            $row{$h} = $row{$h1}
                        }
                        else
                        {
                            $row{$h} = $row{$h1}."+".$tree{$node{$h1}};
                        }
                        $node{$h} = $tree{$node{$h1}};
                        delete $node{$h1};
                        delete $row{$h1};
                    }
                }
                else
                {
                   delete $node{$h1};
                   delete $row{$h1};
                }
                $terminate++;
                if ($terminate > 1500)
                {
                    delete $row{$h1};
                    last TERMINATE;
                } 
            }
        }
        my $g = '1';

        foreach my $row (keys %row_circle)
        {
            foreach my $contig_num (keys %contig_num)
            {
                $row{$row} =~ s/\+$contig_num\+/\+$contig_num{$contig_num}\+/g;
                $row{$row} =~ s/\+$contig_num$/\+$contig_num{$contig_num}/g;
                $row{$row} =~ s/\+$contig_num(R)\+/\+$contig_num{$contig_num}R\+/g;
                $row{$row} =~ s/\+$contig_num(R)$/\+$contig_num{$contig_num}R/g;
            }
            
            my @row = split /\+/,$row{$row};
            my $assembly = "";
            
            foreach my $cont (@row)
            {
                my $check = $cont =~ tr/R//d;
                if (exists($contigs2{$cont}))
                {
                    my $repe2 = substr $assembly, -1;
                    if ($cont eq '1' || $cont eq '01')
                    {
                        $assembly = $contigs2{$cont};
                    }
                    elsif ($check > 0)
                    {
                        $assembly .= "RRRRRRRRRRRRRR".$contigs2{$cont};
                    }
                    else
                    {
                        my $end_assembly = substr $assembly, -30;
                        $end_assembly =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        my $start_next_contig = substr $contigs2{$cont},0 ,500;
                        $start_next_contig =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
                        my $x1 = length($start_next_contig);
                        $start_next_contig =~ s/.*$end_assembly//;
                        my $x2 = length($start_next_contig);
                        if ($x1-$x2 ne '0')
                        {
                            my $assembly_tmp = substr $assembly, 0, -($x1-$x2);
                            $assembly = $assembly_tmp.$contigs2{$cont};
                        }
                        else
                        {
                            $assembly .= "2RRRRRRRRRRRRRR".$contigs2{$cont};
                        }  
                    }
                }
            }
            my $start_assembly = substr $assembly, 50, 30;
            $start_assembly =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
            my $end_last_contig = substr $assembly, -1000;
            $end_last_contig =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
            my $x1 = length($end_last_contig);
            $end_last_contig =~ s/(.)$start_assembly.*/$1/;
            my $x2 = length($end_last_contig);
            my $assembly2;
            if (-$x1+$x2 ne '0')
            {
                $assembly2 = substr $assembly, 50, -$x1+$x2;
            }
            else
            {
                $assembly2 = $assembly;
            }  
            $assembly = $assembly2;
            $assembly =~ tr/\./N/;
            
            if (length($assembly) > 500)
            {
                print OUTPUT7 "\n\nOPTION ".$g."\n-------------------\n";
                print OUTPUT7 "Contig Arrangement : ";
                print OUTPUT7 $row{$row}."\n";
                print OUTPUT7 "Assembly length    : ".length($assembly)." bp\n\n";
                my @contigs0 = split /RRRRRRRRRRRRRR/, $assembly;
                my $size_c = @contigs0;
                my $jj = '0';
                my $t ='0';
                my @order = split /\+/, $row{$row};
                foreach (@contigs0)
                {
                    $jj++;
                    my $fin = $_;
                    $fin =~ tr/L/N/;
                    my $order = "";
                    my $h = '0';
ORDER:              foreach my $tmp (@order)
                    {
                        if ($h < $t)
                        {
                            $h++;
                            next ORDER;
                        } 
                        my $R = $tmp =~ s/RRRRRRRRRRRRRR//g;
                        if ($R ne '1'&& $h > 0)
                        {
                            $order .= "+".$tmp;
                        }
                        elsif ($R ne '1' && $h eq '0')
                        {
                            $order = $tmp;
                        }
                        if ($R > 0 && $h ne '0')
                        {
                            last ORDER;
                        }
                        $t++;
                        $h++;
                    }
                    my $gg = substr $order, 0, 1;
                    if ($gg eq "+")
                    {
                        substr $order, 0, 1, "";
                    }
                    
                    
                    print OUTPUT7 ">Contig ".$order."\n";
                    my $fail = substr $fin, -1;
                    if ($fail eq '2')
                    {
                        substr $fin, -1, 1, "";
                    }
                    print OUTPUT7 $fin."\n";
                    
                    if ($size_c > 1 && $jj < $size_c && $fail ne '2')
                    {
                        print OUTPUT7 "\n-------Repetitive region detected, exact length unknown-------\n\n";
                    }
                    elsif ($size_c > 1 && $jj < $size_c && $fail eq '2')
                    {
                        print OUTPUT7 "\n-------Couldn't merge automatically, try manually-------\n\n";
                    }
                }
                 
                print OUTPUT7 "\n";
                
                my $output_file9  = "Option_".$g."_".$project.".fasta";
                open(OUTPUT9, ">" .$output_file9) or die "Can't open file $output_file9, $!\n";
                $assembly =~ tr/L/N/;
                my @contigs2 = split /RRRRRRRRRRRRRR/, $assembly;
                my $ww = '0';
                foreach (@contigs2)
                {
                    $ww++;
                    my $fin = $_;
                    my $fin2 = $fin;
                    $fin =~ s/(.{1,150})/$1\n/gs;
                                                             
                    print OUTPUT9 ">Contig".$ww."\n";
                    print OUTPUT9 $fin;
                }
                close OUTPUT9;
                $g++;
                
                if (length($assembly) > $genome_range_low && length($assembly) < $genome_range_high)
                {
                    $tree_succes = "yes";
                    $assembly_length = length($assembly);
                }
            }                                                   
        }
        print OUTPUT7 "\n\nEach option has a separate fasta file\n";
    }
    elsif ($count_contigs eq 1)
    {
        foreach my $contig_tmp1 (keys %contigs)
        {
            $assembly_length = length($contigs{$contig_tmp1});
        } 
    }
    if ($finish ne "yes")
    {
                                                my $output_file2  = "Contigs_".$option."_".$project.".fasta";
                                                open(OUTPUT2, ">" .$output_file2) or die "Can't open file $output_file2, $!\n";
                                    
                                                my $l = '0';
                                                my $largest_contig = '0';
                                                my $miminum_contig = '1000000000000000000000000000000000';
                                                
                                                print "\b" x length($progress_before);
                                                print ' ' x length($progress_before);
                                                print "\b" x length($progress_before);
                                                print "\n------------Assembly ".$option." finished: Contigs are automatically merged in Merged_contigs file------------\n\n";
                                                print OUTPUT4 "\n------------Assembly ".$option." finished: Contigs are automatically merged in Merged_contigs file------------\n\n";
                                                print OUTPUT5 "\n------------Assembly ".$option." finished: Contigs are automatically merged in Merged_contigs file------------\n\n";                                     
                                                
                                                
                                                foreach my $contig_tmp (sort keys %contigs)
                                                {
                                                    my $contig_tmp2;
                                                    if ($contig_tmp =~ m/(\d+)\+*\d*/)
                                                    {
                                                        $contig_tmp2 = $1;
                                                    }
                                                    
                                                    $read = $contigs{$contig_tmp};
                                                    $read =~ tr/\./N/;
                                                    $read =~ tr/X//d;
                                                    my @contigs = split /L+/, $read;
                                                    my $jj = '0';
                                                    
                                                    foreach (@contigs)
                                                    {
                                                        $l++;
                                                        $jj++;
                                                        my $fin = $_;
                                                        my $fin2 = $fin;
                                                        $fin =~ s/(.{1,150})/$1\n/gs;
                                                        
                                                        if ($jj > 1)
                                                        {
                                                            my $gap_min = sprintf("%.0f", $contig_gap_min{$id."_".($jj-1)});
                                                            my $gap_max = sprintf("%.0f", $contig_gap_max{$id."_".($jj-1)});
                                                            print "Estimated gap      : ".$gap_min." bp to ".$gap_max." bp";
                                                            print OUTPUT4 "Estimated gap      : ".$gap_min." bp to ".$gap_max." bp";
                                                            if ($contig_gap_min{$id."_".($jj-1)} < 0)
                                                            {
                                                                print " (Check manually if the two contigs overlap to merge them together!)\n";
                                                                print OUTPUT4 " (Check manually if the two contigs overlap to merge them together!)\n";
                                                            }
                                                            else
                                                            {
                                                                print "\n";
                                                                print OUTPUT4 "\n";
                                                            }
                                                        }                                          
                                                        print OUTPUT2 ">Contig".$contig_tmp."\n";
                                                        print OUTPUT2 $fin;
                                                        if (length($fin2) > $largest_contig)
                                                        {
                                                            $largest_contig = length($fin2);
                                                        }
                                                        if (length($fin2) < $miminum_contig)
                                                        {
                                                            $miminum_contig = length($fin2);
                                                        }
                                                        print "Contig ".$contig_tmp2."                  : ".length($fin2)." bp\n";
                                                        print OUTPUT4 "Contig ".$contig_tmp2."                  : ".length($fin2)." bp\n";
                                                    }
                                                    if ($y > $startprint2)
                                                    {                                                                                 
                                                        print OUTPUT5 ">Contig".$contig_tmp."\n";
                                                        print OUTPUT5 $read."\n";
                                                    }
                                                }
                                                
                                                print "\nTotal contigs              : ".$l."\n";
                                                print "Largest contig             : ".$largest_contig." bp\n";
                                                print "Smallest contig            : ".$miminum_contig." bp\n";
                                                print "Average insert size        : ".$insert_size." bp\n\n";
                                                print OUTPUT4 "\nTotal contigs              : ".$l."\n";
                                                print OUTPUT4 "Largest contig             : ".$largest_contig." bp\n";
                                                print OUTPUT4 "Smallest contig            : ".$miminum_contig." bp\n";
                                                print OUTPUT4 "Average insert size        : ".$insert_size." bp\n\n";
                                                $option++;
                                                
                                                close OUTPUT2
    }



                                                my $total_reads_organelle = (keys %count_reads)*2;
                                                my $total_reads_organelle_all = (keys %count_reads_all)*2;
                                                my $total_reads = (keys %hash)*2;
                                                my $organelle_percentage = sprintf("%.2f",($total_reads_organelle_all*100)/$total_reads);
                                                my $average_coverage = sprintf("%.0f",$total_reads_organelle_all*$read_length/$assembly_length);
                                                
                                                print "\n-----------------------------------------Input data metrics-----------------------------------------\n\n";
                                                print "Total reads                : ".$total_reads."\n";
                                                print "Aligned reads              : ".$total_reads_organelle_all."\n";
                                                print "Assembled reads            : ".$total_reads_organelle."\n";
                                                if ($assembly_success eq "yes" ||  $tree_succes eq "yes")
                                                {
                                                   print "Organelle genome %         : ".$organelle_percentage." %\n"; 
                                                }
                                                if ($assembly_length > 5000)
                                                {
                                                   print "Average organelle coverage : ".$average_coverage."\n\n"; 
                                                }
                                                print "----------------------------------------------------------------------------------------------------\n\n";
                                                print OUTPUT4 "\n-----------------------------------------Input data metrics-----------------------------------------\n\n";
                                                print OUTPUT4 "Total reads                : ".$total_reads."\n";
                                                print OUTPUT4 "Aligned reads              : ".$total_reads_organelle_all."\n";
                                                print OUTPUT4 "Assembled reads            : ".$total_reads_organelle."\n";
                                                if ($assembly_success eq "yes" ||  $tree_succes eq "yes")
                                                {
                                                   print OUTPUT4 "Organelle genome %         : ".$organelle_percentage." %\n"; 
                                                }
                                                if ($assembly_length > 5000)
                                                {
                                                   print OUTPUT4 "Average organelle coverage : ".$average_coverage."\n\n"; 
                                                }
                                                print OUTPUT4 "----------------------------------------------------------------------------------------------------\n\n";
print "\nThank you for using NOVOPlasty!\n\n";
close INPUT;
close OUTPUT4;
close OUTPUT5;
close OUTPUT6;
close OUTPUT7;
