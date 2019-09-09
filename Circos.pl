#!/usr/bin/perl -w
use strict;

my $input_file  = "Circos_links_mtDNA_1:100_hp_new3.txt";
my $output_file = "Circos_links_mtDNA_1:100_hp_new3b.txt"; 
my $full = "yes";
my $half = "";
my $not = "yes";
my $range_low = 0;
my $range_high = 17000;


open(INPUT, $input_file) or die "Can't open file $input_file, $!\n";
open(OUTPUT, ">" .$output_file) or die "Can't open file $output_file, $!\n";
my $r = '0';

my %first_selection_full;
my %final_selection_full;
my %first_selection_half;
my %final_selection_half;
my %first_selection_not;
my %final_selection_not;
my $list = "";
my @list = (2167,2335,2458,2557,2668,2737);
my %list;
foreach my $l (@list)
{
    $list{$l} = undef;
}
my $g = '0';

while (my $line = <INPUT>)
{
    chomp ($line);
    $r++;
    if ($line =~ m/^hsM\s(\d+)\s\d+\shsM\s(\d+).*color=240,190,0$/ && $half eq "yes")
    {
        my $first = $1;
        my $second = $2;
        if ($first < $range_high && $first > $range_low && $second < $range_high && $second > $range_low)
        {
            #print OUTPUT $line."\n";
            if (exists($list{$first}))
            {
                $final_selection_half{$g} = $line;
            }
            elsif (exists($list{$second}))
            {
                $final_selection_half{$g} = $line;
            }
            $g++;
        }
    }
    elsif ($line =~ m/^hsM\s(\d+)\s\d+\shsM\s(\d+).*color=0,114,178$/ && $full eq "yes")
    {
        my $first = $1;
        my $second = $2;
        if ($first < $range_high && $first > $range_low && $second < $range_high && $second > $range_low)
        {
            #print $line."\n";
            my $tmp = $second." ".$first;
            my $tmp2 = $first." ".$second;
            #if (exists($first_selection_full{$tmp}))
            if (exists($first_selection_full{$tmp}) && (exists($list{$first}) || exists($list{$second}) || $list eq ""))
            {
               $final_selection_full{$tmp} = $line;
            }
            else
            {
                $first_selection_full{$tmp2} = $line;
            }
        }
    }
    elsif ($line =~ m/^hsM\s(\d+)\s\d+\shsM\s(\d+).*color=213,94,0$/ && $not eq "yes")
    {
        my $first = $1;
        my $second = $2;
        if ($first < $range_high && $first > $range_low && $second < $range_high && $second > $range_low)
        {
            #print $line."\n";
            my $tmp = $second." ".$first;
            my $tmp2 = $first." ".$second;
            if (exists($final_selection_not{$tmp}))
            {
               #$final_selection_not{$tmp} = $line;
            }
            elsif (exists($final_selection_not{$tmp2}))
            {
            }
            else
            {
                #$first_selection_not{$tmp2} = undef;
                $final_selection_not{$tmp} = $line;
            }
        }
    }
}
foreach my $keys (keys %first_selection_full)
{
    my $first = "";
    my $second = "";   
    if ($keys =~ m/^(\d+)\s(\d+)$/)
    {
        $first = $1;
        $second = $2;
    }
    my $tmp = $second." ".$first;
    my $tmp2 = $first." ".$second;
    if (exists($final_selection_full{$tmp}))
    {
    }
    elsif (exists($final_selection_full{$tmp2}))
    {
    }
    else
    {
        $final_selection_full{$keys} = $first_selection_full{$keys};
    }
}
foreach my $keys (keys %final_selection_full)
{
    print OUTPUT $final_selection_full{$keys}."\n";
}
foreach my $keys (keys %final_selection_half)
{
    print OUTPUT $final_selection_half{$keys}."\n";
}
foreach my $keys (keys %final_selection_not)
{
    print OUTPUT $final_selection_not{$keys}."\n";
}
close INPUT;
close OUTPUT;
