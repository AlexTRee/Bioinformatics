#!/usr/bin/perl
use strict;
use warnings;

`rm extra_tail_list.txt`;
`rm maid_oligo_data.fa`;

open (FILE, "maid_oligo_data.txt") || die "Can't open data file!\n";
my @lines = <FILE>;
close (FILE);

foreach my $line (@lines)
{
	my @arr = split ("\t", $line); 
	my $maid = $arr[0];
	my $oligo = $arr[1];
	my $seq = $arr[2];
	chomp $seq;
	my $len = length($seq);
	if ((substr $oligo, 0, 3) eq "HUF")
	{
		my $tail1 = substr $seq, 0, 15;
		my $tail2 = substr $seq, 0, 14;
		my $tail3 = substr $seq, 0, 12;
		my $tail4 = substr $seq, 0, 13;
		my $tail5 = substr $seq, 0, 20;
		my $tail6 = substr $seq, 0, 11;
		my $tail7 = substr $seq, 0, 17;
		
		if (($tail5 eq "AAGCTTACTCGAGCCTGCAG") || ($tail5 eq "TCATAGCTCAGCCTTTCTCC"))
		{
			$seq = substr $seq, 20;
		}
		elsif ($tail7 eq "AGCAAGCTTGCGGCCGC")
		{
			$seq = substr $seq, 17;
		}		
		elsif (($tail1 eq "AGCAAGCTTCTGCAG") || ($tail1 eq "AGCAAGCTTGAGCTC") || ($tail1 eq "AGCAAGCTTATGCAT"))
		{
			$seq = substr $seq, 15;
		}
		elsif (($tail2 eq "GCAAGCTTCTGCAG") || ($tail2 eq "GCAAGCTTGAGCTC") || ($tail2 eq "CAGCTTGCGGCCGC"))
		{
			$seq = substr $seq, 14;
		}
		elsif (($tail4 eq "TCCTTGCGGCCGC") || ($tail4 eq "CCTCCGCGATCGC") || ($tail4 eq "GCCAAGCGGCCGC"))
		{
			$seq = substr $seq, 13;
		}
		elsif ($tail3 eq "AAGCTTCTGCAG")
		{
			$seq = substr $seq, 12;
		}
		elsif ($tail6 eq "CCTCCCTGCAG")
		{
			$seq = substr $seq, 11;
		}
			
		elsif ($len > 26)
		{
			open (OUTPUT, '>>extra_tail_list.txt');
			print OUTPUT ">$maid"."$oligo\n$seq\n";
		}
	}
	
	elsif ((substr $oligo, 0, 3)  eq "HDR")
	{
		my $tail1 = substr $seq, 0, 15;
		my $tail2 = substr $seq, 0, 12;
		my $tail3 = substr $seq, 0, 14;
		my $tail4 = substr $seq, 0, 11;
		my $tail5 = substr $seq, 0, 20;
		
		if ($tail5 eq "GGGCTAGCAAGCTTGAGCTC")
		{
			$seq = substr $seq, 20;
		}
		elsif (($tail1 eq "AGCAAGCTTGAGCTC") || ($tail1 eq "AGCAAGCTTCTGCAG"))
		{
			$seq = substr $seq, 15;
		}
		elsif (($tail3 eq "CCGTGGAGCTCGAG") || ($tail3 eq "GCTTGACGTGAGTC"))
		{
			$seq = substr $seq, 14;
		}
		elsif (($tail2 eq "CTCTTGGCGCGC") || ($tail2 eq "GCCAGTGAGCTC"))
		{
			$seq = substr $seq, 12;
		}
		elsif (($tail4 eq "CCTCCATCGAT") || ($tail4 eq "CCTCCCTCGAG")|| ($tail4 eq "CCTCCGAGCTC"))
		{
			$seq = substr $seq, 11;
		}
		
		elsif ($len > 26)
		{
			open (OUTPUT, '>>extra_tail_list.txt');
			print OUTPUT ">$maid"."$oligo\n$seq\n";
		}
	}
	elsif ((substr $oligo, 0, 3) eq "HUR")
	{
		my $tail1 = substr $seq, 0, 11;
		my $tail2 = substr $seq, 0, 12;
		my $tail3 = substr $seq, 0, 15;
		my $tail4 = substr $seq, 0, 16;
		my $tail5 = substr $seq, 0, 17;
		my $tail6 = substr $seq, 0, 19;
		my $tail7 = substr $seq, 0, 23;
		my $tail8 = substr $seq, 0, 20;
		my $tail9 = substr $seq, 0, 21;
		my $tail10 = substr $seq, 0, 25;
		my $tail11 = substr $seq, 0, 22;
		my $tail12 = substr $seq, -19;
		my $tail13 = substr $seq, 0, 13;
		
		if ($tail12 eq "ATCTGAAGTGCCATGATCC") 		
		{
			$seq = substr $seq, -19;
		}
		elsif ($tail10 eq "CGTGTCTCCTGCAGCCCCAGGTACC")
		{
			$seq = substr $seq, 25;
		}
		elsif ($tail7 eq "TCCTCAACTGGGATGATGGTACC")
		{
			$seq = substr $seq, 23;
		}
		elsif ($tail11 eq "GGCACCAGAATGTACCACTGGG")
		{
			$seq = substr $seq, 22;
		}
		elsif (($tail9 eq "CCTCCCCCGGGGGGCCTGCAG") || ($tail9 eq "CCTCCCTCGAGGGGCCTGCAG"))
		{
			$seq = substr $seq, 21;
		}
		elsif ($tail8 eq "GGATCATTTAAATCGGTACC")
		{
			$seq = substr $seq, 20;
		}
		elsif ($tail6 eq "GGCTGGCAGGGGGGCATGC")
		{
			$seq = substr $seq, 19;
		}
		elsif (($tail5 eq "CGACGTCGTCTCGTCAC") || ($tail5 eq "CCTCCCCCGGGCTCGAG") || ($tail5 eq "GGCACCAGAATGTACCA") || ($tail5 eq "GTCGGAGCACCGCTCAT"))
		{
			$seq = substr $seq, 17;
		}
		elsif (($tail4 eq "GTGACGAGACGACGTC")|| ($tail4 eq "ACGAAGTTATCTCGAG") || ($tail4 eq "ACGAAGTTATGCTAGG"))
		{
			$seq = substr $seq, 16;
		}
		elsif (($tail3 eq "AGCAAGCTTGAGCTC") || ($tail3 eq "GGTCCCCAAACTCAC"))
		{
			$seq = substr $seq, 15;
		}
		elsif (($tail2 eq "CCTTCCGGTACC") || ($tail2 eq "GAACTGCCCGGG") || ($tail2 eq "ACTGGACCCGGG"))
		{
			$seq = substr $seq, 12;
		}
		elsif ($tail13 eq "CCTCCGCGATCGC")
		{
			$seq = substr $seq, 13;
		}
		elsif (($tail1 eq "CCTCCGGTACC") || ($tail1 eq "CCTCCCCTAGG") || ($tail1 eq "CCTCTGGTACC") || ($tail1 eq "GAGGAGGTACC") || ($tail1 eq "CCTCCTGTACA")|| ($tail1 eq "CCTCCCTCGAG") || ($tail1 eq "CCTCCCCCGGG") || ($tail1 eq "CCTCCACGCGT")|| ($tail1 eq "CCTCCGAGCTC")|| ($tail1 eq "CCTCCGCGATC"))
		{
			$seq = substr $seq, 11;
		}
		
		elsif ($len > 26)
		{
			open (OUTPUT, '>>extra_tail_list.txt');
			print OUTPUT ">$maid"."$oligo\n$seq\n";
		}
	}
	elsif ((substr $oligo, 0, 3) eq "HDF")
	{
		my $tail1 = substr $seq, 0, 11;
		my $tail2 = substr $seq, 0, 10;
		my $tail3 = substr $seq, 0, 15;
		my $tail4 = substr $seq, 0, 12;
		my $tail5 = substr $seq, 0, 16;
		my $tail6 = substr $seq, 0, 14;
		
		if (($tail5 eq "ACGAAGTTATCTCGAG") || ($tail5 eq "ACGAAGTTATGCTAGG"))
		{
			$seq = substr $seq, 16;
		}
		elsif (($tail3 eq "AGCAAGCTTCTGCAG") || ($tail3 eq "AGCAAGCTTGAGCTC") || ($tail3 eq "TTCTCTTTCCTACAG"))
		{
			$seq = substr $seq, 15;
		}
		elsif ($tail6 eq "ACGAAGTTATGCTA")
		{
			$seq = substr $seq, 14;
		}
		elsif (($tail4 eq "GGACATCGTCTC") || ($tail4 eq "GAACTGCCCGGG"))
		{
			$seq = substr $seq, 12;
		}
		elsif (($tail1 eq "CCTCCCCTAGG") || ($tail1 eq "CCTCCGCTAGC") || ($tail1 eq "CCTCCCTCGAG") || ($tail1 eq "CCTCCCCCGGG") || ($tail1 eq "CCTCCCTGCAG"))
		{
			$seq = substr $seq, 11;
		}
		elsif ($tail2 eq "CTCCCCTAGG")
		{
			$seq = substr $seq, 10;
		}
			
		elsif ($len > 26)
		{
			open (OUTPUT, '>>extra_tail_list.txt');
			print OUTPUT ">$maid"."$oligo\n$seq\n";
		}
	}
	else
	{
		open (OUTPUT, '>>maid_oligo_data.fa');
		print OUTPUT ">$maid"."$oligo\n$seq\n";
		next;
	}
	open (OUTPUT, '>>maid_oligo_data.fa');
	print OUTPUT ">$maid"."$oligo\n$seq\n";
	
}
