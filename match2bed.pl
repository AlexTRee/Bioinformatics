#!/usr/bin/perl

########################################################
# Author: Tiange Cui                                   #
# Date: August-07-2015                                 #
# Funtion: Convert BLAT parser results into .BED files.# 
# Arguments: None                                      #
# Execute: perl match2bed.pl                           #
# Requirements: None                                   #
########################################################

use strict;
use List::Util qw(first);

`rm taqman.txt`;
`rm recomb_oligos.txt`;

my $filename = "match.txt";
open (FILE, $filename)|| die "Can't open data file!\n";
my @lines = <FILE>;

# sort first by MAID numerically and then oligos alphabetically
my @lines = sort {(split(/\t/,$a))[5]<=>(split(/\t/,$b))[5] || (split(/\t/,$a))[6] cmp (split(/\t/,$b))[6]} @lines;

# use only unique lines, filter out the duplicate ones.
my %unique = ();
foreach my $item (@lines)
{
	$unique{$item}++;
}
my @uniqarr = keys %unique;

# sort first by MAID numerically and then oligos alphabetically
my @lines2 = sort {(split(/\t/,$a))[5]<=>(split(/\t/,$b))[5] || (split(/\t/,$a))[6] cmp (split(/\t/,$b))[6]} @uniqarr;

for (my $n= 0; $n < @lines2; $n++)
{
	my @arr = split ("\t", $lines2[$n]);
	my $oligo = $arr[6];
	chomp $oligo;
	
	# Output results in 2 files, TaqMan assay and Recombination oligos.
	if ((substr $oligo, 0, 1) eq "T")
	{
		open (OUTPUT1, '>>taqman.txt');
		print OUTPUT1 "chr".$arr[0]."\t".$arr[2]."\t".$arr[3]."\t".$arr[5].$oligo."\t".$arr[4]."\t".$arr[1]."\n";
	}
	else
	{
		open (OUTPUT2, '>>recomb_oligos.txt');
		print OUTPUT2 "chr".$arr[0]."\t".$arr[2]."\t".$arr[3]."\t".$arr[5].$oligo."\t".$arr[4]."\t".$arr[1]."\n";
	}
}
close (OUTPUT1);
close (OUTPUT2);
close (FILE);

`rm deletion.txt`;

my $filename = "recomb_oligos.txt";
open (FILE, $filename)|| die "Can't opennnn data file!\n";
my @lines = reverse <FILE>;
my @temparr = ();
my @temparr2 = ();
my @maid = ();
my $line;
my $curlinenum = 0;
my $cnt = 0;
my $cnt2 = 1;
my $index = 0;
my ($hur, $hdf, $hu, $hd);

my @firstline = split ("\t", $lines[0]);
my $chrf = $firstline[0];
my $startf = $firstline[1];
my $endf = $firstline[2];
my @maid_oligof = grep {$_} split /(\d+)/, $firstline[3];
my $maidf = $maid_oligof[0];
my $oligof = $maid_oligof[1];
my $oligof_num = 0;
if (defined $maid_oligof[2])
{
	$oligof_num = $maid_oligof[2];
	$oligof = $oligof.$oligof_num;
}

my $strf = $firstline[5];
chomp $strf;

if (($oligof_num) && ((substr $oligof, 0, 3) eq "HUR"))
{
	if ((!$hur) && ($oligof_num == 4))
	{
		push (@temparr, $startf);
		push (@temparr, $endf);
		$hur = 1;
	}
	if ((!$hur) && ($oligof_num == 3))
	{
		push (@temparr, $startf);
		push (@temparr, $endf);
		$hur = 1;
	}
	if ((!$hur) && ($oligof_num == 2))
	{
		push (@temparr, $startf);
		push (@temparr, $endf);
		$hur = 1;
	}
	if ((!$hur) && ($oligof_num == 1))
	{
		push (@temparr, $startf);
		push (@temparr, $endf);
		$hur = 1;
	}
}
elsif ((substr $oligof, 0, 2) eq "HU")
{
	if ((!$hur) && ($oligof eq "HUR"))
	{
		push (@temparr, $startf);
		push (@temparr, $endf);
		$hur = 1;
	}
	elsif ((!$hu) && ($oligof eq "HU"))
	{
		push (@temparr2, $startf);
		push (@temparr2, $endf);
		$hu = 1;
	}
	elsif ((!$hu) && ($oligof eq "HU2"))
	{
		push (@temparr2, $startf);
		push (@temparr2, $endf);
		$hu = 1;
	}
}

for (my $n = 1; $n < @lines; $n++)
{
	my @arr = split ("\t", $lines[$n]);
	my $chr = $arr[0];
	my $start = $arr[1];
	my $end = $arr[2];
	my @maid_oligo = grep {$_} split /(\d+)/, $arr[3];
	my $maid = $maid_oligo[0];
	my $oligo = $maid_oligo[1];
	my $oligo_num = 0;
	if (defined $maid_oligo[2])
	{
		$oligo_num = $maid_oligo[2];
		$oligo = $oligo.$oligo_num;
	}
	my $score = $arr[4];
	chomp $score;
	my $str = $arr[5];
	chomp $str;
	if ($maid eq $maidf) # 100  ->  100
	{
		$maidf = $maid;
		$oligof = $oligo;
		if (!@temparr2)
		{
			if ((substr $oligo, 0, 3) eq "HUR")
			{
				if ((!$hur) && ($oligo_num == 4))
				{
					push (@temparr, $start);
					push (@temparr, $end);
					$startf = $start;
						 
					$hur = 1;
				}
				elsif ((!$hur) && ($oligo_num == 3))
				{
					push (@temparr, $start);
					push (@temparr, $end);
					$startf = $start;
						
					$hur = 1;
				}
				elsif ((!$hur) && ($oligo_num == 2))
				{
					push (@temparr, $start);
					push (@temparr, $end);
					$startf = $start;
						
					$hur = 1;
				}
				elsif ((!$hur) && ($oligo_num == 1))
				{
					push (@temparr, $start);
					push (@temparr, $end);
					$startf = $start;
						
					$hur = 1;
				}
				elsif ((!$hur) && ($oligo eq "HUR"))
				{
					push (@temparr, $start);
					push (@temparr, $end);
					$startf = $start;
						
					$hur = 1;
				}
			}
			elsif ((!$hu) && ($oligo eq "HU"))
			{
				push (@temparr2, $start);
				push (@temparr2, $end);
				$startf = $start;
					
				$hu = 1;
			}
			elsif ((!$hu) && ($oligo eq "HU2"))
			{
				push (@temparr2, $start);
				push (@temparr2, $end);
				$startf = $start;
					
				$hu = 1;
			}
			elsif ((substr $oligo, 0, 3) eq "HDF")
			{
				if (@temparr)
				{
					if ((!$hdf) && ($oligo_num == 4))
					{
						if ($startf < $start)
						{	
							shift @temparr;
							push (@temparr, $start);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
						elsif ($startf > $start)
						{
							pop @temparr;
							push (@temparr, $end);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
							
						$hdf = 1;
					}
					elsif ((!$hdf) && ($oligo_num == 3))
					{
						if ($startf < $start)
						{	
							shift @temparr;
							push (@temparr, $start);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
						elsif ($startf > $start)
						{
							pop @temparr;
							push (@temparr, $end);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
							
						$hdf = 1;
					}
					elsif ((!$hdf) && ($oligo_num == 2))
					{
						if ($startf < $start)
						{	
							shift @temparr;
							push (@temparr, $start);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
						elsif ($startf > $start)
						{
							pop @temparr;
							push (@temparr, $end);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
							
						$hdf = 1;
					}
					elsif ((!$hdf) && ($oligo_num == 1))
					{
						if ($startf < $start)
						{	
							shift @temparr;
							push (@temparr, $start);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
						elsif ($startf > $start)
						{
							pop @temparr;
							push (@temparr, $end);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
							
						$hdf = 1;
					}
					elsif ((!$hdf) && (!$oligo_num))
					{
						if ($startf < $start)
						{	
							shift @temparr;
							push (@temparr, $start);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
						elsif ($startf > $start)
						{
							pop @temparr;
							push (@temparr, $end);
							push (@temparr, $chr);
							push (@temparr, $score);
						}
							
						$hdf = 1;
					}
				}
				else
				{
					print "Warning: $maid is missing a matching pair.\n";
				}
			}
		}
		else
		{
			if ($oligo eq "HD")
			{
				if ($startf < $start)
				{	
					shift @temparr2;
					push (@temparr2, $start);
					push (@temparr2, $chr);
					push (@temparr2, $score);
				}
				elsif ($startf > $start)
				{
					pop @temparr2;
					push (@temparr2, $end);
					push (@temparr2, $chr);
					push (@temparr2, $score);
				}	
				$hd = 1;
			}
		}
		if ($n == $#lines)
		{
			if (@temparr == 4)
			{
				if ($temparr[0] > $temparr[1])
				{
					$line = $temparr[2]."\t".$temparr[1]."\t".$temparr[0]."\t".$maidf."DEL\t".$temparr[3]."\t-\n";
				}
				else
				{
					$line = $temparr[2]."\t".$temparr[0]."\t".$temparr[1]."\t".$maidf."DEL\t".$temparr[3]."\t+\n";
				}
				$cnt++;
				open (OUTPUT, '>>del.txt');
				print OUTPUT "$line";
				close (OUTPUT);
			}
			elsif (@temparr2 == 4)
			{
				if ($temparr2[0] > $temparr2[1])
				{
					$line = $temparr2[2]."\t".$temparr2[1]."\t".$temparr2[0]."\t".$maidf."DEL\t".$temparr2[3]."\t-\n";
				}
				else
				{
					$line = $temparr2[2]."\t".$temparr2[0]."\t".$temparr2[1]."\t".$maidf."DEL\t".$temparr2[3]."\t+\n";
				}
				$cnt++;
				open (OUTPUT, '>>del.txt');
				print OUTPUT "$line";
				close (OUTPUT);
			}
		}
	} 
	else # 100 -> 101
	{
		$cnt2++;
		if (@temparr == 4)
		{
			if ($temparr[0] > $temparr[1])
			{
				$line = $temparr[2]."\t".$temparr[1]."\t".$temparr[0]."\t".$maidf."DEL\t".$temparr[3]."\t-\n";
			}
			else
			{
				$line = $temparr[2]."\t".$temparr[0]."\t".$temparr[1]."\t".$maidf."DEL\t".$temparr[3]."\t+\n";
			}
			$cnt++;
			open (OUTPUT, '>>del.txt');
			print OUTPUT "$line";
			close (OUTPUT);
		}
		elsif (@temparr2 == 4)
		{
			if ($temparr2[0] > $temparr2[1])
			{
				$line = $temparr2[2]."\t".$temparr2[1]."\t".$temparr2[0]."\t".$maidf."DEL\t".$temparr2[3]."\t-\n";
			}
			else
			{
				$line = $temparr2[2]."\t".$temparr2[0]."\t".$temparr2[1]."\t".$maidf."DEL\t".$temparr2[3]."\t+\n";
			}
			$cnt++;
			open (OUTPUT, '>>del.txt');
			print OUTPUT "$line";
			close (OUTPUT);
		}	
		else
		{
			open (OUTPUT, '>>no_del.txt');
			print OUTPUT $maidf."\n";
			close (OUTPUT);
		}
		@temparr = ();
		@temparr2 = ();
		$maidf = $maid;
		$hur = 0;
		$hdf = 0;
		$hu = 0;
		$hd = 0;
		
		if (($oligo_num) && ((substr $oligo, 0, 3) eq "HUR"))
		{
			if ((!$hur) && ($oligo_num == 4))
			{
				push (@temparr, $start);
				push (@temparr, $end);
				$startf = $start;
				
				$hur = 1;
			}
			if ((!$hur) && ($oligo_num == 3))
			{
				push (@temparr, $start);
				push (@temparr, $end);
				$startf = $start;
				
				$hur = 1;
			}
			if ((!$hur) && ($oligo_num == 2))
			{
				push (@temparr, $start);
				push (@temparr, $end);
				$startf = $start;
				$hur = 1;
			}
			if ((!$hur) && ($oligo_num == 1))
			{
				push (@temparr, $start);
				push (@temparr, $end);
				$startf = $start;
				
				$hur = 1;
			}
		}
		elsif ((substr $oligo, 0, 2) eq "HU")
		{
			if ((!$hur) && ($oligo eq "HUR"))
			{
				push (@temparr, $start);
				push (@temparr, $end);
				$startf = $start;
				
				$hur = 1;
			}
			elsif ((!$hu) && ($oligo eq "HU"))
			{
				push (@temparr2, $start);
				push (@temparr2, $end);
				$startf = $start;
				
				$hu = 1;
			}
			elsif ((!$hu) && ($oligo eq "HU2"))
			{
				push (@temparr2, $start);
				push (@temparr2, $end);
				$startf = $start;		
				$hu = 1;
			}
		}
	}
}

my $filename2 = "del.txt";
open (FILE2, $filename2);
my @lines2 = reverse <FILE2>;
for (my $n = 0; $n < @lines2; $n++)
{
	open (OUTPUT, '>>deletion.txt');
	print OUTPUT "$lines2[$n]";
	close (OUTPUT);
}

print "$cnt2 MAID processed\n";
print "Deletion Count= $cnt\n";

`rm del.txt`;
