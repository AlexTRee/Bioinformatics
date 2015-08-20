#!/usr/bin/perl

########################################################
# Author: Tiange Cui                                   #
# Date: August-07-2015                                 #
# Funtion: Parse BLAT search results, find match pairs.# 
# Arguments: None                                      #
# Execute: perl blat_parser.pl                         #
# Requirements: List::Util, blat_result.txt            #
########################################################

use strict;
use List::Util qw(first);
use List::Util qw[min max];

`rm match.txt`;
`rm maid_list.txt`;

my (%hash, %seen, %count) = ();
my (@no_matching_pair, @chr_list, @matchset1, @matchset2, @matchset3, @bd, @bu, @hd, @hu, @sd, @su) = ();
my ($initlinenum, $matchnum, $index, $max, $toligo1, $toligo2, $toligo3, $cnt, $tempstart, $min_diff);
$initlinenum = $matchnum = $index = $max = 0;
$cnt = 1;
my @list = ();

my $filename = "blat_result.txt";
open (FILE, $filename)|| die "Can't open data file!\n";
my @lines = <FILE>;
close (FILE);

my @firstline = split ("\t", $lines[0]);
my $maidf = $firstline[0];
my $oligof = $firstline[1];
my @oligosf = grep {$_} split /(\d+)/, $oligof;
$oligof = $oligosf[0];
if (defined $oligosf[1])
{
	my $oligof_num = $oligosf[1];
}
my $chrf = $firstline[5];
my $strf = $firstline[6];
my $startf = $firstline[7];
my $endf = $firstline[8];
my $scoref = $firstline[9];
chomp $scoref;
my $keyf = $chrf.".".$strf.".".$startf.".".$endf;
$hash{$keyf} = $lines[0];

push (@matchset1, $chrf."\t".$strf."\t".$startf."\t".$endf."\t".$scoref."\t".$maidf."\t".$oligof."\n");
push (@list, $maidf);

for (my $i = 1; $i < @lines; $i++)
{
	my @arr = split ("\t", $lines[$i]);
	my $maid = $arr[0];
	my $oligo = $arr[1];
	my @oligos = ();
	@oligos = grep {$_} split /(\d+)/, $oligo;
	$oligo = $oligos[0];
	my $oligo_num;
	if (defined $oligos[1])
	{
		$oligo_num = $oligos[1];
	}	
	my $chr = $arr[5];
	chomp $chr;
	my $str = $arr[6];
	
	my $start = $arr[7];
	my $end = $arr[8];
	my $score = $arr[9];
	chomp $score;
	my $key = $chr.".".$str.".".$start.".".$end;
	$hash{$key} = $lines[$i];
	
	if ($i == $#lines)
	{
		if (substr ($oligo, 0, 1) eq "T")
		{
			if ((substr $oligo, -1) eq "F")
			{
				if ($oligo_num)
				{
					push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
				}
				else
				{
					push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
				}
			}
			elsif ((substr $oligo, -1) eq "P")
			{
				if ($oligo_num)
				{
					push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
				}
				else
				{
					push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
				}
			}
			elsif ((substr $oligo, -1) eq "R")
			{
				if ($oligo_num)
				{
					push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
				}
				else
				{
					push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
				}
			}
		}
		else
		{
			if ($oligo_num)
			{
				push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
			}
			else
			{
				push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
			}
		}
				
		if ($matchnum == 0)   # only has one pair HD -> SU
		{
			push (@no_matching_pair, @matchset1);
			if ($oligo_num)
			{
				push (@no_matching_pair, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
			}
			else
			{
				push (@no_matching_pair, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
			}
		}
		elsif ($matchnum >= 1)  # HU  HD -> SU
		{
			if ((substr $oligof, 0, 1) ne "T")
			{
				if ((length($oligof) == 2))
				{
					compare_none_taq_2let_push(\@matchset1, \@matchset2);
				}
				else
				{
					compare_none_taq_3let_print(\@matchset1, \@matchset2);
				}
			}
			elsif ((substr $oligof, 0, 1) eq "T")
			{
				compare_taq(\@matchset1, \@matchset2, \@matchset3);
			}
		}
		
		for(my $c = 0; $c < @chr_list; $c++) 
		{     
			$count{$chr_list[$c]}++;     
			if($count{$chr_list[$c]} > $count{$chr_list[$max]}) 
			{ 
				$max = $c;
			} 
		}

		my %unique = ();
		
		foreach my $item (@no_matching_pair)
		{
			$unique{$item} ++;
		}
		my @uniqarr = keys %unique;
		@uniqarr = sort {(split(/\t/,$a))[6] cmp (split(/\t/,$b))[6]} @uniqarr;
		
		foreach my $p (@uniqarr)
		{
			my @arr = split ("\t", $p);
			if (($arr[0] eq $chr_list[$max]) && ((substr $arr[2], 0, 1) == (substr $tempstart, 0, 1)) && (length($arr[2]) == length($tempstart)))
			{
				if ((substr $arr[6], 0, 2) eq "HD")
				{
					if (((substr $arr[6], 0, 3) eq "HDF") || ((substr $arr[6], 0, 3) eq "HDR"))
					{
						my $distance = ($arr[2] > $tempstart) ? ($arr[2] - $tempstart) : ($tempstart - $arr[2]);
						if ($distance < 100000)
						{
							open (OUTPUT, '>>match.txt');
							print OUTPUT ("$p");
							close (OUTPUT);
						}
					}
					else
					{
						push(@hd, $p);
					}
				}
				elsif ((substr $arr[6], 0, 2) eq "HU")
				{
					if (((substr $arr[6], 0, 3) eq "HUF") || ((substr $arr[6], 0, 3) eq "HUR"))
					{
						my $distance = ($arr[2] > $tempstart) ? ($arr[2] - $tempstart) : ($tempstart - $arr[2]);
						if ($distance < 100000)
						{
							open (OUTPUT, '>>match.txt');
							print OUTPUT ("$p");
							close (OUTPUT);
						}
					}
					else
					{
						push(@hu, $p);
					}					
				}
				elsif ((substr $arr[6], 0, 2) eq "SD")
				{
					push(@sd, $p);
				}
				elsif ((substr $arr[6], 0, 2) eq "SU")
				{
					push(@su, $p);
				}
				elsif ((substr $arr[6], 0, 2) eq "BD")
				{
					push(@bd, $p);
				}
				elsif ((substr $arr[6], 0, 2) eq "BU")
				{
					push(@bu, $p);
				}
				else
				{
					my $distance = ($arr[2] > $tempstart) ? ($arr[2] - $tempstart) : ($tempstart - $arr[2]);
					if ($distance < 100000)
					{
						open (OUTPUT, '>>match.txt');
						print OUTPUT ("$p");
						close (OUTPUT);
					}
				}
			}
		}
		
		my $pickbu;
		$min_diff = 999999999;
		my @set1 = split ("\t", $bu[0]);
		my $bu_oligof = $set1[6];
		chomp $bu_oligof;
		my $bu_startf = $set1[2];
		my $bu_endf = $set1[3];
		my $diff = (($bu_startf > $tempstart) ? ($bu_startf - $tempstart) : ($tempstart - $bu_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$pickbu = $bu[0];
		}
		foreach my $item (@bu)
		{
			my @set1 = split ("\t", $item);
			my $bu_oligo = $set1[6];
			chomp $bu_oligo;
			my $bu_start = $set1[2];
			my $bu_end = $set1[3];
			my $diff = (($bu_start > $tempstart) ? ($bu_start - $tempstart) : ($tempstart - $bu_start));
			if ($bu_oligo eq $bu_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickbu = $item;
					$bu_oligof = $bu_oligo;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$pickbu");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickbu = $item;
					$bu_oligof = $bu_oligo;
				}
				if ($item eq $bu[$#bu])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$pickbu");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$pickbu");
		close (OUTPUT);
		
		my $pickbd;
		$min_diff = 999999999;
		my @set1 = split ("\t", $bd[0]);
		my $bd_oligof = $set1[6];
		chomp $bd_oligof;
		my $bd_startf = $set1[2];
		my $bd_endf = $set1[3];
		my $diff = (($bd_startf > $tempstart) ? ($bd_startf - $tempstart) : ($tempstart - $bd_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$pickbd = $bd[0];
		}
		foreach my $item (@bd)
		{
			my @set1 = split ("\t", $item);
			my $bd_oligo = $set1[6];
			chomp $bd_oligo;
			my $bd_start = $set1[2];
			my $bd_end = $set1[3];
			my $diff = (($bd_start > $tempstart) ? ($bd_start - $tempstart) : ($tempstart - $bd_start));
			if ($bd_oligo eq $bd_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickbd = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$pickbd");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickbd = $item;
					$bd_oligof = $bd_oligo;
				}
				if ($item eq $bd[$#bd])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$pickbd");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$pickbd");
		close (OUTPUT);
		
		my $pickhu;
		$min_diff = 999999999;
		my @set1 = split ("\t", $hu[0]);
		my $hu_oligof = $set1[6];
		chomp $hu_oligof;
		my $hu_startf = $set1[2];
		my $hu_endf = $set1[3];
		my $diff = (($hu_startf > $tempstart) ? ($hu_startf - $tempstart) : ($tempstart - $hu_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$pickhu = $hu[0];
		}
		foreach my $item (@hu)
		{
			my @set1 = split ("\t", $item);
			my $hu_oligo = $set1[6];
			chomp $hu_oligo;
			my $hu_start = $set1[2];
			my $hu_end = $set1[3];
			my $diff = (($hu_start > $tempstart) ? ($hu_start - $tempstart) : ($tempstart - $hu_start));
			if ($hu_oligo eq $hu_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickhu = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$pickhu");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickhu = $item;
					$hu_oligof = $hu_oligo;
				}
				if ($item eq $hu[$#hu])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$pickhu");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$pickhu");
		close (OUTPUT);
		
		my $pickhd;
		$min_diff = 999999999;
		my @set1 = split ("\t", $hd[0]);
		my $hd_oligof = $set1[6];
		chomp $hd_oligof;
		my $hd_startf = $set1[2];
		my $hd_endf = $set1[3];
		my $diff = (($hd_startf > $tempstart) ? ($hd_startf - $tempstart) : ($tempstart - $hd_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$pickhd = $hd[0];
		}
		foreach my $item (@hd)
		{
			my @set1 = split ("\t", $item);
			my $hd_oligo = $set1[6];
			chomp $hd_oligo;
			my $hd_start = $set1[2];
			my $hd_end = $set1[3];
			my $diff = (($hd_start > $tempstart) ? ($hd_start - $tempstart) : ($tempstart - $hd_start));
			if ($hd_oligo eq $hd_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickhd = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$pickhd");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickhd = $item;
					$hd_oligof = $hd_oligo;
				}
				if ($item eq $hd[$#hd])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$pickhd");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$pickhd");
		close (OUTPUT);
		
		my $picksu;
		$min_diff = 999999999;
		my @set1 = split ("\t", $su[0]);
		my $su_oligof = $set1[6];
		chomp $su_oligof;
		my $su_startf = $set1[2];
		my $su_endf = $set1[3];
		my $diff = (($su_startf > $tempstart) ? ($su_startf - $tempstart) : ($tempstart - $su_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$picksu = $su[0];
		}
		foreach my $item (@su)
		{
			my @set1 = split ("\t", $item);
			my $su_oligo = $set1[6];
			chomp $su_oligo;
			my $su_start = $set1[2];
			my $su_end = $set1[3];
			my $diff = (($su_start > $tempstart) ? ($su_start - $tempstart) : ($tempstart - $su_start));
			if ($su_oligo eq $su_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$picksu = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$picksu");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$picksu = $item;
					$su_oligof = $su_oligo;
				}
				if ($item eq $su[$#su])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$picksu");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$picksu");
		close (OUTPUT);
		
		my $picksd;
		$min_diff = 999999999;
		my @set1 = split ("\t", $sd[0]);
		my $sd_oligof = $set1[6];
		chomp $sd_oligof;
		my $sd_startf = $set1[2];
		my $sd_endf = $set1[3];
		my $diff = (($sd_startf > $tempstart) ? ($sd_startf - $tempstart) : ($tempstart - $sd_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$picksd = $sd[0];
		}
		foreach my $item (@sd)
		{
			my @set1 = split ("\t", $item);
			my $sd_oligo = $set1[6];
			chomp $sd_oligo;
			my $sd_start = $set1[2];
			my $sd_end = $set1[3];
			my $diff = (($sd_start > $tempstart) ? ($sd_start - $tempstart) : ($tempstart - $sd_start));
			if ($sd_oligo eq $sd_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$picksd = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$picksd");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$picksd = $item;
					$sd_oligof = $sd_oligo;
				}
				if ($item eq $sd[$#sd])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$picksd");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$picksd");
		close (OUTPUT);
	}
	
	if ($maid eq $maidf) # 100  ->  100
	{		
		if ($oligo eq $oligof)  # HUF  ->  HUF
		{
			if (@matchset3)
			{
				if ($oligo_num)
				{
					push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
				}
				else
				{
					push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
				}
			}
			elsif (@matchset2)
			{
				if ($oligo_num)
				{
					push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
				}
				else
				{
					push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
				}
			}
			else
			{
				if ($oligo_num)
				{
					push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
				}
				else
				{
					push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
				}
			}
			$oligof = $oligo;
			$maidf = $maid;
			$chrf = $chr;
			$strf = $str;
			$startf = $start;
			$endf = $end;
		}
		else  # HUF -> HDR
		{
			%seen = ();
			#################### CASE 1			
			if ((length($oligo) == length($oligof)) && (length($oligo) == 2) && ($scoref == 100))   # XX ->  XX
			{
				if ((substr $oligo, 0, 1) eq (substr $oligof, 0, 1))   #  HU ->  HD
				{
					if ($oligo_num)
					{
						push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
					}
					else
					{
						push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
					}
					$matchnum++;
					$oligof = $oligo;
					$maidf = $maid;
					$chrf = $chr;
					$strf = $str;
					$startf = $start;
					$endf = $end;
				}
				else  #  HD -> SU
				{
					%seen = ();
					if ($matchnum == 0)   # only has one pair HD -> SU
					{
						push (@no_matching_pair, @matchset1);
					}
					elsif ($matchnum >= 1)  # HU  HD -> SU
					{
						compare_none_taq_2let_push(\@matchset1, \@matchset2);
					}
					
					%seen = ();
					$matchnum = 0;
					$initlinenum = $i;
					(@matchset1, @matchset2) = ();
					if ($oligo_num)
					{
						push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
					}
					else
					{
						push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
					}
					push (@no_matching_pair, @matchset1);
					$oligof = $oligo;
					$maidf = $maid;
					$chrf = $chr;
					$strf = $str;
					$startf = $start;
					$endf = $end;
				}
			}
			#################### CASE 2
			elsif ((length($oligo) != length($oligof)) && (length($oligo) == 3) && ($scoref == 100))  # XX -> XXX
			{
				if ($matchnum == 0)   # only has one pair HD -> TDF
				{
					push (@no_matching_pair, @matchset1);
				}
				elsif ($matchnum >= 1)  # HU  HD -> TDF
				{
					compare_none_taq_2let_push(\@matchset1, \@matchset2);
				}
				%seen = ();
				$matchnum = 0;
				$initlinenum = $i;
				(@matchset1, @matchset2) = ();
				if ((substr $oligo, 0, 1) eq "T")
				{
					if ((substr $oligo, -1) eq "F")
					{
						if ($oligo_num)
						{
							push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
						}
						else
						{
							push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
						}
					}
					elsif ((substr $oligo, -1) eq "P")
					{
						if ($oligo_num)
						{
							push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
						}
						else
						{
							push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
						}
					}
					elsif ((substr $oligo, -1) eq "R")
					{
						if ($oligo_num)
						{
							push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
						}
						else
						{
							push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
						}
					}
				}
				else
				{
					if ($oligo_num)
					{
						push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
					}
					else
					{
						push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
					}
				}
				$oligof = $oligo;
				$maidf = $maid;
				$chrf = $chr;
				$strf = $str;
				$startf = $start;
				$endf = $end;
			}
			#################### CASE 3
			elsif ((length($oligo) == length($oligof)) && (length($oligo) == 3) && ($scoref == 100))  # XXX -> XXX
			{
				if (((substr $oligo, 0, 2) eq (substr $oligof, 0, 2)) && ((substr $oligo, -1) ne (substr $oligof, -1)))   #  HUF ->  HUR
				{
					if (substr ($oligo, 0, 1) eq "T")
					{
						if ((substr $oligo, -1) eq "F")
						{	
							if ($oligo_num)
							{
								push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
							}
							else
							{
								push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
							}
						}
						elsif ((substr $oligo, -1) eq "P")
						{
							if ($oligo_num)
							{
								push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
							}
							else
							{
								push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
							}
						}
						elsif ((substr $oligo, -1) eq "R")
						{
							if ($oligo_num)
							{
								push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
							}
							else
							{
								push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
							}
						}
					}
					else
					{
						if ($oligo_num)
						{
							push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
						}
						else
						{
							push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
						}
					}
					$matchnum++;
					$oligof = $oligo;
					$maidf = $maid;
					$chrf = $chr;
					$strf = $str;
					$startf = $start;
					$endf = $end;
					$initlinenum = $i;
				}
				else   #  HUF -> HDR
				{
					if ($matchnum == 0)   # only has one pair HUR -> HDF
					{
						push (@no_matching_pair, @matchset1);						
					}
					
					elsif ($matchnum >= 1)  # HUF  HUR -> HDF
					{
						if ((substr $oligof, 0, 1) ne "T")
						{
							compare_none_taq_3let_print(\@matchset1, \@matchset2);
						}
						elsif ((substr $oligof, 0, 1) eq "T")
						{
							compare_taq(\@matchset1, \@matchset2, \@matchset3);
						}
					}
					%seen = ();
					$matchnum = 0;
					$initlinenum = $i;
					(@matchset1, @matchset2, @matchset3) = ();
					if ((substr $oligo, 0, 1) ne "T")
					{
						if ($oligo_num)
						{
							push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
						}
						else
						{
							push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
						}
					}
					else
					{
						if ((substr $oligo, 0, 1) eq "T")
						{
							if ((substr $oligo, -1) eq "F")
							{
								if ($oligo_num)
								{
									push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
								}
								else
								{
									push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
								}
							}
							elsif ((substr $oligo, -1) eq "P")
							{
								if ($oligo_num)
								{
									push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
								}
								else
								{
									push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
								}
							}
							elsif ((substr $oligo, -1) eq "R")
							{
								if ($oligo_num)
								{
									push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
								}
								else
								{
									push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
								}
							}
						}
						else
						{
							push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
						}
					}
					$oligof = $oligo;
					$maidf = $maid;
					$chrf = $chr;
					$strf = $str;
					$startf = $start;
					$endf = $end;
				}
			}	
			#################### CASE 4
			elsif ((length($oligo) != length($oligof)) && (length($oligo) == 2) && ($scoref == 100))  # XXX -> XX
			{
				if ($matchnum == 0)   # only has on pair HUR -> SD
				{
					push (@no_matching_pair, @matchset1);
				}
				elsif ($matchnum >= 1)  # HUF  HDF -> SD
				{
					if ((substr $oligof, 0, 1) ne "T")
					{
						compare_none_taq_3let_print(\@matchset1, \@matchset2);
					}
					elsif ((substr $oligof, 0, 1) eq "T")
					{
						compare_taq(\@matchset1, \@matchset2, \@matchset3);
					}
				}
				
				%seen = ();
				$matchnum = 0;
				$initlinenum = $i;
				(@matchset1, @matchset2, @matchset3) = ();
				
				if ($oligo_num)
				{
					push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
				}
				else
				{
					push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
				}
				$oligof = $oligo;
				$maidf = $maid;
				$chrf = $chr;
				$strf = $str;
				$startf = $start;
				$endf = $end;
			}
		}
	}
	else # 100 -> 101****************************************************************
	{
		$cnt++;
		push (@list, $maid);
		if (($matchnum == 0) && (length($oligof) == 2))  # only has one pair HD -> SU
		{
			push (@no_matching_pair, @matchset1);
		}
		elsif ($matchnum >= 1)  # HU  HD -> SU
		{
			if ((substr $oligof, 0, 1) ne "T")
			{
				if ((length($oligof) == 2))
				{
					compare_none_taq_2let_push(\@matchset1, \@matchset2);
				}
				else
				{
					compare_none_taq_3let_print(\@matchset1, \@matchset2);
				}
			}
			elsif ((substr $oligof, 0, 1) eq "T")
			{
				compare_taq(\@matchset1, \@matchset2, \@matchset3);
			}
		}

		%seen = ();
		(@matchset1, @matchset2, @matchset3) = ();
		$matchnum = 0;
		if ((substr $oligo, 0, 1) ne "T")
		{
			if ($oligo_num)
			{
				push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
			}
			else
			{
				push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
			}
		}
		else
		{
			$initlinenum = $i;
			if ((substr $oligo, 0, 1) eq "T")
			{
				if ((substr $oligo, -1) eq "F")
				{
					if ($oligo_num)
					{
						push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
					}
					else
					{
						push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
					}
				}
				elsif ((substr $oligo, -1) eq "P")
				{
					if ($oligo_num)
					{
						push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
					}
					else
					{
						push (@matchset2, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
					}
				}
				elsif ((substr $oligo, -1) eq "R")
				{
					if ($oligo_num)
					{
						push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
					}
					else
					{
						push (@matchset3, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
					}
				}
			}
			else
			{
				if ($oligo_num)
				{
					push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo.$oligo_num."\n");
				}
				else
				{
					push (@matchset1, $chr."\t".$str."\t".$start."\t".$end."\t".$score."\t".$maid."\t".$oligo."\n");
				}
			}
		}
		$oligof = $oligo;
		$maidf = $maid;
		$chrf = $chr;
		$strf = $str;
		$startf = $start;
		$endf = $end;
		
		for(my $c = 0; $c < @chr_list; $c++) 
		{     
			$count{$chr_list[$c]}++;     
			if($count{$chr_list[$c]} > $count{$chr_list[$max]}) 
			{ 
				$max = $c;
			} 
		}
		my %unique = ();
		
		foreach my $item (@no_matching_pair)
		{
			$unique{$item} ++;
		}
		my @uniqarr = keys %unique;
		
		@uniqarr = sort {(split(/\t/,$a))[6] cmp (split(/\t/,$b))[6]} @uniqarr;
		
		foreach my $p (@uniqarr)
		{
			my @arr = split ("\t", $p);
			if (($arr[0] eq $chr_list[$max]) && ((substr $arr[2], 0, 1) == (substr $tempstart, 0, 1)) && (length($arr[2]) == length($tempstart)))
			{
				if ((substr $arr[6], 0, 2) eq "HD")
				{
					if (((substr $arr[6], 0, 3) eq "HDF") || ((substr $arr[6], 0, 3) eq "HDR"))
					{
						my $distance = ($arr[2] > $tempstart) ? ($arr[2] - $tempstart) : ($tempstart - $arr[2]);
						if ($distance < 100000)
						{
							open (OUTPUT, '>>match.txt');
							print OUTPUT ("$p");
							close (OUTPUT);
						}
					}
					else
					{
						push(@hd, $p);
					}
				}
				elsif ((substr $arr[6], 0, 2) eq "HU")
				{
					if (((substr $arr[6], 0, 3) eq "HUF") || ((substr $arr[6], 0, 3) eq "HUR"))
					{
						my $distance = ($arr[2] > $tempstart) ? ($arr[2] - $tempstart) : ($tempstart - $arr[2]);
						if ($distance < 100000)
						{
							open (OUTPUT, '>>match.txt');
							print OUTPUT ("$p");
							close (OUTPUT);
						}
					}
					else
					{
						push(@hu, $p);
					}					
				}
				elsif ((substr $arr[6], 0, 2) eq "SD")
				{
					push(@sd, $p);
				}
				elsif ((substr $arr[6], 0, 2) eq "SU")
				{
					push(@su, $p);
				}
				elsif ((substr $arr[6], 0, 2) eq "BD")
				{
					push(@bd, $p);
				}
				elsif ((substr $arr[6], 0, 2) eq "BU")
				{
					push(@bu, $p);
				}
				else
				{
					my $distance = ($arr[2] > $tempstart) ? ($arr[2] - $tempstart) : ($tempstart - $arr[2]);
					if ($distance < 100000)
					{
						open (OUTPUT, '>>match.txt');
						print OUTPUT ("$p");
						close (OUTPUT);
					}
				}
			}
		}
		
		my $pickbu;
		$min_diff = 999999999;
		my @set1 = split ("\t", $bu[0]);
		my $bu_oligof = $set1[6];
		chomp $bu_oligof;
		my $bu_startf = $set1[2];
		my $bu_endf = $set1[3];
		my $diff = (($bu_startf > $tempstart) ? ($bu_startf - $tempstart) : ($tempstart - $bu_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$pickbu = $bu[0];
		}
		foreach my $item (@bu)
		{
			my @set1 = split ("\t", $item);
			my $bu_oligo = $set1[6];
			chomp $bu_oligo;
			my $bu_start = $set1[2];
			my $bu_end = $set1[3];
			my $diff = (($bu_start > $tempstart) ? ($bu_start - $tempstart) : ($tempstart - $bu_start));
			if ($bu_oligo eq $bu_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickbu = $item;
					$bu_oligof = $bu_oligo;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$pickbu");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickbu = $item;
					$bu_oligof = $bu_oligo;
				}
				if ($item eq $bu[$#bu])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$pickbu");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$pickbu");
		close (OUTPUT);
		
		my $pickbd;
		$min_diff = 999999999;
		my @set1 = split ("\t", $bd[0]);
		my $bd_oligof = $set1[6];
		chomp $bd_oligof;
		my $bd_startf = $set1[2];
		my $bd_endf = $set1[3];
		my $diff = (($bd_startf > $tempstart) ? ($bd_startf - $tempstart) : ($tempstart - $bd_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$pickbd = $bd[0];
		}
		foreach my $item (@bd)
		{
			my @set1 = split ("\t", $item);
			my $bd_oligo = $set1[6];
			chomp $bd_oligo;
			my $bd_start = $set1[2];
			my $bd_end = $set1[3];
			my $diff = (($bd_start > $tempstart) ? ($bd_start - $tempstart) : ($tempstart - $bd_start));
			if ($bd_oligo eq $bd_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickbd = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$pickbd");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickbd = $item;
					$bd_oligof = $bd_oligo;
				}
				if ($item eq $bd[$#bd])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$pickbd");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$pickbd");
		close (OUTPUT);
		
		my $pickhu;
		$min_diff = 999999999;
		my @set1 = split ("\t", $hu[0]);
		my $hu_oligof = $set1[6];
		chomp $hu_oligof;
		my $hu_startf = $set1[2];
		my $hu_endf = $set1[3];
		my $diff = (($hu_startf > $tempstart) ? ($hu_startf - $tempstart) : ($tempstart - $hu_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$pickhu = $hu[0];
		}
		foreach my $item (@hu)
		{
			my @set1 = split ("\t", $item);
			my $hu_oligo = $set1[6];
			chomp $hu_oligo;
			my $hu_start = $set1[2];
			my $hu_end = $set1[3];
			my $diff = (($hu_start > $tempstart) ? ($hu_start - $tempstart) : ($tempstart - $hu_start));
			if ($hu_oligo eq $hu_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickhu = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$pickhu");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickhu = $item;
					$hu_oligof = $hu_oligo;
				}
				if ($item eq $hu[$#hu])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$pickhu");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$pickhu");
		close (OUTPUT);
		
		my $pickhd;
		$min_diff = 999999999;
		my @set1 = split ("\t", $hd[0]);
		my $hd_oligof = $set1[6];
		chomp $hd_oligof;
		my $hd_startf = $set1[2];
		my $hd_endf = $set1[3];
		my $diff = (($hd_startf > $tempstart) ? ($hd_startf - $tempstart) : ($tempstart - $hd_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$pickhd = $hd[0];
		}
		foreach my $item (@hd)
		{
			my @set1 = split ("\t", $item);
			my $hd_oligo = $set1[6];
			chomp $hd_oligo;
			my $hd_start = $set1[2];
			my $hd_end = $set1[3];
			my $diff = (($hd_start > $tempstart) ? ($hd_start - $tempstart) : ($tempstart - $hd_start));
			if ($hd_oligo eq $hd_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickhd = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$pickhd");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$pickhd = $item;
					$hd_oligof = $hd_oligo;
				}
				if ($item eq $hd[$#hd])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$pickhd");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$pickhd");
		close (OUTPUT);
		
		my $picksu;
		$min_diff = 999999999;
		my @set1 = split ("\t", $su[0]);
		my $su_oligof = $set1[6];
		chomp $su_oligof;
		my $su_startf = $set1[2];
		my $su_endf = $set1[3];
		my $diff = (($su_startf > $tempstart) ? ($su_startf - $tempstart) : ($tempstart - $su_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$picksu = $su[0];
		}
		foreach my $item (@su)
		{
			my @set1 = split ("\t", $item);
			my $su_oligo = $set1[6];
			chomp $su_oligo;
			my $su_start = $set1[2];
			my $su_end = $set1[3];
			my $diff = (($su_start > $tempstart) ? ($su_start - $tempstart) : ($tempstart - $su_start));
			if ($su_oligo eq $su_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$picksu = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$picksu");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$picksu = $item;
					$su_oligof = $su_oligo;
				}
				if ($item eq $su[$#su])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$picksu");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$picksu");
		close (OUTPUT);
		
		my $picksd;
		$min_diff = 999999999;
		my @set1 = split ("\t", $sd[0]);
		my $sd_oligof = $set1[6];
		chomp $sd_oligof;
		my $sd_startf = $set1[2];
		my $sd_endf = $set1[3];
		my $diff = (($sd_startf > $tempstart) ? ($sd_startf - $tempstart) : ($tempstart - $sd_startf));
		if ($diff < $min_diff) 
		{
				$min_diff = $diff;
				$picksd = $sd[0];
		}
		foreach my $item (@sd)
		{
			my @set1 = split ("\t", $item);
			my $sd_oligo = $set1[6];
			chomp $sd_oligo;
			my $sd_start = $set1[2];
			my $sd_end = $set1[3];
			my $diff = (($sd_start > $tempstart) ? ($sd_start - $tempstart) : ($tempstart - $sd_start));
			if ($sd_oligo eq $sd_oligof)
			{
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$picksd = $item;
				}
			}
			else
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$picksd");
				close (OUTPUT);
				$min_diff = 999999999;
				if ($diff <= $min_diff) 
				{
					$min_diff = $diff;
					$picksd = $item;
					$sd_oligof = $sd_oligo;
				}
				if ($item eq $sd[$#sd])
				{
					open (OUTPUT, '>>match.txt');
					print OUTPUT ("$picksd");
					close (OUTPUT);
				}
			}
		}
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$picksd");
		close (OUTPUT);
		
		%count = ();
		(@uniqarr, @no_matching_pair, @chr_list, @bd, @bu, @hd, @hu, @sd, @su) = ();
		$max = 0;
		$tempstart = 0;
	}
}

#None-TaqMan, HU, HD, SU, SD, not output results
sub compare_none_taq_2let_push
{
	my (@index1, @index2, @index3, @index4, @matchset1, @matchset2, @match_maid_chr_1, @match_maid_chr_2, @match_chr_1, @match_chr_2, @temp_index1, @temp_index2, @temp_index3, @temp_index4) = ();
	my $matchset1ref = shift;
	my $matchset2ref = shift;
	@matchset1 = @{$matchset1ref};
	@matchset2 = @{$matchset2ref};
	
	push (@no_matching_pair, @matchset1);
	push (@no_matching_pair, @matchset2);
}

#None-TaqMan, HUF, HUR, HDF, HDR, output results, the distance between a pair (eg. HDF & HDR) < 600bp)
sub compare_none_taq_3let_print
{
	my (@index1, @index2, @index3, @index4, @matchset1, @matchset2, @match_maid_chr_1, @match_maid_chr_2, @match_chr_1, @match_chr_2, @temp_index1, @temp_index2, @temp_index3, @temp_index4) = ();
	my $matchset1ref = shift;
	my $matchset2ref = shift;
	@matchset1 = @{$matchset1ref};
	@matchset2 = @{$matchset2ref};
	
	if (scalar(@matchset1) >= scalar(@matchset2))
	{
		foreach my $item (@matchset2)
		{
			my @set2 = split ("\t", $item);
			my $m = $set2[0];
			my $s = $set2[1];
			push (@match_maid_chr_2, $m.$s);
			push (@match_chr_2, $m);
			$seen{$m} = 1;
			$seen{$m.$s} = 1;
		}

		for (my $n = 0; $n < @matchset1; $n++)
		{
			my @set1 = split ("\t", $matchset1[$n]);
			my $m = $set1[0];
			my $s = $set1[1];
			if ($s eq "+") 
			{
				$s = "-";
			}
			elsif ($s eq "-") 
			{
				$s = "+";
			}
			push (@match_maid_chr_1, $m.$s);
			push (@match_chr_1, $m);

			if ($seen{$m.$s}) # 5 +  ->  5 -	
			{
				@index1= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
				@index2= grep {$match_maid_chr_2[$_] eq $m.$s} 0..$#match_maid_chr_2;
				
				if ((@index1 == 1) && (@index2 > 1)) # more than one match
				{
					my @a1 = split ("\t", $matchset1[$index1[0]]);
					$toligo1 = $a1[6];
					chomp $toligo1;
					my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
					$toligo1 = $toligos1[0];
					if (defined $toligos1[1])
					{
						my $toligo1_num = $toligos1[1];
					}
					my $tstart1 = $a1[2];
					my $tend1 = $a1[3];
					my @sorted_index2 = sort {$a <=> $b} @index2;
					foreach my $k (@sorted_index2)
					{	
						my @a2 = split ("\t", $matchset2[$k]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my $seqlen1 = abs($tend2 - $tstart1);
						my $seqlen2 = abs($tend1 - $tstart2);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
						if ($seqlen > 600)
						{
							@index2 = grep {!/$k/} @index2;
						}
						else
						{
							if (! grep(/$k/, @temp_index2)) 
							{
								push @temp_index2, $k;
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index1[0]/, @temp_index1))
							{
								push @temp_index1, $index1[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
					}
					if (@index2 == 0)
					{
						@index1 = grep {!/$index1[0]/} @index1;
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1 > 1) && (@index2 == 1))
				{
					my @a2 = split ("\t", $matchset2[$index2[0]]);
					$toligo2 = $a2[6];
					chomp $toligo2;
					my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
					$toligo2 = $toligos2[0];
					if (defined $toligos2[1])
					{
						my $toligo2_num = $toligos2[1];
					}
					my $tstart2 = $a2[2];
					my $tend2 = $a2[3];
					my @sorted_index1 = sort {$a <=> $b} @index1;
					foreach my $k (@sorted_index1)
					{
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my $seqlen1 = abs($tend1 - $tstart2);
						my $seqlen2 = abs($tend2 - $tstart1);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 600)
						{
							@index1 = grep {!/$k/} @index1;
						}
						else
						{
							if (! grep(/$k/, @temp_index1)) 
							{
								push @temp_index1, $k;
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index2[0]/, @temp_index2))
							{
								push @temp_index2, $index2[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
					}
					if (@index1 == 0)
					{
						@index2 = grep {!/$index2[0]/} @index2;
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1== 1) && (@index2 == 1))
				{
					my @a1 = split ("\t", $matchset1[$index1[0]]);
					$toligo1 = $a1[6];
					chomp $toligo1;
					my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
					$toligo1 = $toligos1[0];
					if (defined $toligos1[1])
					{
						my $toligo1_num = $toligos1[1];
					}
					my $tstart1 = $a1[2];
					my $tend1 = $a1[3];

					my @a2 = split ("\t", $matchset2[$index2[0]]);
					$toligo2 = $a2[6];
					chomp $toligo2;
					my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
					$toligo2 = $toligos2[0];
					if (defined $toligos2[1])
					{
						my $toligo2_num = $toligos2[1];
					}
					my $tstart2 = $a2[2];
					my $tend2 = $a2[3];

					my $seqlen1 = abs($tend2 - $tstart1);
					my $seqlen2 = abs($tend1 - $tstart2);
					my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
					if ($seqlen > 600)
					{
						shift @index1;
						shift @index2;
					}
					else
					{
						if (! grep(/$index1[0]/, @temp_index1))
						{
							push @temp_index1, $index1[0];
							if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
						}
						if (! grep(/$index2[0]/, @temp_index2))
						{
							push @temp_index2, $index2[0];
							if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
						}
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1> 1) && (@index2 > 1))
				{
					my @tindex1 = @index1;
					my @tindex2 = @index2;
					my @sorted_index1 = sort {$a <=> $b} @index1;
					foreach my $k (@sorted_index1)
					{
						@tindex2 = @index2;
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index2 = sort {$a <=> $b} @index2;
						foreach my $j (@sorted_index2)
						{
							my @a2 = split ("\t", $matchset2[$j]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 600)
							{
								@tindex2 = grep {!/$j/} @tindex2;
							}
							else
							{
								if (! grep(/$j/, @temp_index2)) 
								{
									push @temp_index2, $j;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$k/, @temp_index1)) 
								{
									push @temp_index1, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (!@tindex2)
						{
							@tindex1 = grep {!/$k/} @tindex1;
						}								
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}
			
				if ($n == $#matchset1)
				{
					if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
					{
						push (@no_matching_pair, @matchset1);
						push (@no_matching_pair, @matchset2);
					}
				}
			}
			else # 5 -  -> 5 -
			{
				if ($seen{$m})
				{
					@index3= grep {$match_chr_1[$_] eq $m} 0..$#match_chr_1;
					if ($s eq "+") {$s = "-";}
					elsif ($s eq "-") {$s = "+";}
					@index4= grep {$match_chr_2[$_] eq $m} 0..$#match_chr_2;
					my @sorted_index3 = sort {$a <=> $b} @index3;
					foreach my $k (@sorted_index3)
					{
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index4 = sort {$a <=> $b} @index4;
						foreach my $j (@sorted_index4)
						{
							my @a2 = split ("\t", $matchset2[$j]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 600)
							{
								@index4 = grep {!/$j/} @index4;
							}
							else
							{
								if (! grep(/$j/, @temp_index4)) 
								{
									push @temp_index4, $j;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$k/, @temp_index3)) 
								{
									push @temp_index3, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index4 == 0)
						{
							@index3 = grep {!/$k/} @index3;
						}
					}
					@index3 = @temp_index3;
					@index4 = @temp_index4;
				}
				elsif ($n == $#matchset1)
				{
					if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
					{
						push (@no_matching_pair, @matchset1);
						push (@no_matching_pair, @matchset2);
					}
				}
			}
		}
	}
	else
	{
		foreach my $item (@matchset1)
		{
			my @set1 = split ("\t", $item);
			my $m = $set1[0];
			my $s = $set1[1];
			push (@match_maid_chr_1, $m.$s);
			push (@match_chr_1, $m);
			$seen{$m} = 1;
			$seen{$m.$s} = 1;
		}
		for (my $n = 0; $n < @matchset2; $n++)
		{
			my @set2 = split ("\t", $matchset2[$n]);
			my $m = $set2[0];
			my $s = $set2[1];
			if ($s eq "+") {$s = "-";}
			elsif ($s eq "-") {$s = "+";}
			push (@match_maid_chr_2, $m.$s);
			push (@match_chr_2, $m);
			if ($seen{$m.$s}) # 5 +  ->  5 -
			{
				@index1= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
				@index2= grep {$match_maid_chr_2[$_] eq $m.$s} 0..$#match_maid_chr_2;
				
				if ((@index1 == 1) && (@index2 > 1)) # more than one match
				{
					my @a1 = split ("\t", $matchset1[$index1[0]]);
					$toligo1 = $a1[6];
					chomp $toligo1;
					my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
					$toligo1 = $toligos1[0];
					if (defined $toligos1[1])
					{
						my $toligo1_num = $toligos1[1];
					}
					my $tstart1 = $a1[2];
					my $tend1 = $a1[3];
					my @sorted_index2 = sort {$a <=> $b} @index2;
					foreach my $k (@sorted_index2)
					{	
						my @a2 = split ("\t", $matchset2[$k]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my $seqlen1 = abs($tend2 - $tstart1);
						my $seqlen2 = abs($tend1 - $tstart2);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
						if ($seqlen > 600)
						{
							@index2 = grep {!/$k/} @index2;
						}
						else
						{
							if (! grep(/$k/, @temp_index2)) 
							{
								push @temp_index2, $k;
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index1[0]/, @temp_index1))
							{
								push @temp_index1, $index1[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
					}
					if (@index2 == 0)
					{
						@index1 = grep {!/$index1[0]/} @index1;
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1 > 1) && (@index2 == 1))
				{
					my @a2 = split ("\t", $matchset2[$index2[0]]);
					$toligo2 = $a2[6];
					chomp $toligo2;
					my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
					$toligo2 = $toligos2[0];
					if (defined $toligos2[1])
					{
						my $toligo2_num = $toligos2[1];
					}
					my $tstart2 = $a2[2];
					my $tend2 = $a2[3];
					my @sorted_index1 = sort {$a <=> $b} @index1;
					foreach my $k (@sorted_index1)
					{
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my $seqlen1 = abs($tend1 - $tstart2);
						my $seqlen2 = abs($tend2 - $tstart1);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 600)
						{
							@index1 = grep {!/$k/} @index1;
						}
						else
						{
							if (! grep(/$k/, @temp_index1)) 
							{
								push @temp_index1, $k;
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index2[0]/, @temp_index2))
							{
								push @temp_index2, $index2[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
					}
					if (@index1 == 0)
					{
						@index2 = grep {!/$index2[0]/} @index2;
					}	
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1== 1) && (@index2 == 1))
				{
					my @a1 = split ("\t", $matchset1[$index1[0]]);
					$toligo1 = $a1[6];
					chomp $toligo1;
					my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
					$toligo1 = $toligos1[0];
					if (defined $toligos1[1])
					{
						my $toligo1_num = $toligos1[1];
					}
					my $tstart1 = $a1[2];
					my $tend1 = $a1[3];

					my @a2 = split ("\t", $matchset2[$index2[0]]);
					$toligo2 = $a2[6];
					chomp $toligo2;
					my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
					$toligo2 = $toligos2[0];
					if (defined $toligos2[1])
					{
						my $toligo2_num = $toligos2[1];
					}
					my $tstart2 = $a2[2];
					my $tend2 = $a2[3];

					my $seqlen1 = abs($tend2 - $tstart1);
					my $seqlen2 = abs($tend1 - $tstart2);
					my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
					if ($seqlen > 600)
					{
						shift @index1;
						shift @index2;
					}
					else
					{
						if (! grep(/$index1[0]/, @temp_index1))
						{
							push @temp_index1, $index1[0];
							if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
				
						}
						if (! grep(/$index2[0]/, @temp_index2))
						{
							push @temp_index2, $index2[0];
							if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
						}
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1> 1) && (@index2 > 1))
				{
					my @tindex1 = @index1;
					my @tindex2 = @index2;
					my @sorted_index1 = sort {$a <=> $b} @index1;
					foreach my $k (@sorted_index1)
					{
						@tindex2 = @index2;
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index2 = sort {$a <=> $b} @index2;
						foreach my $j (@sorted_index2)
						{
							my @a2 = split ("\t", $matchset2[$j]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 600)
							{
								@tindex2 = grep {!/$j/} @tindex2;
							}
							else
							{
								if (! grep(/$j/, @temp_index2)) 
								{
									push @temp_index2, $j;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$k/, @temp_index1)) 
								{
									push @temp_index1, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (!@tindex2)
						{
							@tindex1 = grep {!/$k/} @tindex1;
						}										
					}
					 @index1 = @temp_index1;
					 @index2 = @temp_index2;
				}
			
				if ($n == $#matchset2)
				{
					if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
					{
						push (@no_matching_pair, @matchset1);
						push (@no_matching_pair, @matchset2);
					}
				}
			}
			else # 5 -  -> 5 -
			{
				if ($seen{$m})
				{
					@index3= grep {$match_chr_1[$_] eq $m} 0..$#match_chr_1;
					if ($s eq "+") {$s = "-";}
					elsif ($s eq "-") {$s = "+";}
					@index4= grep {$match_chr_2[$_] eq $m} 0..$#match_chr_2;
					my @sorted_index3 = sort {$a <=> $b} @index3;
					foreach my $k (@sorted_index3)
					{
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index4 = sort {$a <=> $b} @index4;
						foreach my $j (@sorted_index4)
						{
							my @a2 = split ("\t", $matchset2[$j]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 600)
							{
								@index4 = grep {!/$j/} @index4;
							}
							else
							{
								if (! grep(/$j/, @temp_index4)) 
								{
									push @temp_index4, $j;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$k/, @temp_index3)) 
								{
									push @temp_index3, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index4 == 0)
						{
							@index3 = grep {!/$k/} @index3;
						}
					}
					@index3 = @temp_index3;
					@index4 = @temp_index4;
				
					if ($n == $#matchset2)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset2);
						}
					}
				}
				elsif ($n == $#matchset2)
				{
					if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
					{
						push (@no_matching_pair, @matchset1);
						push (@no_matching_pair, @matchset2);
					}
				}
			}
		}
	}
	
	@index1 = @temp_index1;
	@index2 = @temp_index2;
	@index3 = @temp_index3;
	@index4 = @temp_index4;
	
	if (@temp_index1 != 0)
	{
		my @temp = split ("\t", $matchset1[$temp_index1[$#temp_index1]]);
		$tempstart = $temp[2];
	}
	elsif (@temp_index2 != 0)
	{
		my @temp = split ("\t", $matchset2[$temp_index2[$#temp_index2]]);
		$tempstart = $temp[2];
	}
	
	if ((@temp_index1 > 1) || (@temp_index2 > 1))
	{
		foreach my $t (@temp_index1)
		{
			push (@no_matching_pair, $matchset1[$t]);
		}
		foreach my $v (@temp_index2)
		{
			push (@no_matching_pair, $matchset2[$v]);
		}
		(@index1, @index2) = ();
	}
	
	if (@temp_index3 != 0)
	{
		my @temp = split ("\t", $matchset1[$temp_index3[$#temp_index3]]);
		$tempstart = $temp[2];
	}
	elsif (@temp_index4 != 0)
	{
		my @temp = split ("\t", $matchset2[$temp_index4[$#temp_index4]]);
		$tempstart = $temp[2];
	}
	
	if ((@temp_index3 > 1) || (@temp_index4 > 1))
	{
		foreach my $t (@temp_index3)
		{
			push (@no_matching_pair, $matchset1[$t]);
		}
		foreach my $v (@temp_index4)
		
		{
			
			push (@no_matching_pair, $matchset2[$v]);
		}
		
		(@index3, @index4) = ();
	}
	
	foreach my $j (@index1)
	{
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$matchset1[$j]");
	}
	foreach my $p (@index2)
	{
		open (OUTPUT, '>>match.txt');
		print OUTPUT ("$matchset2[$p]");
		close (OUTPUT);
	}
		
	if (!@index1)
	{
		foreach my $k (@index3)
		{
			open (OUTPUT2, '>>match.txt');
			print OUTPUT2 ("$matchset1[$k]");
		}
	}
	if (!@index2)
	{
		foreach my $n (@index4)
		{
			open (OUTPUT2, '>>match.txt');
			print OUTPUT2 ("$matchset2[$n]");
			close (OUTPUT2);
		}
	}
}

#None-TaqMan, HUF, HUR, HDF, HDR, not output results, the distance between a pair (eg. HDF & HDR) < 600bp)
sub compare_none_taq_3let_push
{
	my (@index1, @index2, @index3, @index4, @matchset1, @matchset2, @match_maid_chr_1, @match_maid_chr_2, @match_chr_1, @match_chr_2, @temp_index1, @temp_index2, @temp_index3, @temp_index4) = ();
	my $matchset1ref = shift;
	my $matchset2ref = shift;
	@matchset1 = @{$matchset1ref};
	@matchset2 = @{$matchset2ref};
	
	if (scalar(@matchset1) >= scalar(@matchset2))
	{
		foreach my $item (@matchset2)
		{
			my @set2 = split ("\t", $item);
			my $m = $set2[0];
			my $s = $set2[1];
			push (@match_maid_chr_2, $m.$s);
			push (@match_chr_2, $m);
			$seen{$m} = 1;
			$seen{$m.$s} = 1;
		}

		for (my $n = 0; $n < @matchset1; $n++)
		{
			my @set1 = split ("\t", $matchset1[$n]);
			my $m = $set1[0];
			my $s = $set1[1];
			if ($s eq "+") 
			{
				$s = "-";
			}
			elsif ($s eq "-") 
			{
				$s = "+";
			}
			push (@match_maid_chr_1, $m.$s);
			push (@match_chr_1, $m);

			if ($seen{$m.$s}) # 5 +  ->  5 -	
			{
				@index1= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
				@index2= grep {$match_maid_chr_2[$_] eq $m.$s} 0..$#match_maid_chr_2;
				
				if ((@index1 == 1) && (@index2 > 1)) # more than one match
				{
					my @a1 = split ("\t", $matchset1[$index1[0]]);
					$toligo1 = $a1[6];
					chomp $toligo1;
					my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
					$toligo1 = $toligos1[0];
					if (defined $toligos1[1])
					{
						my $toligo1_num = $toligos1[1];
					}
					my $tstart1 = $a1[2];
					my $tend1 = $a1[3];
					my @sorted_index2 = sort {$a <=> $b} @index2;
					foreach my $k (@sorted_index2)
					{	
						my @a2 = split ("\t", $matchset2[$k]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my $seqlen1 = abs($tend2 - $tstart1);
						my $seqlen2 = abs($tend1 - $tstart2);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
						if ($seqlen > 600)
						{
							@index2 = grep {!/$k/} @index2;
						}
						else
						{
							if (! grep(/$k/, @temp_index2)) 
							{
								push @temp_index2, $k;
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index1[0]/, @temp_index1))
							{
								push @temp_index1, $index1[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
					}
					if (@index2 == 0)
					{
						@index1 = grep {!/$index1[0]/} @index1;
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1 > 1) && (@index2 == 1))
				{
					my @a2 = split ("\t", $matchset2[$index2[0]]);
					$toligo2 = $a2[6];
					chomp $toligo2;
					my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
					$toligo2 = $toligos2[0];
					if (defined $toligos2[1])
					{
						my $toligo2_num = $toligos2[1];
					}
					my $tstart2 = $a2[2];
					my $tend2 = $a2[3];
					my @sorted_index1 = sort {$a <=> $b} @index1;
					foreach my $k (@sorted_index1)
					{
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my $seqlen1 = abs($tend1 - $tstart2);
						my $seqlen2 = abs($tend2 - $tstart1);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 600)
						{
							@index1 = grep {!/$k/} @index1;
						}
						else
						{
							if (! grep(/$k/, @temp_index1)) 
							{
								push @temp_index1, $k;
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index2[0]/, @temp_index2))
							{
								push @temp_index2, $index2[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
					}
					if (@index1 == 0)
					{
						@index2 = grep {!/$index2[0]/} @index2;
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1== 1) && (@index2 == 1))
				{
					my @a1 = split ("\t", $matchset1[$index1[0]]);
					$toligo1 = $a1[6];
					chomp $toligo1;
					my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
					$toligo1 = $toligos1[0];
					if (defined $toligos1[1])
					{
						my $toligo1_num = $toligos1[1];
					}
					my $tstart1 = $a1[2];
					my $tend1 = $a1[3];

					my @a2 = split ("\t", $matchset2[$index2[0]]);
					$toligo2 = $a2[6];
					chomp $toligo2;
					my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
					$toligo2 = $toligos2[0];
					if (defined $toligos2[1])
					{
						my $toligo2_num = $toligos2[1];
					}
					my $tstart2 = $a2[2];
					my $tend2 = $a2[3];

					my $seqlen1 = abs($tend2 - $tstart1);
					my $seqlen2 = abs($tend1 - $tstart2);
					my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
					if ($seqlen > 600)
					{
						shift @index1;
						shift @index2;
					}
					else
					{
						if (! grep(/$index1[0]/, @temp_index1))
						{
							push @temp_index1, $index1[0];
							if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
						}
						if (! grep(/$index2[0]/, @temp_index2))
						{
							push @temp_index2, $index2[0];
							if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
						}
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1> 1) && (@index2 > 1))
				{
					my @tindex1 = @index1;
					my @tindex2 = @index2;
					my @sorted_index1 = sort {$a <=> $b} @index1;
					foreach my $k (@sorted_index1)
					{
						@tindex2 = @index2;
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index2 = sort {$a <=> $b} @index2;
						foreach my $j (@sorted_index2)
						{
							my @a2 = split ("\t", $matchset2[$j]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 600)
							{
								@tindex2 = grep {!/$j/} @tindex2;
							}
							else
							{
								if (! grep(/$j/, @temp_index2)) 
								{
									push @temp_index2, $j;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$k/, @temp_index1)) 
								{
									push @temp_index1, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (!@tindex2)
						{
							@tindex1 = grep {!/$k/} @tindex1;
						}								
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}
			
				if ($n == $#matchset1)
				{
					if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
					{
						push (@no_matching_pair, @matchset1);
						push (@no_matching_pair, @matchset2);
					}
				}
			}
			else # 5 -  -> 5 -
			{
				if ($seen{$m})
				{
					@index3= grep {$match_chr_1[$_] eq $m} 0..$#match_chr_1;
					if ($s eq "+") {$s = "-";}
					elsif ($s eq "-") {$s = "+";}
					@index4= grep {$match_chr_2[$_] eq $m} 0..$#match_chr_2;
					my @sorted_index3 = sort {$a <=> $b} @index3;
					foreach my $k (@sorted_index3)
					{
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index4 = sort {$a <=> $b} @index4;
						foreach my $j (@sorted_index4)
						{
							my @a2 = split ("\t", $matchset2[$j]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 600)
							{
								@index4 = grep {!/$j/} @index4;
							}
							else
							{
								if (! grep(/$j/, @temp_index4)) 
								{
									push @temp_index4, $j;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$k/, @temp_index3)) 
								{
									push @temp_index3, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index4 == 0)
						{
							@index3 = grep {!/$k/} @index3;
						}
					}
					@index3 = @temp_index3;
					@index4 = @temp_index4;
				
					if ($n == $#matchset1)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset2);
						}
					}
				}
				elsif ($n == $#matchset1)
				{
					if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
					{
						push (@no_matching_pair, @matchset1);
						push (@no_matching_pair, @matchset2);
					}
				}
			}
		}
	}
	else
	{
		foreach my $item (@matchset1)
		{
			my @set1 = split ("\t", $item);
			my $m = $set1[0];
			my $s = $set1[1];
			push (@match_maid_chr_1, $m.$s);
			push (@match_chr_1, $m);
			$seen{$m} = 1;
			$seen{$m.$s} = 1;
		}
		for (my $n = 0; $n < @matchset2; $n++)
		{
			my @set2 = split ("\t", $matchset2[$n]);
			my $m = $set2[0];
			my $s = $set2[1];
			if ($s eq "+") {$s = "-";}
			elsif ($s eq "-") {$s = "+";}
			push (@match_maid_chr_2, $m.$s);
			push (@match_chr_2, $m);
			if ($seen{$m.$s}) # 5 +  ->  5 -
			{
				@index1= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
				@index2= grep {$match_maid_chr_2[$_] eq $m.$s} 0..$#match_maid_chr_2;
				
				if ((@index1 == 1) && (@index2 > 1)) # more than one match
				{
					my @a1 = split ("\t", $matchset1[$index1[0]]);
					$toligo1 = $a1[6];
					chomp $toligo1;
					my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
					$toligo1 = $toligos1[0];
					if (defined $toligos1[1])
					{
						my $toligo1_num = $toligos1[1];
					}
					my $tstart1 = $a1[2];
					my $tend1 = $a1[3];
					my @sorted_index2 = sort {$a <=> $b} @index2;
					foreach my $k (@sorted_index2)
					{	
						my @a2 = split ("\t", $matchset2[$k]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my $seqlen1 = abs($tend2 - $tstart1);
						my $seqlen2 = abs($tend1 - $tstart2);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
						if ($seqlen > 600)
						{
							@index2 = grep {!/$k/} @index2;
						}
						else
						{
							if (! grep(/$k/, @temp_index2)) 
							{
								push @temp_index2, $k;
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index1[0]/, @temp_index1))
							{
								push @temp_index1, $index1[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
					}
					if (@index2 == 0)
					{
						@index1 = grep {!/$index1[0]/} @index1;
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1 > 1) && (@index2 == 1))
				{
					my @a2 = split ("\t", $matchset2[$index2[0]]);
					$toligo2 = $a2[6];
					chomp $toligo2;
					my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
					$toligo2 = $toligos2[0];
					if (defined $toligos2[1])
					{
						my $toligo2_num = $toligos2[1];
					}
					my $tstart2 = $a2[2];
					my $tend2 = $a2[3];
					my @sorted_index1 = sort {$a <=> $b} @index1;
					foreach my $k (@sorted_index1)
					{
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my $seqlen1 = abs($tend1 - $tstart2);
						my $seqlen2 = abs($tend2 - $tstart1);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 600)
						{
							@index1 = grep {!/$k/} @index1;
						}
						else
						{
							if (! grep(/$k/, @temp_index1)) 
							{
								push @temp_index1, $k;
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index2[0]/, @temp_index2))
							{
								push @temp_index2, $index2[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
					}
					if (@index1 == 0)
					{
						@index2 = grep {!/$index2[0]/} @index2;
					}	
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1== 1) && (@index2 == 1))
				{
					my @a1 = split ("\t", $matchset1[$index1[0]]);
					$toligo1 = $a1[6];
					chomp $toligo1;
					my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
					$toligo1 = $toligos1[0];
					if (defined $toligos1[1])
					{
						my $toligo1_num = $toligos1[1];
					}
					my $tstart1 = $a1[2];
					my $tend1 = $a1[3];

					my @a2 = split ("\t", $matchset2[$index2[0]]);
					$toligo2 = $a2[6];
					chomp $toligo2;
					my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
					$toligo2 = $toligos2[0];
					if (defined $toligos2[1])
					{
						my $toligo2_num = $toligos2[1];
					}
					my $tstart2 = $a2[2];
					my $tend2 = $a2[3];

					my $seqlen1 = abs($tend2 - $tstart1);
					my $seqlen2 = abs($tend1 - $tstart2);
					my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
					if ($seqlen > 600)
					{
						shift @index1;
						shift @index2;
					}
					else
					{
						if (! grep(/$index1[0]/, @temp_index1))
						{
							push @temp_index1, $index1[0];
							if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
				
						}
						if (! grep(/$index2[0]/, @temp_index2))
						{
							push @temp_index2, $index2[0];
							if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
						}
					}
					@index1 = @temp_index1;
					@index2 = @temp_index2;
				}

				elsif ((@index1> 1) && (@index2 > 1))
				{
					my @tindex1 = @index1;
					my @tindex2 = @index2;
					my @sorted_index1 = sort {$a <=> $b} @index1;
					foreach my $k (@sorted_index1)
					{
						@tindex2 = @index2;
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index2 = sort {$a <=> $b} @index2;
						foreach my $j (@sorted_index2)
						{
							my @a2 = split ("\t", $matchset2[$j]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 600)
							{
								@tindex2 = grep {!/$j/} @tindex2;
							}
							else
							{
								if (! grep(/$j/, @temp_index2)) 
								{
									push @temp_index2, $j;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$k/, @temp_index1)) 
								{
									push @temp_index1, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (!@tindex2)
						{
							@tindex1 = grep {!/$k/} @tindex1;
						}										
					}
					 @index1 = @temp_index1;
					 @index2 = @temp_index2;
				}
			
				if ($n == $#matchset2)
				{
					if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
					{
						push (@no_matching_pair, @matchset1);
						push (@no_matching_pair, @matchset2);
					}
				}
			}
			else # 5 -  -> 5 -
			{
				if ($seen{$m})
				{
					@index3= grep {$match_chr_1[$_] eq $m} 0..$#match_chr_1;
					if ($s eq "+") {$s = "-";}
					elsif ($s eq "-") {$s = "+";}
					@index4= grep {$match_chr_2[$_] eq $m} 0..$#match_chr_2;
					my @sorted_index3 = sort {$a <=> $b} @index3;
					foreach my $k (@sorted_index3)
					{
						my @a1 = split ("\t", $matchset1[$k]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index4 = sort {$a <=> $b} @index4;
						foreach my $j (@sorted_index4)
						{
							my @a2 = split ("\t", $matchset2[$j]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 600)
							{
								@index4 = grep {!/$j/} @index4;
							}
							else
							{
								if (! grep(/$j/, @temp_index4)) 
								{
									push @temp_index4, $j;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$k/, @temp_index3)) 
								{
									push @temp_index3, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index4 == 0)
						{
							@index3 = grep {!/$k/} @index3;
						}
					}
					@index3 = @temp_index3;
					@index4 = @temp_index4;
				
					if ($n == $#matchset2)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset2);
						}
					}
				}
				elsif ($n == $#matchset2)
				{
					if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
					{
						push (@no_matching_pair, @matchset1);
						push (@no_matching_pair, @matchset2);
					}
				}
			}
		}
	}
	
	@index1 = @temp_index1;
	@index2 = @temp_index2;
	@index3 = @temp_index3;
	@index4 = @temp_index4;
	
	if (@temp_index1 != 0)
	{
		my @temp = split ("\t", $matchset1[$temp_index1[$#temp_index1]]);
		$tempstart = $temp[2];
	}
	elsif (@temp_index2 != 0)
	{
		my @temp = split ("\t", $matchset2[$temp_index2[$#temp_index2]]);
		$tempstart = $temp[2];
	}
	
	if ((@temp_index1 > 1) || (@temp_index2 > 1))
	{
		foreach my $t (@temp_index1)
		{
			push (@no_matching_pair, $matchset1[$t]);
		}
		foreach my $v (@temp_index2)
		{
			push (@no_matching_pair, $matchset2[$v]);
		}
		(@index1, @index2) = ();
	}
	
	if (@temp_index3 != 0)
	{
		my @temp = split ("\t", $matchset1[$temp_index3[$#temp_index3]]);
		$tempstart = $temp[2];
	}
	elsif (@temp_index4 != 0)
	{
		my @temp = split ("\t", $matchset2[$temp_index4[$#temp_index4]]);
		$tempstart = $temp[2];
	}
	
	if ((@temp_index3 > 1) || (@temp_index4 > 1))
	{
		foreach my $t (@temp_index3)
		{
			push (@no_matching_pair, $matchset1[$t]);
		}
		foreach my $v (@temp_index4)
		{
			push (@no_matching_pair, $matchset2[$v]);
		}
		(@index3, @index4) = ();
	}
	
	foreach my $j (@index1)
	{
		push (@no_matching_pair, $matchset1[$j]);
	}
	foreach my $p (@index2)
	{
		push (@no_matching_pair, $matchset2[$p]);
	}
	
	if(!@index1)
	{
		foreach my $k (@index3)
		{
			push (@no_matching_pair, $matchset1[$k]);
		}
	}
	if (!@index2)
	{
		foreach my $n (@index4)
		{
			push (@no_matching_pair, $matchset2[$n]);
		}
	}
}

#TDF TDP TDR TUF TUP TUR 3 letters, the distance between a pair (eg. TDF & TDR) < 500bp)
sub compare_taq
{
	my (@index1, @index2, @index3, @index4, @index5, @index6, @matchset1, @matchset2, @matchset3, @match_maid_chr_1, @match_maid_chr_2, @match_maid_chr_3, @match_chr_1, @match_chr_2, @match_chr_3, @temp_index1, @temp_index2, @temp_index3, @temp_index4, @temp_index5, @temp_index6) = ();
	my $matchset1ref = shift;
	my $matchset2ref = shift;
	my $matchset3ref = shift;
	@matchset1 = @{$matchset1ref};
	@matchset2 = @{$matchset2ref};
	@matchset3 = @{$matchset3ref};
	
	if (scalar(@matchset1) >= scalar(@matchset3))
	{
		if (scalar(@matchset2) >= scalar(@matchset3))
		{
			foreach my $item (@matchset3)
			{
				my @set3 = split ("\t", $item);
				my $m = $set3[0];
				my $s = $set3[1];
				push (@match_maid_chr_3, $m.$s);
				push (@match_chr_3, $m);
				$seen{$m} = 1;
				$seen{$m.$s} = 1;
			}	

			for (my $n = 0; $n < @matchset1; $n++)
			{
				my @set1 = split ("\t", $matchset1[$n]);
				my $m = $set1[0];
				my $s = $set1[1];
				if ($s eq "+") {$s = "-";}
				elsif ($s eq "-") {$s = "+";}
				push (@match_maid_chr_1, $m.$s);
				push (@match_chr_1, $m);
				if ($seen{$m.$s}) # 5 +  ->  5 -
				{
					@index1= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
					@index2= grep {$match_maid_chr_3[$_] eq $m.$s} 0..$#match_maid_chr_3;

					if ((@index1 == 1) && (@index2 > 1)) # more than one match
					{
						my @a1 = split ("\t", $matchset1[$index1[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index2 = sort {$a <=> $b} @index2;
						foreach my $k (@sorted_index2)
						{
							my @a3 = split ("\t", $matchset3[$k]);
							$toligo3 = $a3[6];
							chomp $toligo3;
							my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
							$toligo3 = $toligos3[0];
							if (defined $toligos3[1])
							{
								my $toligo3_num = $toligos3[1];
							}
							my $tstart3 = $a3[2];
							my $tend3 = $a3[3];
							my $seqlen1 = abs($tend3 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart3);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index2 = grep {!/$k/} @index2;
							}
							else
							{
								if (! grep(/$k/, @temp_index2)) 
								{
									push @temp_index2, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index1[0]/, @temp_index1))
								{
									push @temp_index1, $index1[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index2 == 0)
						{
							@index1 = grep {!/$index1[0]/} @index1;
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1 > 1) && (@index2 == 1))
					{
						my @a3 = split ("\t", $matchset3[$index2[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my @sorted_index1 = sort {$a <=> $b} @index1;
						foreach my $k (@sorted_index1)
						{
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my $seqlen1 = abs($tend1 - $tstart3);
							my $seqlen2 = abs($tend3 - $tstart1);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

							if ($seqlen > 600)
							{
								@index1 = grep {!/$k/} @index1;
							}
							else
							{
								if (! grep(/$k/, @temp_index1)) 
								{
									push @temp_index1, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index2[0]/, @temp_index2))
								{
									push @temp_index2, $index2[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index1 == 0)
						{
							@index2 = grep {!/$index2[0]/} @index2;
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1== 1) && (@index2 == 1))
					{
						my @a1 = split ("\t", $matchset1[$index1[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];

						my @a3 = split ("\t", $matchset3[$index2[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];

						my $seqlen1 = abs($tend3 - $tstart1);
						my $seqlen2 = abs($tend1 - $tstart3);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 600)
						{
							shift @index1;
							shift @index2;
						}
						else
						{
							if (! grep(/$index1[0]/, @temp_index1))
							{
								push @temp_index1, $index1[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index2[0]/, @temp_index2))
							{
								push @temp_index2, $index2[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1> 1) && (@index2 > 1))
					{
						my @tindex1 = @index1;
						my @tindex2 = @index2;
						my @sorted_index1 = sort {$a <=> $b} @index1;
						foreach my $k (@sorted_index1)
						{
							@tindex2 = @index2;
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my @sorted_index2 = sort {$a <=> $b} @index2;
							foreach my $j (@sorted_index2)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart1);
								my $seqlen2 = abs($tend1 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 600)
								{
									@tindex2 = grep {!/$j/} @tindex2;
								}
								else
								{
									if (! grep(/$j/, @temp_index2)) 
									{
										push @temp_index2, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index1)) 
									{
										push @temp_index1, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (!@tindex2)
							{
								@tindex1 = grep {!/$k/} @tindex1;
							}										
						}
						 @index1 = @temp_index1;
						 @index2 = @temp_index2;
					}

					if ($n == $#matchset1)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset3);
						}
					}

				}
				else # 5 -  -> 5 -
				{
					if ($seen{$m})
					{
						@index3= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
						if ($s eq "+") {$s = "-";}
						elsif ($s eq "-") {$s = "+";}
						@index4= grep {$match_maid_chr_3[$_] eq $m.$s} 0..$#match_maid_chr_3;
						my @sorted_index3 = sort {$a <=> $b} @index3;
						foreach my $k (@sorted_index3)
						{
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my @sorted_index4 = sort {$a <=> $b} @index4;
							foreach my $j (@sorted_index4)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart1);
								my $seqlen2 = abs($tend1 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@index4 = grep {!/$j/} @index4;
								}
								else
								{
									if (! grep(/$j/, @temp_index4)) 
									{
										push @temp_index4, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index3)) 
									{
										push @temp_index3, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@index4 == 0)
							{
								@index3 = grep {!/$k/} @index3;
							}
						}
						@index3 = @temp_index3;
						@index4 = @temp_index4;

						if (@temp_index3 == 1)
						{
							my @temp = split ("\t", $matchset1[$temp_index3[0]]);
							$tempstart = $temp[2];
						}
						elsif (@temp_index4 == 1)
						{
							my @temp = split ("\t", $matchset3[$temp_index4[0]]);
							$tempstart = $temp[2];
						}

						if ((@temp_index3 > 1) || (@temp_index4 > 1))
						{
							foreach my $t (@temp_index3)
							{
								push (@no_matching_pair, $matchset1[$t]);
							}
							foreach my $v (@temp_index4)
							{
								push (@no_matching_pair, $matchset3[$v]);
							}
							(@index3, @index4) = ();
						}

						if ($n == $#matchset1)
						{
							if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
							{
								push (@no_matching_pair, @matchset1);
								push (@no_matching_pair, @matchset3);
							}
						}
					}
					elsif ($n == $#matchset1)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset3);
						}
					}
				}
			}

			for (my $n = 0; $n < @matchset2; $n++)
			{
				my @set2 = split ("\t", $matchset2[$n]);
				my $m = $set2[0];
				my $s = $set2[1];
				push (@match_maid_chr_2, $m.$s);
				push (@match_chr_2, $m);
				if ($seen{$m}) # 5   ->  5 
				{
					@index5= grep {$match_chr_2[$_] eq $m} 0..$#match_chr_2;
					@index6= grep {$match_chr_3[$_] eq $m} 0..$#match_chr_3;

					if ((@index5 == 1) && (@index6 > 1)) # more than one match
					{
						my @a2 = split ("\t", $matchset2[$index5[0]]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my @sorted_index6 = sort {$a <=> $b} @index6;
						foreach my $k (@sorted_index6)
						{
							my @a3 = split ("\t", $matchset3[$k]);
							$toligo3 = $a3[6];
							chomp $toligo3;
							my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
							$toligo3 = $toligos3[0];
							if (defined $toligos3[1])
							{
								my $toligo3_num = $toligos3[1];
							}
							my $tstart3 = $a3[2];
							my $tend3 = $a3[3];
							my $seqlen1 = abs($tend3 - $tstart2);
							my $seqlen2 = abs($tend2 - $tstart3);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index6 = grep {!/$k/} @index6;
							}
							else
							{
								if (! grep(/$k/, @temp_index6)) 
								{
									push @temp_index6, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index5[0]/, @temp_index5))
								{
									push @temp_index5, $index5[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}

						}
						if (@index6 == 0)
						{
							@index5 = grep {!/$index5[0]/} @index5;
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5 > 1) && (@index6 == 1))
					{
						my @a3 = split ("\t", $matchset3[$index6[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my @sorted_index5 = sort {$a <=> $b} @index5;
						foreach my $k (@sorted_index5)
						{
							my @a2 = split ("\t", $matchset2[$k]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart3);
							my $seqlen2 = abs($tend3 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index5 = grep {!/$k/} @index5;
							}
							else
							{
								if (! grep(/$k/, @temp_index5)) 
								{
									push @temp_index5, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index6[0]/, @temp_index6))
								{
									push @temp_index6, $index6[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index5 == 0)
						{
							@index6 = grep {!/$index6[0]/} @index6;
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5 == 1) && (@index6 == 1))
					{
						my @a2 = split ("\t", $matchset2[$index5[0]]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my @a3 = split ("\t", $matchset3[$index6[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my $seqlen1 = abs($tend3 - $tstart2);
						my $seqlen2 = abs($tend2 - $tstart3);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 500)
						{
							shift @index5;
							shift @index6;
						}
						else
						{
							if (! grep(/$index5[0]/, @temp_index5))
							{
								push @temp_index5, $index5[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index6[0]/, @temp_index6))
							{
								push @temp_index6, $index6[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5> 1) && (@index6 > 1))
					{
						my @tindex5 = @index5;
						my @tindex6 = @index6;
						my @sorted_index5 = sort {$a <=> $b} @index5;
						foreach my $k (@sorted_index5)
						{
							@tindex6 = @index6;
							my @a2 = split ("\t", $matchset2[$k]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my @sorted_index6 = sort {$a <=> $b} @index6;
							foreach my $j (@sorted_index6)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart2);
								my $seqlen2 = abs($tend2 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@tindex6 = grep {!/$j/} @tindex6;
								}
								else
								{

									if (! grep(/$j/, @temp_index6)) 
									{
										push @temp_index6, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index5)) 
									{
										push @temp_index5, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@tindex6 == 0)
							{
								@tindex5 = grep {!/$k/} @tindex5;
							}
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					if ($n == $#matchset2)
					{
						if ((!@index5) && (!@index6))
						{
							push (@no_matching_pair, @matchset2);
							push (@no_matching_pair, @matchset3);
						}
					}
				}
				elsif ($n == $#matchset2)
				{
					if ((!@index5) && (!@index6))
					{
						push (@no_matching_pair, @matchset2);
						push (@no_matching_pair, @matchset3);
					}
				}
			}
		}
		else
		{
			foreach my $item (@matchset3)
			{
				my @set3 = split ("\t", $item);
				my $m = $set3[0];
				my $s = $set3[1];
				push (@match_maid_chr_3, $m.$s);
				push (@match_chr_3, $m);
				$seen{$m} = 1;
				$seen{$m.$s} = 1;
			}	

			for (my $n = 0; $n < @matchset1; $n++)
			{
				my @set1 = split ("\t", $matchset1[$n]);
				my $m = $set1[0];
				my $s = $set1[1];
				if ($s eq "+") {$s = "-";}
				elsif ($s eq "-") {$s = "+";}
				push (@match_maid_chr_1, $m.$s);
				push (@match_chr_1, $m);
				if ($seen{$m.$s}) # 5 +  ->  5 -
				{
					@index1= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
					@index2= grep {$match_maid_chr_3[$_] eq $m.$s} 0..$#match_maid_chr_3;

					if ((@index1 == 1) && (@index2 > 1)) # more than one match
					{
						my @a1 = split ("\t", $matchset1[$index1[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index2 = sort {$a <=> $b} @index2;
						foreach my $k (@sorted_index2)
						{
							my @a3 = split ("\t", $matchset3[$k]);
							$toligo3 = $a3[6];
							chomp $toligo3;
							my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
							$toligo3 = $toligos3[0];
							if (defined $toligos3[1])
							{
								my $toligo3_num = $toligos3[1];
							}
							my $tstart3 = $a3[2];
							my $tend3 = $a3[3];
							my $seqlen1 = abs($tend3 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart3);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index2 = grep {!/$k/} @index2;
							}
							else
							{
								if (! grep(/$k/, @temp_index2)) 
								{
									push @temp_index2, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index1[0]/, @temp_index1))
								{
									push @temp_index1, $index1[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index2 == 0)
						{
							@index1 = grep {!/$index1[0]/} @index1;
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1 > 1) && (@index2 == 1))
					{
						my @a3 = split ("\t", $matchset3[$index2[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my @sorted_index1 = sort {$a <=> $b} @index1;
						foreach my $k (@sorted_index1)
						{
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my $seqlen1 = abs($tend1 - $tstart3);
							my $seqlen2 = abs($tend3 - $tstart1);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

							if ($seqlen > 500)
							{
								@index1 = grep {!/$k/} @index1;
							}
							else
							{
								if (! grep(/$k/, @temp_index1)) 
								{
									push @temp_index1, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index2[0]/, @temp_index2))
								{
									push @temp_index2, $index2[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index1 == 0)
						{
							@index2 = grep {!/$index2[0]/} @index2;
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1== 1) && (@index2 == 1))
					{
						my @a1 = split ("\t", $matchset1[$index1[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];

						my @a3 = split ("\t", $matchset3[$index2[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];

						my $seqlen1 = abs($tend3 - $tstart1);
						my $seqlen2 = abs($tend1 - $tstart3);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 500)
						{
							shift @index1;
							shift @index2;
						}
						else
						{
							if (! grep(/$index1[0]/, @temp_index1))
							{
								push @temp_index1, $index1[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index2[0]/, @temp_index2))
							{
								push @temp_index2, $index2[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1> 1) && (@index2 > 1))
					{
						my @tindex1 = @index1;
						my @tindex2 = @index2;
						my @sorted_index1 = sort {$a <=> $b} @index1;
						foreach my $k (@sorted_index1)
						{
							@tindex2 = @index2;
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my @sorted_index2 = sort {$a <=> $b} @index2;
							foreach my $j (@sorted_index2)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart1);
								my $seqlen2 = abs($tend1 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@tindex2 = grep {!/$j/} @tindex2;
								}
								else
								{
									if (! grep(/$j/, @temp_index2)) 
									{
										push @temp_index2, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index1)) 
									{
										push @temp_index1, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@tindex2 == 0)
							{
								@tindex1 = grep {!/$k/} @tindex1;
							}
						}
						 @index1 = @temp_index1;
						 @index2 = @temp_index2;
					}

					if ($n == $#matchset1)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset3);
						}
					}

				}
				else # 5 -  -> 5 -
				{
					if ($seen{$m})
					{
						@index3= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
						if ($s eq "+") {$s = "-";}
						elsif ($s eq "-") {$s = "+";}
						@index4= grep {$match_maid_chr_3[$_] eq $m.$s} 0..$#match_maid_chr_3;
						my @sorted_index3 = sort {$a <=> $b} @index3;
						foreach my $k (@sorted_index3)
						{
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my @sorted_index4 = sort {$a <=> $b} @index4;
							foreach my $j (@sorted_index4)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart1);
								my $seqlen2 = abs($tend1 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@index4 = grep {!/$j/} @index4;
								}
								else
								{
									if (! grep(/$j/, @temp_index4)) 
									{
										push @temp_index4, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index3)) 
									{
										push @temp_index3, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@index4 == 0)
							{
								@index3 = grep {!/$k/} @index3;
							}
						}
						@index3 = @temp_index3;
						@index4 = @temp_index4;

						if ($n == $#matchset1)
						{
							if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
							{
								push (@no_matching_pair, @matchset1);
								push (@no_matching_pair, @matchset3);
							}
						}
					}
					elsif ($n == $#matchset1)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset3);
						}
					}
				}
			}

			@match_maid_chr_3 = ();
			@match_chr_3 = ();

			foreach my $item (@matchset2)
			{
				my @set2 = split ("\t", $item);
				my $m = $set2[0];
				my $s = $set2[1];
				push (@match_maid_chr_2, $m.$s);
				push (@match_chr_2, $m);
				$seen{$m} = 1;
				$seen{$m.$s} = 1;
			}

			for (my $n = 0; $n < @matchset3; $n++)
			{
				my @set3 = split ("\t", $matchset3[$n]);
				my $m = $set3[0];
				my $s = $set3[1];
				push (@match_maid_chr_3, $m.$s);
				push (@match_chr_3, $m);
				if ($seen{$m}) # 5   ->  5 
				{
					@index5= grep {$match_chr_2[$_] eq $m} 0..$#match_chr_2;
					@index6= grep {$match_chr_3[$_] eq $m} 0..$#match_chr_3;

					if ((@index5 == 1) && (@index6 > 1)) # more than one match
					{
						my @a2 = split ("\t", $matchset2[$index5[0]]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my @sorted_index6 = sort {$a <=> $b} @index6;
						foreach my $k (@sorted_index6)
						{
							my @a3 = split ("\t", $matchset3[$k]);
							$toligo3 = $a3[6];
							chomp $toligo3;
							my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
							$toligo3 = $toligos3[0];
							if (defined $toligos3[1])
							{
								my $toligo3_num = $toligos3[1];
							}
							my $tstart3 = $a3[2];
							my $tend3 = $a3[3];
							my $seqlen1 = abs($tend3 - $tstart2);
							my $seqlen2 = abs($tend2 - $tstart3);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index6 = grep {!/$k/} @index6;
							}
							else
							{
								if (! grep(/$k/, @temp_index6)) 
								{
									push @temp_index6, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index5[0]/, @temp_index5))
								{
									push @temp_index5, $index5[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index6 == 0)
						{
							@index5 = grep {!/$index5[0]/} @index5;
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5 > 1) && (@index6 == 1))
					{
						my @a3 = split ("\t", $matchset3[$index6[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my @sorted_index5 = sort {$a <=> $b} @index5;
						foreach my $k (@sorted_index5)
						{
							my @a2 = split ("\t", $matchset2[$k]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart3);
							my $seqlen2 = abs($tend3 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index5 = grep {!/$k/} @index5;
							}
							else
							{
								if (! grep(/$k/, @temp_index5)) 
								{
									push @temp_index5, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index6[0]/, @temp_index6))
								{
									push @temp_index6, $index6[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index5 == 0)
						{
							@index6 = grep {!/$index6[0]/} @index6;
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5 == 1) && (@index6 == 1))
					{
						my @a2 = split ("\t", $matchset2[$index5[0]]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my @a3 = split ("\t", $matchset3[$index6[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my $seqlen1 = abs($tend3 - $tstart2);
						my $seqlen2 = abs($tend2 - $tstart3);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 500)
						{
							shift @index5;
							shift @index6;
						}
						else
						{
							if (! grep(/$index5[0]/, @temp_index5))
							{
								push @temp_index5, $index5[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index6[0]/, @temp_index6))
							{
								push @temp_index6, $index6[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5> 1) && (@index6 > 1))
					{
						my @tindex5 = @index5;
						my @tindex6 = @index6;
						my @sorted_index5 = sort {$a <=> $b} @index5;
						foreach my $k (@sorted_index5)
						{
							@tindex6 = @index6;
							my @a2 = split ("\t", $matchset2[$k]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my @sorted_index6 = sort {$a <=> $b} @index6;
							foreach my $j (@sorted_index6)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart2);
								my $seqlen2 = abs($tend2 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@tindex6 = grep {!/$j/} @tindex6;
								}
								else
								{
									if (! grep(/$j/, @temp_index6)) 
									{
										push @temp_index6, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index5)) 
									{
										push @temp_index5, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@tindex6 == 0)
							{
								@tindex5 = grep {!/$k/} @tindex5;
							}
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					if ($n == $#matchset3)
					{
						if ((!@index5) && (!@index6))
						{
							push (@no_matching_pair, @matchset2);
							push (@no_matching_pair, @matchset3);
						}
					}
				}
				elsif ($n == $#matchset3)
				{
					if ((!@index5) && (!@index6))
					{
						push (@no_matching_pair, @matchset2);
						push (@no_matching_pair, @matchset3);
					}
				}
			}
		}

		@index1 = @temp_index1;
		@index2 = @temp_index2;
		@index3 = @temp_index3;
		@index4 = @temp_index4;
		@index5 = @temp_index5;
		@index6 = @temp_index6;

		if (@temp_index1 != 0)
		{
			my @temp = split ("\t", $matchset1[$temp_index1[$#temp_index1]]);
			$tempstart = $temp[2];
		}
		elsif (@temp_index2 != 0)
		{
			my @temp = split ("\t", $matchset3[$temp_index2[$#temp_index2]]);
			$tempstart = $temp[2];
		}

		if ((@temp_index1 > 1) || (@temp_index2 > 1))
		{
			foreach my $t (@temp_index1)
			{
				push (@no_matching_pair, $matchset1[$t]);
			}
			foreach my $v (@temp_index2)
			{
				push (@no_matching_pair, $matchset3[$v]);
			}
			(@index1, @index2) = ();
		}

		if (@temp_index3 != 0)
		{
			my @temp = split ("\t", $matchset1[$temp_index3[$#temp_index3]]);
			$tempstart = $temp[2];
		}
		elsif (@temp_index4 != 0)
		{
			my @temp = split ("\t", $matchset3[$temp_index4[$#temp_index4]]);
			$tempstart = $temp[2];
		}

		if ((@temp_index3 > 1) || (@temp_index4 > 1))
		{
			foreach my $t (@temp_index3)
			{
				push (@no_matching_pair, $matchset1[$t]);
			}
			foreach my $v (@temp_index4)
			{
				push (@no_matching_pair, $matchset3[$v]);
			}
			(@index3, @index4) = ();
		}

		if (@temp_index5 != 0)
		{
			my @temp = split ("\t", $matchset2[$temp_index5[$#temp_index5]]);
			$tempstart = $temp[2];
		}
		elsif (@temp_index6 != 0)
		{
			my @temp = split ("\t", $matchset3[$temp_index6[$#temp_index6]]);
			$tempstart = $temp[2];
		}

		if ((@temp_index5 > 1) || (@temp_index6 > 1))
		{
			foreach my $t (@temp_index5)
			{
				push (@no_matching_pair, $matchset2[$t]);
			}
			foreach my $v (@temp_index6)
			{
				push (@no_matching_pair, $matchset3[$v]);
			}
			(@index5, @index6) = ();
		}

		foreach my $j (@index1)
		{
			open (OUTPUT, '>>match.txt');
			print OUTPUT ("$matchset1[$j]");
		}
		foreach my $p (@index2)
		{		
			open (OUTPUT, '>>match.txt');
			print OUTPUT ("$matchset3[$p]");
			close (OUTPUT);
		}
		if (!@index1)
		{
			foreach my $k (@index3)
			{
				open (OUTPUT2, '>>match.txt');
				print OUTPUT2 ("$matchset1[$k]");
			}
		}
		if (!@index2)
		{
			foreach my $n (@index4)
			{
				open (OUTPUT2, '>>match.txt');
				print OUTPUT2 ("$matchset3[$n]");
				close (OUTPUT2);
			}
		}						
		foreach my $j (@index5)
		{
			push (@no_matching_pair, $matchset2[$j]);
		}
		if (!@index2)
		{
			foreach my $j (@index6)
			{
				push (@no_matching_pair, $matchset3[$j]);
			}
		}
	}							
	else
	{
		if (scalar(@matchset2) >= scalar(@matchset3))
		{
			foreach my $item (@matchset1)
			{
				my @set1 = split ("\t", $item);
				my $m = $set1[0];
				my $s = $set1[1];
				push (@match_maid_chr_1, $m.$s);
				push (@match_chr_1, $m);
				$seen{$m} = 1;
				$seen{$m.$s} = 1;
			}	

			for (my $n = 0; $n < @matchset3; $n++)
			{
				my @set3 = split ("\t", $matchset3[$n]);
				my $m = $set3[0];
				my $s = $set3[1];
				if ($s eq "+") {$s = "-";}
				elsif ($s eq "-") {$s = "+";}
				push (@match_maid_chr_3, $m.$s);
				push (@match_chr_3, $m);
				if ($seen{$m.$s}) # 5 +  ->  5 -
				{
					@index1= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
					@index2= grep {$match_maid_chr_3[$_] eq $m.$s} 0..$#match_maid_chr_3;

					if ((@index1 == 1) && (@index2 > 1)) # more than one match
					{
						my @a1 = split ("\t", $matchset1[$index1[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index2 = sort {$a <=> $b} @index2;
						foreach my $k (@sorted_index2)
						{
							my @a3 = split ("\t", $matchset3[$k]);
							$toligo3 = $a3[6];
							chomp $toligo3;
							my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
							$toligo3 = $toligos3[0];
							if (defined $toligos3[1])
							{
								my $toligo3_num = $toligos3[1];
							}
							my $tstart3 = $a3[2];
							my $tend3 = $a3[3];
							my $seqlen1 = abs($tend3 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart3);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index2 = grep {!/$k/} @index2;
							}
							else
							{
								if (! grep(/$k/, @temp_index2)) 
								{
									push @temp_index2, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index1[0]/, @temp_index1))
								{
									push @temp_index1, $index1[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index2 == 0)
						{
							@index1 = grep {!/$index1[0]/} @index1;
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1 > 1) && (@index2 == 1))
					{
						my @a3 = split ("\t", $matchset3[$index2[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my @sorted_index1 = sort {$a <=> $b} @index1;
						foreach my $k (@sorted_index1)
						{
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my $seqlen1 = abs($tend1 - $tstart3);
							my $seqlen2 = abs($tend3 - $tstart1);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

							if ($seqlen > 500)
							{
								@index1 = grep {!/$k/} @index1;
							}
							else
							{
								if (! grep(/$k/, @temp_index1)) 
								{
									push @temp_index1, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index2[0]/, @temp_index2))
								{
									push @temp_index2, $index2[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index1 == 0)
						{
							@index2 = grep {!/$index2[0]/} @index2;
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1== 1) && (@index2 == 1))
					{
						my @a1 = split ("\t", $matchset1[$index1[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];

						my @a3 = split ("\t", $matchset3[$index2[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];

						my $seqlen1 = abs($tend3 - $tstart1);
						my $seqlen2 = abs($tend1 - $tstart3);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 500)
						{
							shift @index1;
							shift @index2;
						}
						else
						{
							if (! grep(/$index1[0]/, @temp_index1))
							{
								push @temp_index1, $index1[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index2[0]/, @temp_index2))
							{
								push @temp_index2, $index2[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1> 1) && (@index2 > 1))
					{
						my @tindex1 = @index1;
						my @tindex2 = @index2;
						my @sorted_index1 = sort {$a <=> $b} @index1;
						foreach my $k (@sorted_index1)
						{
							@tindex2 = @index2;
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my @sorted_index2 = sort {$a <=> $b} @index2;
							foreach my $j (@sorted_index2)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart1);
								my $seqlen2 = abs($tend1 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@tindex2 = grep {!/$j/} @tindex2;
								}
								else
								{
									if (! grep(/$j/, @temp_index2)) 
									{
										push @temp_index2, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index1)) 
									{
										push @temp_index1, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@tindex2 == 0)
							{
								@tindex1 = grep {!/$k/} @tindex1;
							}
						}
						 @index1 = @temp_index1;
						 @index2 = @temp_index2;
					}

					if ($n == $#matchset3)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset3);
						}
					}
				}
				else # 5 -  -> 5 -
				{
					if ($seen{$m})
					{
						@index3= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
						if ($s eq "+") {$s = "-";}
						elsif ($s eq "-") {$s = "+";}
						@index4= grep {$match_maid_chr_3[$_] eq $m.$s} 0..$#match_maid_chr_3;
						my @sorted_index3 = sort {$a <=> $b} @index3;
						foreach my $k (@sorted_index3)
						{
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my @sorted_index4 = sort {$a <=> $b} @index4;
							foreach my $j (@sorted_index4)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart1);
								my $seqlen2 = abs($tend1 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@index4 = grep {!/$j/} @index4;
								}
								else
								{
									if (! grep(/$j/, @temp_index4)) 
									{
										push @temp_index4, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index3)) 
									{
										push @temp_index3, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@index4 == 0)
							{
								@index3 = grep {!/$k/} @index3;
							}
						}
						@index3 = @temp_index3;
						@index4 = @temp_index4;

						if ($n == $#matchset3)
						{
							if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
							{
								push (@no_matching_pair, @matchset1);
								push (@no_matching_pair, @matchset3);
							}
						}
					}
					elsif ($n == $#matchset3)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset3);
						}
					}
				}
			}

			for (my $n = 0; $n < @matchset2; $n++)
			{
				my @set2 = split ("\t", $matchset2[$n]);
				my $m = $set2[0];
				my $s = $set2[1];
				push (@match_maid_chr_2, $m.$s);
				push (@match_chr_2, $m);
				if ($seen{$m}) # 5   ->  5 
				{
					@index5= grep {$match_chr_2[$_] eq $m} 0..$#match_chr_2;
					@index6= grep {$match_chr_1[$_] eq $m} 0..$#match_chr_1;											

					if ((@index5 == 1) && (@index6 > 1)) # more than one match
					{
						my @a2 = split ("\t", $matchset2[$index5[0]]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my @sorted_index6 = sort {$a <=> $b} @index6;
						foreach my $k (@sorted_index6)
						{
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index6 = grep {!/$k/} @index6;
							}
							else
							{
								if (! grep(/$k/, @temp_index6)) 
								{
									push @temp_index6, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index5[0]/, @temp_index5))
								{
									push @temp_index5, $index5[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index6 == 0)
						{
							@index5 = grep {!/$index5[0]/} @index5;
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5 > 1) && (@index6 == 1))
					{
						my @a1 = split ("\t", $matchset1[$index6[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index5 = sort {$a <=> $b} @index5;
						foreach my $k (@sorted_index5)
						{
							my @a2 = split ("\t", $matchset2[$k]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index5 = grep {!/$k/} @index5;
							}
							else
							{
								if (! grep(/$k/, @temp_index5)) 
								{
									push @temp_index5, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index6[0]/, @temp_index6))
								{
									push @temp_index6, $index6[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index5 == 0)
						{
							@index6 = grep {!/$index6[0]/} @index6;
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5 == 1) && (@index6 == 1))
					{
						my @a2 = split ("\t", $matchset2[$index5[0]]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my @a1 = split ("\t", $matchset1[$index6[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my $seqlen1 = abs($tend2 - $tstart1);
						my $seqlen2 = abs($tend1 - $tstart2);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 500)
						{
							shift @index5;
							shift @index6;
						}
						else
						{
							if (! grep(/$index5[0]/, @temp_index5))
							{
								push @temp_index5, $index5[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index6[0]/, @temp_index6))
							{
								push @temp_index6, $index6[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5> 1) && (@index6 > 1))
					{
						my @tindex5 = @index5;
						my @tindex6 = @index6;
						my @sorted_index5 = sort {$a <=> $b} @index5;
						foreach my $k (@sorted_index5)
						{
							@tindex6 = @index6;
							my @a2 = split ("\t", $matchset2[$k]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my @sorted_index6 = sort {$a <=> $b} @index6;
							foreach my $j (@sorted_index6)
							{
								my @a1 = split ("\t", $matchset1[$j]);
								$toligo1 = $a1[6];
								chomp $toligo1;
								my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
								$toligo1 = $toligos1[0];
								if (defined $toligos1[1])
								{
									my $toligo1_num = $toligos1[1];
								}
								my $tstart1 = $a1[2];
								my $tend1 = $a1[3];
								my $seqlen1 = abs($tend2 - $tstart1);
								my $seqlen2 = abs($tend1 - $tstart2);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@tindex6 = grep {!/$j/} @tindex6;
								}
								else
								{
									if (! grep(/$j/, @temp_index6)) 
									{
										push @temp_index6, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index5)) 
									{
										push @temp_index5, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@tindex6 == 0)
							{
								@tindex5 = grep {!/$k/} @tindex5;
							}
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					if ($n == $#matchset2)
					{
						if ((!@index5) && (!@index6))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset2);
						}
					}
				}
				elsif ($n == $#matchset2)
				{
					if ((!@index5) && (!@index6))
					{
						push (@no_matching_pair, @matchset1);
						push (@no_matching_pair, @matchset2);
					}
				}
			}

			@index1 = @temp_index1;
			@index2 = @temp_index2;
			@index3 = @temp_index3;
			@index4 = @temp_index4;
			@index5 = @temp_index5;
			@index6 = @temp_index6;

			if (@temp_index1 != 0)
			{
				my @temp = split ("\t", $matchset1[$temp_index1[$#temp_index1]]);
				$tempstart = $temp[2];
			}
			elsif (@temp_index2 != 0)
			{
				my @temp = split ("\t", $matchset3[$temp_index2[$#temp_index2]]);
				$tempstart = $temp[2];
			}

			if ((@temp_index1 > 1) || (@temp_index2 > 1))
			{
				foreach my $t (@temp_index1)
				{
					push (@no_matching_pair, $matchset1[$t]);
				}
				foreach my $v (@temp_index2)
				{
					push (@no_matching_pair, $matchset3[$v]);
				}
				(@index1, @index2) = ();
			}

			if (@temp_index3 != 0)
			{
				my @temp = split ("\t", $matchset1[$temp_index3[$#temp_index3]]);
				$tempstart = $temp[2];
			}
			elsif (@temp_index4 != 0)
			{
				my @temp = split ("\t", $matchset3[$temp_index4[$#temp_index4]]);
				$tempstart = $temp[2];
			}

			if ((@temp_index3 > 1) || (@temp_index4 > 1))
			{
				foreach my $t (@temp_index3)
				{
					push (@no_matching_pair, $matchset1[$t]);
				}
				foreach my $v (@temp_index4)
				{
					push (@no_matching_pair, $matchset3[$v]);
				}
				(@index3, @index4) = ();
			}

			if (@temp_index5 != 0)
			{
				my @temp = split ("\t", $matchset2[$temp_index5[$#temp_index5]]);
				$tempstart = $temp[2];
			}
			elsif (@temp_index6 != 0)
			{
				my @temp = split ("\t", $matchset1[$temp_index6[$#temp_index6]]);
				$tempstart = $temp[2];
			}

			if ((@temp_index5 > 1) || (@temp_index6 > 1))
			{
				foreach my $t (@temp_index5)
				{
					push (@no_matching_pair, $matchset2[$t]);
				}
				foreach my $v (@temp_index6)
				{
					push (@no_matching_pair, $matchset1[$v]);
				}
				(@index5, @index6) = ();
			}

			foreach my $j (@index1)
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$matchset1[$j]");
			}
			foreach my $p (@index2)
			{		
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$matchset3[$p]");
				close (OUTPUT);
			}
			if (!@index1)
			{
				foreach my $k (@index3)
				{
					open (OUTPUT2, '>>match.txt');
					print OUTPUT2 ("$matchset1[$k]");
				}
			}
			if (!@index2)
			{
				foreach my $n (@index4)
				{
					open (OUTPUT2, '>>match.txt');
					print OUTPUT2 ("$matchset3[$n]");
					close (OUTPUT2);
				}
			}						
			foreach my $j (@index5)
			{
				push (@no_matching_pair, $matchset2[$j]);
			}
			if (!@index2)
			{
				foreach my $j (@index6)
				{
					push (@no_matching_pair, $matchset1[$j]);
				}
			}
		}
		else
		{
			foreach my $item (@matchset1)
			{
				my @set1 = split ("\t", $item);
				my $m = $set1[0];
				my $s = $set1[1];
				push (@match_maid_chr_1, $m.$s);
				push (@match_chr_1, $m);
				$seen{$m} = 1;
				$seen{$m.$s} = 1;
			}	

			for (my $n = 0; $n < @matchset3; $n++)
			{
				my @set3 = split ("\t", $matchset3[$n]);
				my $m = $set3[0];
				my $s = $set3[1];
				if ($s eq "+") {$s = "-";}
				elsif ($s eq "-") {$s = "+";}
				push (@match_maid_chr_3, $m.$s);
				push (@match_chr_3, $m);
				if ($seen{$m.$s}) # 5 +  ->  5 -
				{
					@index1= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
					@index2= grep {$match_maid_chr_3[$_] eq $m.$s} 0..$#match_maid_chr_3;

					if ((@index1 == 1) && (@index2 > 1)) # more than one match
					{
						my @a1 = split ("\t", $matchset1[$index1[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];
						my @sorted_index2 = sort {$a <=> $b} @index2;
						foreach my $k (@sorted_index2)
						{
							my @a3 = split ("\t", $matchset3[$k]);
							$toligo3 = $a3[6];
							chomp $toligo3;
							my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
							$toligo3 = $toligos3[0];
							if (defined $toligos3[1])
							{
								my $toligo3_num = $toligos3[1];
							}
							my $tstart3 = $a3[2];
							my $tend3 = $a3[3];
							my $seqlen1 = abs($tend3 - $tstart1);
							my $seqlen2 = abs($tend1 - $tstart3);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index2 = grep {!/$k/} @index2;
							}
							else
							{
								if (! grep(/$k/, @temp_index2)) 
								{
									push @temp_index2, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index1[0]/, @temp_index1))
								{
									push @temp_index1, $index1[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index2 == 0)
						{
							@index1 = grep {!/$index1[0]/} @index1;
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1 > 1) && (@index2 == 1))
					{
						my @a3 = split ("\t", $matchset3[$index2[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my @sorted_index1 = sort {$a <=> $b} @index1;
						foreach my $k (@sorted_index1)
						{
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my $seqlen1 = abs($tend1 - $tstart3);
							my $seqlen2 = abs($tend3 - $tstart1);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

							if ($seqlen > 500)
							{
								@index1 = grep {!/$k/} @index1;
							}
							else
							{
								if (! grep(/$k/, @temp_index1)) 
								{
									push @temp_index1, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index2[0]/, @temp_index2))
								{
									push @temp_index2, $index2[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index1 == 0)
						{
							@index2 = grep {!/$index2[0]/} @index2;
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1== 1) && (@index2 == 1))
					{
						my @a1 = split ("\t", $matchset1[$index1[0]]);
						$toligo1 = $a1[6];
						chomp $toligo1;
						my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
						$toligo1 = $toligos1[0];
						if (defined $toligos1[1])
						{
							my $toligo1_num = $toligos1[1];
						}
						my $tstart1 = $a1[2];
						my $tend1 = $a1[3];

						my @a3 = split ("\t", $matchset3[$index2[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];

						my $seqlen1 = abs($tend3 - $tstart1);
						my $seqlen2 = abs($tend1 - $tstart3);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 500)
						{
							shift @index1;
							shift @index2;
						}
						else
						{
							if (! grep(/$index1[0]/, @temp_index1))
							{
								push @temp_index1, $index1[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index2[0]/, @temp_index2))
							{
								push @temp_index2, $index2[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
						@index1 = @temp_index1;
						@index2 = @temp_index2;
					}

					elsif ((@index1> 1) && (@index2 > 1))
					{
						my @tindex1 = @index1;
						my @tindex2 = @index2;
						my @sorted_index1 = sort {$a <=> $b} @index1;
						foreach my $k (@sorted_index1)
						{
							@tindex2 = @index2;
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my @sorted_index2 = sort {$a <=> $b} @index2;
							foreach my $j (@sorted_index2)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart1);
								my $seqlen2 = abs($tend1 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@tindex2 = grep {!/$j/} @tindex2;
								}
								else
								{
									if (! grep(/$j/, @temp_index2)) 
									{
										push @temp_index2, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index1)) 
									{
										push @temp_index1, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@tindex2 == 0)
							{
								@tindex1 = grep {!/$k/} @tindex1;
							}
						}
						 @index1 = @temp_index1;
						 @index2 = @temp_index2;
					}

					if ($n == $#matchset3)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset3);
						}
					}
				}
				else # 5 -  -> 5 -
				{
					if ($seen{$m})
					{
						@index3= grep {$match_maid_chr_1[$_] eq $m.$s} 0..$#match_maid_chr_1;
						if ($s eq "+") {$s = "-";}
						elsif ($s eq "-") {$s = "+";}
						@index4= grep {$match_maid_chr_3[$_] eq $m.$s} 0..$#match_maid_chr_3;
						my @sorted_index3 = sort {$a <=> $b} @index3;
						foreach my $k (@sorted_index3)
						{
							my @a1 = split ("\t", $matchset1[$k]);
							$toligo1 = $a1[6];
							chomp $toligo1;
							my @toligos1 = grep {$_} split /(\d+)/, $toligo1;
							$toligo1 = $toligos1[0];
							if (defined $toligos1[1])
							{
								my $toligo1_num = $toligos1[1];
							}
							my $tstart1 = $a1[2];
							my $tend1 = $a1[3];
							my @sorted_index4 = sort {$a <=> $b} @index4;
							foreach my $j (@sorted_index4)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart1);
								my $seqlen2 = abs($tend1 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@index4 = grep {!/$j/} @index4;
								}
								else
								{
									if (! grep(/$j/, @temp_index4)) 
									{
										push @temp_index4, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index3)) 
									{
										push @temp_index3, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@index4 == 0)
							{
								@index3 = grep {!/$k/} @index3;
							}
						}
						@index3 = @temp_index3;
						@index4 = @temp_index4;

						if ($n == $#matchset3)
						{
							if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
							{
								push (@no_matching_pair, @matchset1);
								push (@no_matching_pair, @matchset3);
							}
						}														
					}
					elsif ($n == $#matchset3)
					{
						if ((!@index1) && (!@index2) && (!@index3) && (!@index4))
						{
							push (@no_matching_pair, @matchset1);
							push (@no_matching_pair, @matchset3);
						}
					}
				}
			}

			@match_maid_chr_3 = ();
			@match_chr_3 = ();

			foreach my $item (@matchset2)
			{
				my @set2 = split ("\t", $item);
				my $m = $set2[0];
				my $s = $set2[1];
				push (@match_maid_chr_2, $m.$s);
				push (@match_chr_2, $m);
				$seen{$m} = 1;
				$seen{$m.$s} = 1;
			}

			for (my $n = 0; $n < @matchset3; $n++)
			{
				my @set3 = split ("\t", $matchset3[$n]);
				my $m = $set3[0];
				my $s = $set3[1];
				push (@match_maid_chr_3, $m.$s);
				push (@match_chr_3, $m);
				if ($seen{$m}) # 5   ->  5 
				{
					@index5= grep {$match_chr_2[$_] eq $m} 0..$#match_chr_2;
					@index6= grep {$match_chr_3[$_] eq $m} 0..$#match_chr_3;
		
					if ((@index5 == 1) && (@index6 > 1)) # more than one match	
					{
						my @a2 = split ("\t", $matchset2[$index5[0]]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my @sorted_index6 = sort {$a <=> $b} @index6;
						foreach my $k (@sorted_index6)
						{
							my @a3 = split ("\t", $matchset3[$k]);
							$toligo3 = $a3[6];
							chomp $toligo3;
							my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
							$toligo3 = $toligos3[0];
							if (defined $toligos3[1])
							{
								my $toligo3_num = $toligos3[1];
							}
							my $tstart3 = $a3[2];
							my $tend3 = $a3[3];
							my $seqlen1 = abs($tend3 - $tstart2);
							my $seqlen2 = abs($tend2 - $tstart3);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index6 = grep {!/$k/} @index6;
							}
							else
							{
								if (! grep(/$k/, @temp_index6)) 
								{
									push @temp_index6, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index5[0]/, @temp_index5))
								{
									push @temp_index5, $index5[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index6 == 0)
						{
							@index5 = grep {!/$index5[0]/} @index5;
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5 > 1) && (@index6 == 1))
					{
						my @a3 = split ("\t", $matchset3[$index6[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my @sorted_index5 = sort {$a <=> $b} @index5;
						foreach my $k (@sorted_index5)
						{
							my @a2 = split ("\t", $matchset2[$k]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my $seqlen1 = abs($tend2 - $tstart3);
							my $seqlen2 = abs($tend3 - $tstart2);
							my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
							if ($seqlen > 500)
							{
								@index5 = grep {!/$k/} @index5;
							}
							else
							{
								if (! grep(/$k/, @temp_index5)) 
								{
									push @temp_index5, $k;
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
								if (! grep(/$index6[0]/, @temp_index6))
								{
									push @temp_index6, $index6[0];
									if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
								}
							}
						}
						if (@index5 == 0)
						{
							@index6 = grep {!/$index6[0]/} @index6;
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5 == 1) && (@index6 == 1))
					{
						my @a2 = split ("\t", $matchset2[$index5[0]]);
						$toligo2 = $a2[6];
						chomp $toligo2;
						my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
						$toligo2 = $toligos2[0];
						if (defined $toligos2[1])
						{
							my $toligo2_num = $toligos2[1];
						}
						my $tstart2 = $a2[2];
						my $tend2 = $a2[3];
						my @a3 = split ("\t", $matchset3[$index6[0]]);
						$toligo3 = $a3[6];
						chomp $toligo3;
						my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
						$toligo3 = $toligos3[0];
						if (defined $toligos3[1])
						{
							my $toligo3_num = $toligos3[1];
						}
						my $tstart3 = $a3[2];
						my $tend3 = $a3[3];
						my $seqlen1 = abs($tend3 - $tstart2);
						my $seqlen2 = abs($tend2 - $tstart3);
						my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;

						if ($seqlen > 500)
						{
							shift @index5;
							shift @index6;
						}
						else
						{
							if (! grep(/$index5[0]/, @temp_index5))
							{
								push @temp_index5, $index5[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
							if (! grep(/$index6[0]/, @temp_index6))
							{
								push @temp_index6, $index6[0];
								if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
							}
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					elsif ((@index5> 1) && (@index6 > 1))
					{
						my @tindex5 = @index5;
						my @tindex6 = @index6;
						my @sorted_index5 = sort {$a <=> $b} @index5;
						foreach my $k (@sorted_index5)
						{
							@tindex6 = @index6;
							my @a2 = split ("\t", $matchset2[$k]);
							$toligo2 = $a2[6];
							chomp $toligo2;
							my @toligos2 = grep {$_} split /(\d+)/, $toligo2;
							$toligo2 = $toligos2[0];
							if (defined $toligos2[1])
							{
								my $toligo2_num = $toligos2[1];
							}
							my $tstart2 = $a2[2];
							my $tend2 = $a2[3];
							my @sorted_index6 = sort {$a <=> $b} @index6;
							foreach my $j (@sorted_index6)
							{
								my @a3 = split ("\t", $matchset3[$j]);
								$toligo3 = $a3[6];
								chomp $toligo3;
								my @toligos3 = grep {$_} split /(\d+)/, $toligo3;
								$toligo3 = $toligos3[0];
								if (defined $toligos3[1])
								{
									my $toligo3_num = $toligos3[1];
								}
								my $tstart3 = $a3[2];
								my $tend3 = $a3[3];
								my $seqlen1 = abs($tend3 - $tstart2);
								my $seqlen2 = abs($tend2 - $tstart3);
								my $seqlen = ($seqlen1 > $seqlen2) ? $seqlen1 : $seqlen2;
								if ($seqlen > 500)
								{
									@tindex6 = grep {!/$j/} @tindex6;
								}
								else
								{
									if (! grep(/$j/, @temp_index6)) 
									{
										push @temp_index6, $j;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
									if (! grep(/$k/, @temp_index5)) 
									{
										push @temp_index5, $k;
										if (! grep(/$m/, @chr_list))
									{ 
										push (@chr_list, $m);
									}
									}
								}
							}
							if (@tindex6 == 0)
							{
								@tindex5 = grep {!/$k/} @tindex5;
							}
						}
						@index5 = @temp_index5;
						@index6 = @temp_index6;
					}

					if ($n == $#matchset3)
					{
						if ((!@index5) && (!@index6))
						{
							push (@no_matching_pair, @matchset2);
							push (@no_matching_pair, @matchset3);
						}
					}
				}
				elsif ($n == $#matchset3)
				{
					if ((!@index5) && (!@index6))
					{
						push (@no_matching_pair, @matchset2);
						push (@no_matching_pair, @matchset3);
					}
				}
			}

			@index1 = @temp_index1;
			@index2 = @temp_index2;
			@index3 = @temp_index3;
			@index4 = @temp_index4;
			@index5 = @temp_index5;
			@index6 = @temp_index6;

			if (@temp_index1 != 0)
			{
				my @temp = split ("\t", $matchset1[$temp_index1[$#temp_index1]]);
				$tempstart = $temp[2];
			}
			elsif (@temp_index2 != 0)
			{
				my @temp = split ("\t", $matchset3[$temp_index2[$#temp_index2]]);
				$tempstart = $temp[2];
			}

			if ((@temp_index1 > 1) || (@temp_index2 > 1))
			{
				foreach my $t (@temp_index1)
				{
					push (@no_matching_pair, $matchset1[$t]);
				}
				foreach my $v (@temp_index2)
				{
					push (@no_matching_pair, $matchset3[$v]);
				}
				(@index1, @index2) = ();
			}

			if (@temp_index3 != 0)
			{
				my @temp = split ("\t", $matchset1[$temp_index3[$#temp_index3]]);
				$tempstart = $temp[2];
			}
			elsif (@temp_index4 != 0)
			{
				my @temp = split ("\t", $matchset3[$temp_index4[$#temp_index4]]);
				$tempstart = $temp[2];
			}

			if ((@temp_index3 > 1) || (@temp_index4 > 1))
			{
				foreach my $t (@temp_index3)
				{
					push (@no_matching_pair, $matchset1[$t]);
				}
				foreach my $v (@temp_index4)
				{
					push (@no_matching_pair, $matchset3[$v]);
				}
				(@index3, @index4) = ();
			}

			if (@temp_index5 != 0)
			{
				my @temp = split ("\t", $matchset2[$temp_index5[$#temp_index5]]);
				$tempstart = $temp[2];
			}
			elsif (@temp_index6 != 0)
			{
				my @temp = split ("\t", $matchset3[$temp_index6[$#temp_index6]]);
				$tempstart = $temp[2];
			}

			if ((@temp_index5 > 1) || (@temp_index6 > 1))
			{
				foreach my $t (@temp_index5)
				{
					push (@no_matching_pair, $matchset2[$t]);
				}
				foreach my $v (@temp_index6)
				{
					push (@no_matching_pair, $matchset3[$v]);
				}
				(@index5, @index6) = ();
			}

			foreach my $j (@index1)
			{
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$matchset1[$j]");
			}
			foreach my $p (@index2)
			{		
				open (OUTPUT, '>>match.txt');
				print OUTPUT ("$matchset3[$p]");
				close (OUTPUT);
			}
			if (!@index1)
			{
				foreach my $k (@index3)
				{
					open (OUTPUT2, '>>match.txt');
					print OUTPUT2 ("$matchset1[$k]");
				}
			}
			if (!@index2)
			{
				foreach my $n (@index4)
				{
					open (OUTPUT2, '>>match.txt');
					print OUTPUT2 ("$matchset3[$n]");
					close (OUTPUT2);
				}
			}						
			foreach my $j (@index5)
			{
				push (@no_matching_pair, $matchset2[$j]);
			}
			if (!@index2)
			{
				foreach my $j (@index6)
				{
					push (@no_matching_pair, $matchset3[$j]);
				}
			}
		}
	}	
}

for my $n (@list)
{	
	open (OUTPUT3, '>>maid_list.txt');
	print OUTPUT3 ("$n\n");
	close (OUTPUT3);
}
print "Total number of MAID being processed: $cnt\n";
