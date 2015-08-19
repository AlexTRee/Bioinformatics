#!/usr/bin/perl
use strict;
use warnings;

`rm maid_oligo_list.txt`;
`rm maid_oligo_data.txt`;
`rm result.txt`;

my ($filename, $noinfo_maid, $maid, $oligo, @maid_oligo, @unsorted_noinfo_maid, @sorted_noinfo_maid, @unsorted_maid, @sorted_maid, @oligo);
$filename = 'data.txt';
open (FILE, $filename)|| die "Can't open data file!\n";
my @lines = sort <FILE>;

foreach my $line (@lines)
{
	my @arr = split ("\t", $line); 
	my $id = $arr[1];
	my $seq= $arr[2];
	my $string = substr ($seq, 0, 2);
	@maid_oligo = grep {$_} split /(\d+)/, $id;
	my $maid = $maid_oligo[0];
	chomp $maid;
	my $oligo = $maid_oligo[1];
	chomp $oligo;
	if (defined $maid_oligo[2])
	{
		$oligo = $oligo.$maid_oligo[2];
	}
	if ($string ~~ [qw(no No Ne ne Sa sa N/ n/ NA)])
	{
		@maid_oligo = grep {$_} split /(\d+)/, $id;
		push (@unsorted_noinfo_maid, $maid_oligo[0]) unless grep{$_ eq $maid_oligo[0]} @unsorted_noinfo_maid;
	}
	else
	{
		$seq =~ tr/ //ds;
		$seq =~ s/\r\n\z//;
		open (OUTPUT, '>>maid_oligo_list_temp.txt');
		print OUTPUT "$maid\t$oligo\t$seq\n";
	}
	push (@unsorted_maid, $maid) unless grep{$_ eq $maid} @unsorted_maid;
	push (@oligo, $oligo) unless grep{$_ eq $oligo} @oligo;
}
close (OUTPUT);

@sorted_noinfo_maid = sort {$a <=> $b} @unsorted_noinfo_maid;
@sorted_maid = sort {$a <=> $b} @unsorted_maid;

open (OUTPUT, '>>maid_oligo_list.txt');
print OUTPUT "MAID\n\nNo information\n";

foreach $noinfo_maid (@sorted_noinfo_maid)
{
	open (OUTPUT, '>>maid_oligo_list.txt');
	print OUTPUT $noinfo_maid."\n";
}

open (OUTPUT, '>>maid_oligo_list.txt');
print OUTPUT "\n\nWith information\n";

foreach $maid (@sorted_maid)
{
	open (OUTPUT, '>>maid_oligo_list.txt');
	print OUTPUT $maid."\n";
}

open (OUTPUT, '>>maid_oligo_list.txt');
print OUTPUT "\n\n";
@oligo = sort @oligo;

foreach $oligo (@oligo)
{
	open (OUTPUT, '>>maid_oligo_list.txt');
	print OUTPUT $oligo."\n";
	close (OUTPUT);
}

open (FILE2, "maid_oligo_list_temp.txt") || die "Can't open data file!\n";
my @lines2 = <FILE2>;
@lines2 = sort {(split(/\t/,$a))[0]<=>(split(/\t/,$b))[0] || (split(/\t/,$a))[1] cmp (split(/\t/,$b))[1]} @lines2;
close (FILE2);

foreach my $line2 (@lines2)
{
	my @arr = split ("\t", $line2); 
	my $maid = $arr[0];
	my $oligo = $arr[1];
	my $seq = $arr[2];
	open (OUTPUT, '>>maid_oligo_data.txt');
	print OUTPUT ">$maid$oligo\n$seq";
	close (OUTPUT);
}
`rm maid_oligo_list_temp.txt`;
