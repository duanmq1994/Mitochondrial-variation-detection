#!usr/perl/bin -w
use strict;
use Getopt::Long;
my ($in,$out,$help);

GetOptions(
		"in=s"=>\$in,
		"out=s"=>\$out,
		"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
	-i <file>	input .fasta file of mitochondrial genome
	-o <file>	output modified .fasta file which statrs from 8001
INFO
die $usage if ($help || !$in || !$out);
open IN,"<$in";
open OUT,">$out";

my @gene=();
while(<IN>){
	chomp;
	if ($_ =~ />/){
		print OUT "$_\n";
		next;
	}
	else{
		my @temp=split("",$_);
		push(@gene,@temp);
	}
}
my $i=1;
my $count;
for $count(8000..$#gene){
	print OUT $gene[$count];
	print OUT "\n" if $i % 70 == 0;
	$i++;
}
my $j=1;
for $count(0..7999){
    print OUT $gene[$count];
	print OUT "\n" if ($i-1+$j) % 70 == 0;
	$j++;
}

close IN;
close OUT;	
	
