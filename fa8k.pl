#!usr/perl/bin -w
use strict;
use Getopt::Long;
my ($in,$help);

GetOptions(
		"in=s"=>\$in,
		"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]

	This script will creat a new .fasta file starting from the 8001 position of mitochondrial genome, to eliminate the influence on the head and tail reads because of a circular reference.

Options:
	-i <file>	input .fasta file of mitochondrial genome
INFO
die $usage if ($help || !$in);
open IN,"<$in";
open OUT,">MT-8k.fasta";

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
	
