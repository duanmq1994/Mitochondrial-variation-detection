#!usr/perl/bin -w
#output all the reads containing one or more hetero-loci(s)
use strict;
use Getopt::Long;
my ($in,$help,$out);

GetOptions(
	"i=s"=>\$in,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
	-i <file>	input "sequence_id\tvariations" file which has been sorted by sequence_id
INFO
die $usage if ($help || !$in ||!$out);

open IN,"<$in";
open OUT,">linkages-count.tsv";
print OUT "number_of_linkages\tall_linkages\n";

my $readp="null";
my @hetero_everyid;
my %hetero_linkage;
my $hetero_merged;
my %count;
my @uniq_hetero;
while(<IN>){
	chomp;
	my @temp=split/\s+/,$_;
	next if($temp[1] =~ /N/);
	if(($temp[0] ne $readp) and ($readp ne "null")){
		%count=();
		@uniq_hetero = grep { ++$count{ $_ } < 2; } @hetero_everyid;
		if($#uniq_hetero >= 1){	
			$hetero_merged=join"\t",@uniq_hetero;
			$hetero_linkage{$hetero_merged}++;
		}
		@hetero_everyid=();
	}
	push(@hetero_everyid,$temp[1]);
	$readp=$temp[0];
}	

%count=();
@uniq_hetero = grep { ++$count{ $_ } < 2; } @hetero_everyid;
if($#uniq_hetero >= 1){
	$hetero_merged=join"\t",@uniq_hetero;
	$hetero_linkage{$hetero_merged}++;
}

my @keys=keys%hetero_linkage;
for my$i(@keys){
	print OUT "$hetero_linkage{$i}\t$i\n";
}

close IN;
close OUT;
