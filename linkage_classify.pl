#!usr/bin/perl -w
use Getopt::Long;
use strict;
my ($in,$out1,$out2,$out3,$help);

GetOptions(
	"in=s"=>\$in,
	"o1=s"=>\$out1,
	"o2=s"=>\$out2,
	"o3=s"=>\$out3,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
	-i <file>	input heteroplasmic variations linkage file
	-o1 <file>	output indels linkage file
	-o2 <file>	output indels-point variations linkage file
	-o3 <file>	output point variations linkage file
INFO
die $usage if ($help || !$in || !$out1 || !$out2 || !$out3);

open IN,"<$in";
open OUT1,">$out1";
open OUT2,">$out2";
open OUT3,">$out3";

print OUT1 "number_of_linkages\tindels_linkage\n";
print OUT2 "number_of_linkages\tindel-point_variations_linkage\n";
print OUT3 "number_of_linkages\tpoint_variations_linkage\n";

my %indel;
my %indel_link;
my %snp_link;
while(<IN>){
	chomp;
	next if $_ =~ "linkage";
	my @temp=split/\s+/,$_;
	my @indel_temp=();
	my @snp_temp=();
	my $ind=0;
	my $snp=0;
	for(my $i=1;$i<=$#temp;$i++){
		if($temp[$i] =~ /\_/){
			$ind++;
			push(@indel_temp,$temp[$i]);
		}
		else{
			$snp++;
			push(@snp_temp,$temp[$i]);
		}
	}
	if($ind > 1){
		my $indel_line=join("\t",@indel_temp);
		$indel_link{$indel_line}+=$temp[0];
	}
	elsif($ind < $#temp and $ind != 0){
		print OUT2 "$_\n";
	}
	if($snp > 1){
		my $snp_line=join("\t",@snp_temp);
		$snp_link{$snp_line}+=$temp[0];
	}
}

foreach my $indel_links(keys %indel_link){
	my $link_num=$indel_link{$indel_links};
	print OUT1 "$link_num\t$indel_links\n";
}

foreach my $snp_links(keys %snp_link){
	my $snp_lnum=$snp_link{$snp_links};
	print OUT3 "$snp_lnum\t$snp_links\n";
}

close IN;
close OUT1;
close OUT2;
close OUT3;
