#!usr/bin/perl -w
use Getopt::Long;
use strict;
my ($in1,$in2,$out,$glen,$row,$help);
GetOptions
(
	"1=s"=>\$in1,
	"2=s"=>\$in2,
	"out=s"=>\$out,
	"gl=i"=>\$glen,
	"r=i"=>\$row,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
	-1 <file>	input the original file whose positions start from 1
	-2 <file>	input the modified file whose positions start from 8001
	-gl <int>	the length of your mitochondrial genome, depending on the length of reference
	-r <int>	the row containing positions' information, default 0
	-o <file>	output the merged file which should be sorted by positions after this step
INFO
die $usage if ($help || !$in1 || !$out || !$in2 || !$glen );
if(!$row){
	$row=0;
}

open OUT,">$out";

open IN1,"<$in1";
while(<IN1>){
	chomp;
	next if $_ =~ "number";
	my @temp1=split/\s+/,$_;
	if($temp1[$row] =~ /([0-9]+)/){
		if($1 >= 4001 and $1 <= 12000){
			print OUT "$_\n";
		}
	}
}
close IN1;

open IN2,"<$in2";
while(<IN2>){
	chomp;
	next if $_ =~ "number";
	my @temp2=split/\s+/,$_;
	my $s=0;
	my $line;
	if($temp2[$row] =~ /([0-9]+)/){
		if($1 >= 4001 and $1 <= ($glen-8000)){
			$s=$1+8000;
			if($temp2[$row] =~ s/$1/$s/){
				$line=join("\t",@temp2);
				print OUT "$line\n";
			}
		}
		if($1 >= ($glen-7999) and $1 <=($glen-4000)){
			$s=$1-($glen-8000);
			if($temp2[$row] =~ s/$1/$s/){
				$line=join("\t",@temp2);
				print OUT "$line\n";
			}
		}
	}
}

close IN2;
close OUT;
