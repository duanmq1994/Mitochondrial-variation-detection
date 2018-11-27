#!usr/bin/perl -w
#count the base-type and indel-type file together
use Getopt::Long;
use strict;
my($in,$out1,$out2,$out3,$cut,$HF,$glen,$help);
GetOptions
(
	"i=s"=>\$in,
	"ob=s"=>\$out1,
	"oi=s"=>\$out2,
	"ov=s"=>\$out3,
	"c=i"=>\$cut,
	"HF=f"=>\$HF,
	"g=i"=>\$glen,
	"help|?"=>\$help
);
my $usage=<<INFO;
Usage:
	perl $0	[options]
Options:
	-i <file>		input .sam/.bam file after alignment and sort
	-ob	<file>	output base-type count file
	-oi	<file>	output indel-type count file
	-c <int>	Phred quality value cutoff of every base, default=20
	-HF <float>	the threshold value of MAF(minor allele frequency), defalut with no threshold
	-ov <file>	output the heteroplasmic point variations according to the -HF parament, default no output
	-g <int>	length of mitochondrial genome
INFO

die "You lose some parameters!\n" if ($help || !$in || !$out1 ||!$out2|| !$glen);

if($in =~ /bam/){
	open (IN,"samtools view $in |");
}
else{
	open IN,"<$in";
}

$cut=20 if !$cut;
if($out3){
	die "You should assign a threshold value of MAF to -HF!\n" if !$HF;
}else{
	die "You should output the heteroplasmic point variations into a file!\n" if $HF;
}

open OUT1,">$out1";
open OUT2,">$out2";
open OUT3,">$out3" if $out3;

my%indelcount=();
my @count;
for(my $a=0;$a<=$glen;$a++){
	for(my $b=1;$b<=6;$b++){
		$count[$a]->[$b]=0;
	}
}

my($aloci,$cloci,$num_insertion,$num_deletion)=(0,0,0,0);
while(<IN>){
	chomp;
	my @temp=split/\s+/,$_;
	next if $#temp<17 or $temp[2] eq "*" or $temp[5] eq "*";
	my @seq=split//,$temp[9];
	my @qua=split//,$temp[10];
	$aloci=$temp[3];
		my @c=split//,$temp[5];
		my @md=split//,$temp[17];
		my $e=-1;my @num=();
		for(my $i=0;$i<=$#c;$i++){
			if($c[$i] =~ /[0-9]/){
				push(@num,$c[$i]);
				next;
			}
			if($c[$i] =~ /[A-Z]/){
				my $n=join('',@num);
				@num=();
				if ($c[$i] ne "D"){
					my $s=$e+1;$e+=$n;
					if($c[$i] eq "M"){
						for my$m($s..$e){
							if(ord($qua[$m]) >= ($cut+33)){
								$cloci=$aloci+($m-$s);
								$count[$cloci]->[1]++ if $seq[$m] eq "A";
								$count[$cloci]->[2]++ if $seq[$m] eq "T";
								$count[$cloci]->[3]++ if $seq[$m] eq "C";
								$count[$cloci]->[4]++ if $seq[$m] eq "G";
							}
						}
						$aloci=$aloci+$n-1;
						next;
					}
					if($c[$i] eq "I"){
						my $in;
						for my$q($s..$e){
							if(ord($qua[$q]) < ($cut+33)){
								$in = "null";
								last;
							}
							$in=join("_",$aloci,"I",join('',@seq[$s..$e]));
							$count[$aloci]->[5]++;
							if($in ne "null"){
								$indelcount{$in}++;
							}
						}
						$aloci+=1;
						next;
					}
				}
				else{
					$aloci+=1;
					for(my$k=0;$k<=$#md;$k++){
						my $de;
						if($md[$k] eq "^"){
							my $de_s=$k+1;
							my $de_e=$k+$n;
							$de="null";
							if(ord($qua[$e]) > ($cut+33) and ord($qua[$e+1]) > ($cut+33)){
								$de=join("_",$aloci,"D",join('',@md[$de_s..$de_e]));
								$count[$aloci]->[6]++;
							}
							if($de ne "null" and $de !~ /3107/){
								$indelcount{$de}++;
							}
							$aloci+=$n;
							$md[$k]="N";
							last;
						}
						else{
							next;
						}
					}
				}
			}
		}
}

print OUT1 "pos\tref_allele\tnumber_of_ref_allele\t1st_alteration\tnumber_of_the_1st_alteration\t2nd_alteration\tnumber_of_the_2nd_alteration\t3rd_alteration\tnumber_of_the_3rd_alteration\tInsertion\tDeletion\tDepthi\tMAF\n";
print OUT3 "pos\tref_allele\tnumber_of_ref_allele\t1st_alteration\tnumber_of_the_1st_alteration\t2nd_alteration\tnumber_of_the_2nd_alteration\t3rd_alteration\tnumber_of_the_3rd_alteration\tInsertion\tDeletion\tDepthi\tMAF\n" if $out3;
for(my$i=1;$i<=$glen;$i++){
	print OUT1 "$i\t";
	my %hash=(
		"A"=>$count[$i]->[1],
		"T"=>$count[$i]->[2],
		"C"=>$count[$i]->[3],
		"G"=>$count[$i]->[4]
	);
	my $depth=0;
	my @key=sort{$hash{$b} <=> $hash{$a}}keys%hash;
	my $minor_base=$key[2];
	my $minor_depth=$hash{$minor_base};
	for (@key){
			print OUT1 "$_\t$hash{$_}\t";
			print OUT3 "$_\t$hash{$_}\t" if $out3;
			$depth+=$hash{$_};
	}
	$depth+=$count[$i]->[6];
	my $MAF=sprintf "%.2f",($minor_depth/$depth);
	my $insertion=$count[$i]->[5];
	my $deletion=$count[$i]->[6];
	print OUT1 "$insertion\t$deletion\t$depth\t$MAF\n";
	print OUT3 "$insertion\t$deletion\t$depth\t$MAF\n" if $out3 and $MAF >= $HF;		
}

my @indels=keys%indelcount;
for my$indel (@indels){
	print OUT2 "$indelcount{$indel}\t$indel\n";
}

close IN;
close OUT2 if $out2;
close OUT1;
close OUT3 if $out3;
