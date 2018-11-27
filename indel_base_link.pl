#!usr/bin/perl -w
#print the heteroplasmic point mutations and indels linkage
use Getopt::Long;
use strict;
my($in1,$in2,$cut,$help);
GetOptions
(
	"1=s"=>\$in1,
	"2=s"=>\$in2,
	"c=i"=>\$cut,
	"help|?"=>\$help
);
my $usage=<<INFO;
Usage:
	perl $0	[options]
Options:
	-1 <file>	input the merged point-variation.tsv file
	-2 <file>	input .sam/.bam mapping file, should be sorted
	-c <int>	Phred quality score cutoff value of each base, defalut=20 
INFO

die $usage if ($help || !$in1 || !$in2  || !$cut);

$cut=20 if !$cut;

open IN1,"<$in1";
my @snp_pos;

while(<IN1>){
	chomp;
	next if $_ =~ "pos";
	my @temp1=split/\s+/,$_;
	my $mutations=join("\t",@temp1[1..$#temp1]);
	@snp_pos[$temp1[0]]=$mutations;
}

if($in2 =~ /bam/){
	open (IN2,"samtools view $in2 |");
}
else{
	open IN2,"<$in2";
}

open OUT,">id-variation.tsv";
my @seq;
my @qua;
my @num;

while(<IN2>){
	chomp;
	my @temp2=split/\s+/,$_;
	next if $#temp2<4 or $temp2[2] eq "*" or $temp2[5] eq "*";
	@seq=split//,$temp2[9];
	@qua=split//,$temp2[10];
	my $aloci=$temp2[3];
	my @m_seq=();
	my @m_qua=();
	if($temp2[5] ne "120M"){
		my @c=split//,$temp2[5];
		my @md=split//,$temp2[17];
		my $e=-1;
		for(my$i=0;$i<=$#c;$i++){
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
						push(@m_seq,@seq[$s..$e]);
						push(@m_qua,@qua[$s..$e]);
						$aloci=$aloci+$n-1;
						next;
					}
					if($c[$i] eq "I"){
						my $insertion;
						for my $q($s..$e){
							if(ord($qua[$q]) < ($cut+33)){
								$insertion = "null";
								last;
							}
							$insertion=join("_",$aloci,"I",join('',@seq[$s..$e]));
						}
						if($insertion ne "null"){
							print OUT "$temp2[0]\t$insertion\n";
						}
						$aloci+=1;
						next;
					}
				}
				else{
					$aloci+=1;
					for(my $k=0;$k<=$#md;$k++){
						my $deletion;
						if($md[$k] eq "^"){
							my $de_s=$k+1;
							my $de_e=$k+$n;
							$deletion="null";
							if(ord($qua[$e]) > ($cut+33) and ord($qua[$e+1]) > ($cut+33)){
								$deletion=join("_",$aloci,"D",join('',@md[$de_s..$de_e]));
							}
							if($deletion ne "null" and $deletion !~ /3107/){
								my $N="N" x ($de_e-$de_s+1);
								push(@m_seq,$N);
								push(@m_qua,$N);
								print OUT "$temp2[0]\t$deletion\n";
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
	else{
		@m_seq=@seq;
		@m_qua=@qua;
	}
	for(my$i=0;$i<=$#m_seq;$i++){
		my $locus=$temp2[3]+$i;
		if(defined($snp_pos[$locus])){
			my @genotypes=split/\s+/,$snp_pos[$locus];		
			for my$base(@genotypes){
				if($m_seq[$i] =~ $base and ord($m_qua[$i])>=($cut+33)){
					print OUT "$temp2[0]\t$locus$base\n";
				}
			}
		}
	}
}

close IN1;
close IN2;
close OUT;
