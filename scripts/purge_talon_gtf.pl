#!/net/isi-software/server/bin/perl
use strict;
use warnings;

my $project=$ARGV[0];
CREATE_LOG(\@ARGV);

my ($gtf, $outpath, $output)=@ARGV[1 .. $#ARGV];

# $gtf: TALON annotation file in gtf format including te full path
# $outpath: path to the output
# $output: file name of the output file


my %flagged;
my %strand;
open(IN, "$gtf")||die"IN $gtf\n";
while(<IN>){
	chomp;
	if($_=~/TALON/){
		my $line=$_;
		$_=~s/\"//g;
		$_=~s/;//g;
		$_=~s/[\s]+/\t/g;
		my @split=split /\t/, $_;
		if($line=~/antisense_gene "TRUE"/){
			$flagged{$split[9]}=$line;
		}
	}
	if($_=~/\tgene\t/ && $_!~/TALON/){
		$_=~s/\"//g;
		$_=~s/;//g;
		$_=~s/[\s]+/\t/g;
		my @split=split /\t/, $_;
		if($split[6]=~/-/){
			$strand{$split[9]}=-1;
		}
		else{
			$strand{$split[9]}=1;
		}
	}
}

open(OUT, ">$outpath/$output");
open(IN, "$gtf")||die"IN $gtf\n";
while(<IN>){
	chomp;
	if($_=~/TALON/){
		my $line=$_;
		$_=~s/\"//g;
		$_=~s/;//g;
		$_=~s/[\s]+/\t/g;
		my @split=split /\t/, $_;
		if(!$flagged{$split[9]}){
			if($strand{$split[9]}){
				if(($strand{$split[9]}==-1 && $split[6]=~/-/) || $strand{$split[9]}==1 && $split[6]=~/\+/){
					print OUT "$line\n";
				}
			}
			else{
				print OUT "$line\n";
			}
		}
	}
	else{
		print OUT "$_\n";
	}
}
close IN;
close OUT;

print join ("\t", "TOTAL NUMBER TALON GENES FLAGGED:", scalar keys %flagged),"\n";

open(OUT, ">$outpath/TALONG_gene_flagged_overlapping.gtf");
foreach my $loc (keys %flagged){
	print OUT "$flagged{$loc}\n";
}
close OUT;


#########################################################
#########################################################
sub CREATE_LOG{
	my ($array)=(@_);
	my $date=localtime();
	if(-e "$project"){
		open ("IN", "$project")||die"IN $project\n";
		my @project_log=<IN>;
		close IN;
		open(LOG, ">$project");
		for my $n (0 .. $#project_log){
			print LOG $project_log[$n];
		}
	}
	else{
		open(LOG, ">$project");
	}
	print LOG "\n\n#########################################################\n#########################################################\n
	$date\n
#########################################################\n#########################################################\n
	\n";
	print LOG "Script written to parse a TALON gtf output and remove novel loci that are perfect antisense mapping\n";
	print LOG "$0\n";
	for my $n (1 .. $#{$array}){
		print LOG join ("\t", $array->[$n]),"\n";
	}
	close LOG;
}

