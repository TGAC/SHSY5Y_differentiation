#!/net/isi-software/server/bin/perl
use strict;
use warnings;

my $project=$ARGV[0];
CREATE_LOG(\@ARGV);

my ($gtf1, $gtf2, $outpath, $output)=@ARGV[1 .. $#ARGV];

#$gtf1: TALON annotation file in gtf format including the full path to it
#$gtf2: reference annotation file in gtf format including the full path to it
#$outpath: path to the output file
#$output: output file name


my %genes;
my %transcripts;
open(OUT, ">$outpath/$output");
open(IN, "$gtf1")||die"IN $gtf1\n";
while(<IN>){
	chomp;
	my $line=$_;
	if($_=~/\ttranscript\t/){
		$_=~s/;//g;
		$_=~s/\"//g;
		$_=~s/[\s]+/\t/g;
		my @split=split /\t/, $_;
		$genes{$split[9]}=1;
		$transcripts{$split[11]}=1;
	}
	print OUT "$line\n";
}
close IN;

open(IN, "$gtf2")||die"IN $gtf2\n";
while(<IN>){
	chomp;
	if($_!~/\#/){
		my $line=$_;
		if($_=~/\tgene\t/){
			$_=~s/;//g;
			$_=~s/\"//g;
			$_=~s/[\s]+/\t/g;
			my @split=split /\t/, $_;
			if(!$genes{$split[9]}){
				print OUT "$line\n";
			}
		}
		else{
			$_=~s/;//g;
			$_=~s/\"//g;
			$_=~s/[\s]+/\t/g;
			my @split=split /\t/, $_;
			if(!$transcripts{$split[11]}){
				print OUT "$line\n";
			}
		}
	}
}
close IN;

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
	print LOG "### Script written to compare TALON gtf file after reassignment of transcripts to an existing reference annotation and add the novel genes and transcripts to the reference annotation to produce a complete annotation as not all genes are reported duing the TALON process \n";
	print LOG "$0\n";
	for my $n (1 .. $#{$array}){
		print LOG join ("\t", $array->[$n]),"\n";
	}
	close LOG;
}

