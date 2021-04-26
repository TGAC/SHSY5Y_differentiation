#!/net/isi-software/server/bin/perl
use strict;
use warnings;

my $project=$ARGV[0];
CREATE_LOG(\@ARGV);

my ($gtf, $genome, $outpath, $output)=@ARGV[1 .. $#ARGV];

#$gtf: annotation file in gtf format including the full path of to the file
#$genome: genome sequence in fasta format including the full path to the file
#$outpath: path to the output
#$output: name of the output file

my %cds;

PARSE_GTF($gtf);

open(OUT, ">$outpath/$output");
GET_SEQ($genome);
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
	print LOG "### Script generated to extract transcript sequences from a genome file in fasta format using a gtf file \n";
	print LOG "$0\n";
	for my $n (1 .. $#{$array}){
		print LOG join ("\t", $array->[$n]),"\n";
	}
	close LOG;
}
#########################################################
#########################################################
sub PARSE_GTF{
	my $file=shift;
	open(IN, "$file")||die"IN $file\n";
	while(<IN>){
		if($_=~/\texon\t/){
			chomp;
			$_=~s/[\s]+/\t/g;
			$_=~s/;//g;
			$_=~s/"//g;
			my @split=split /\t/, $_;
			push @{$cds{$split[0]}{$split[11]}}, [$split[3], $split[4], $split[6]];
		}
	}
	close IN;
	
}	
#########################################################
#########################################################
sub GET_SEQ{
	my ($gene)=(@_);
	my $seq="";
	my $head;
	my $c=0;
	open(SEQ, "$gene")||die"SEQ $gene\n";
	while(<SEQ>){
		if($c==0){
			if($_=~/>/){
				chomp;
				$_=~s/[\s]+/\t/g;
				my @split=split /\t/, $_;
				$split[0]=~s/>//;
				$head=$split[0];
				$c=1;
				$_="";
			}
		}
		if($c==1){
			unless($_=~/>/){
				$_=~s/\s+//g;
				$seq.=$_;
			}
			if($_=~/>/ or eof(SEQ)){
				if($cds{$head}){
					foreach my $loc (keys %{$cds{$head}}){
						@{$cds{$head}{$loc}}=sort{$a->[0]<=>$b->[0]||$a->[1]<=>$b->[1]}@{$cds{$head}{$loc}};
						my $temp_cds;
						for my $i (0 .. $#{$cds{$head}{$loc}}){
							$temp_cds.=substr($seq, $cds{$head}{$loc}[$i][0]-1, ($cds{$head}{$loc}[$i][1]-$cds{$head}{$loc}[$i][0])+1);
						}
						if($cds{$head}{$loc}[0][2]=~/-/){
							$temp_cds=~tr/ATGCNatgcn/TACGNtacgn/;
							my @bases=split //, $temp_cds;
							@bases=reverse @bases;
							$temp_cds=join("", @bases);
							@bases=();
						}
						print OUT join ("\n", ">$loc", $temp_cds),"\n";
					}
				}
				$head="";
				$seq="";
				if($_=~/>/){
					chomp;
					$_=~s/[\s]+/\t/g;
					my @split=split /\t/, $_;
					$split[0]=~s/>//;
					$head=$split[0];
					$c=1;
					$_="";
				}
				else{
					$c=0;
				}
			}
		}
	}
	close SEQ;
}
