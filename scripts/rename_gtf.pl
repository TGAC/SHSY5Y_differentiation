#!/net/isi-software/server/bin/perl
use strict;
use warnings;

my $project=$ARGV[0];

my ($reference, $talon, $outpath, $output)=@ARGV[1 .. $#ARGV];

my %genes;
PARSE_GTF($reference);

my %transcripts;
PARSE_TALON($talon);

my %flag;
foreach my $trans (keys %transcripts){
	@{$transcripts{$trans}}=sort{$b->[3]<=>$a->[3]}@{$transcripts{$trans}};
	@{$flag{$trans}}=($transcripts{$trans}[0][0], $transcripts{$trans}[0][1], $transcripts{$trans}[0][2]);
}

open(OUT, ">$outpath/$output");
open (IN, "$talon")||die"IN $talon\n";
while(<IN>){
	chomp;
	if($_=~/TALON/){
		my @var=split /\t/, $_;
		my $line=$var[8];
		$_=~s/;//g;
		$_=~s/\"//g;
		$_=~s/[\s]+/\t/g;
		my @split=split /\t/, $_;
		if($flag{$split[11]}){
			my $name1;
			my $name2;
			for my $i (8 .. $#split){
				if($split[$i]=~/gene_id/){
					$name1=$split[$i+1];
				}
				if($split[$i]=~/gene_name/){
					$name2=$split[$i+1];
				}
			}
			$line=~s/$name1/$flag{$split[11]}[0]/;
			$line=~s/$name2/$flag{$split[11]}[1]/;

			my $strand2=1;
			if($var[6]=~/-/){
				$strand2=-1;
			}
			if($flag{$split[11]}[2]==$strand2){
				print OUT join ("\t", @var[0 .. 6], ".", $line),"\n";
			}
			else{
				if($strand2==1){
					print OUT join ("\t", @var[0 .. 5], "-", ".", $line),"\n";
				}
				else{
					print OUT join ("\t", @var[0 .. 5], "+", ".", $line),"\n";
				}
			}
		}
		else{
			print OUT join ("\t", @var[0 .. 6], ".", $line),"\n";
		}
	}
	else{
		print OUT "$_\n";
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
	print LOG "#### Script writen to rename TALON transcripts that match an existing genome model from a reference gtf file\n";
	print LOG "$0\n";
	for my $n (1 .. $#{$array}){
		print LOG join ("\t", $array->[$n]),"\n";
	}
	close LOG;
}

#########################################################
#########################################################
sub PARSE_GTF{
	my ($file)=(@_);
	open(IN, "$file")||die"IN $file\n";
	while(<IN>){
		chomp;
		if($_=~/\tgene\t/){
			$_=~s/;//g;
			$_=~s/\"//g;
			$_=~s/[\s]+/\t/g;
			my @split=split /\t/, $_;
			my $strand=1;
			if($split[6]=~/-/){
				$strand=-1;
			}
			push @{$genes{$split[0]}}, [$split[3], $split[4], $strand, $split[9], $split[13]];
		}
	}
	close IN;
}


#########################################################
#compare Talon genes models to reference gene models to enable reassignments of transcripts based on overlap with existing annotations
#########################################################
sub PARSE_TALON{
	my ($file)=(@_);
	my %temp;
	open(IN, "$file")||die"IN $file\n";
	while(<IN>){
		chomp;
		if($_=~/\ttranscript\t/ && $_=~/TALON/){
			$_=~s/;//g;
			$_=~s/\"//g;
			$_=~s/[\s]+/\t/g;
			my @split=split /\t/, $_;
			push @{$temp{$split[0]}}, [$split[3], $split[4], $split[11]];
		}
	}
	close IN;
	
	foreach my $ele (keys %temp){
		@{$genes{$ele}}=sort{$a->[0]<=>$b->[0]||$a->[1]<=>$b->[1]}@{$genes{$ele}};
		@{$temp{$ele}}=sort{$a->[0]<=>$b->[0]||$a->[1]<=>$b->[1]}@{$temp{$ele}};
		for my $i (0 .. $#{$temp{$ele}}){
			my $s1=$temp{$ele}[$i][0];
			my $e1=$temp{$ele}[$i][1];
			for my $j (0 .. $#{$genes{$ele}}){
				my $s2=$genes{$ele}[$j][0];
				my $e2=$genes{$ele}[$j][1];
				if($s2>$e1){
					last;
				}
				elsif(($s1>=$s2 && $e1<=$e2)||($s1<=$s2 && $e1>=$e2)||($s1>=$s2 && $s1<$e2 && $e1>=$e2)||($s1<=$s2 && $s2<$e1 && $e1<=$e2)){
					my $dif=0;
					if($s1>=$s2 && $e1<=$e2){
						$dif=1;
					}
					elsif($s1<=$s2 && $e1>=$e2){
						$dif=(($e2-$s2)+1)/(($e1-$s1)+1);
					}
					elsif($s1>=$s2 && $s1<$e2 && $e1>=$e2){
						$dif=(($e2-$s1)+1)/(($e1-$s1)+1);
					}
					elsif($s1<=$s2 && $s2<$e1 && $e1<=$e2){
						$dif=(($e1-$s2)+1)/(($e1-$s1)+1);
					}
					push @{$transcripts{$temp{$ele}[$i][2]}}, [$genes{$ele}[$j][3], $genes{$ele}[$j][4], $genes{$ele}[$j][2], $dif];
				}
			}
		}
	}
}

