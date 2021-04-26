#!/net/isi-software/server/bin/perl
use strict;
use warnings;

my $project=$ARGV[0];
CREATE_LOG(\@ARGV);

my ($script_path, $talon_gtf, $reference_gtf, $genome, $outpath)=@ARGV[1 .. $#ARGV];

# $script_path: path to the folder with the scripts
# $talon_gtf: TALON output in gtf format including the full path
# $reference_gtf: Reference annotation in gtf format including the full path
# $genome: genome file in fasta format including the full path
# $outpath: output directory


system("perl $script_path/purge_talon_gtf.pl $project $talon_gtf $outpath talon_purged.gtf \n")==0||die"FAIL purge_talon_gtf.pl\n";

system("perl $script_path/rename_gtf.pl $project $reference_gtf $outpath/talon_purged.gtf $outpath talon_purged_renamed.gtf \n")==0||die"FAIL purge_talon_gtf_renamed.pl\n";

system ("perl $script_path/compare_complete_gtf.pl $project $outpath/talon_purged_renamed.gtf $reference_gtf $outpath talon_purged_renamed_complete.gtf \n")==0||die"FAIL talon_purged_renamed_complete.gtf\n";

system("perl $script_path/extract_transcripts_sequences_gtf.pl $project $outpath/talon_purged_renamed_complete.gtf $genome $outpath talon_transcripts.fa \n")==0||die"FAIL talon_purged_renamed_complete.gtf \n";


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
	print LOG "$0\n";
	for my $n (1 .. $#{$array}){
		print LOG join ("\t", $array->[$n]),"\n";
	}
	close LOG;
}

