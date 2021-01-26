#!/usr/bin/env perl
# ==============================================================
# Tomas Bruna, Alex Lomsadze
# Copyright 2019, Georgia Institute of Technology, USA
#
# This script runs large splicing alignment jobs on cluster using PBS.
# Script should be adjusted for cluster specific configuration.
# ===================================

use strict;
use warnings;

use Getopt::Long;
use File::Spec;
use Cwd qw( abs_path cwd );
use File::Temp qw/ tempfile /;
use FindBin qw( $RealBin );
use Data::Dumper;

# ------------------------------------------------
my $v = 0;
my $debug = 0;
my $cfg;
my $log;
my $bin = $RealBin;
my $work_dir = cwd;
# ------------------------------------------------

my $N  = 90;            # number of jobs
my $K  = 8;             # number of cores per job
my $db = '';            # path to blast formatted database ; plain FASTA formatted file with sequence and blast db name is expected at this location
my $aligner= '';        # Alignment engine
my $seq_file = '';      # input sequence
my $list_file = '';     # file with pairs of gene_id and protein_id (blast output)
my $out_file_gff = '';  # output in GFF
my $min_exon_score = 25;

my $node_dir = '/tmp';  # use this directory on for holding tmp data on node
my $SPALN_OUT = "spaln.gff";
# ------------------------------------------------

Usage() if $#ARGV == -1;
ParseCMD();
my $store = $aligner; # put files with sequence and outputs here
CheckBeforeRun();

print STDERR "[" . localtime() . "] Starting spliced alignment on PBS\n" if $v;
print STDERR "[" . localtime() . "] Splitting the input file\n" if $v;
my @files = PrepareSequences();
my %pairs = PreparePairs( \@files );
my %jobs;   # keep job id's here

foreach my $file ( @files )
{
	print STDERR "[" . localtime() . "] submitting $file\n" if $v;
	RunOnPBS($file, $pairs{$file});
}

print STDERR "[" . localtime() . "] All jobs submitted. Waiting for jobs to finish.\n" if $v;
WaitForIt(1);
print STDERR "[" . localtime() . "] Finished, writing output\n" if $v;
CreateOutput();

print STDERR Dumper(\%jobs) if $debug;
print STDERR "[" . localtime() . "] Spliced alignment on PBS finished\n" if $v;
exit 0;


# ================== sub =====================
sub RunOnPBS
{

	my ($name, $list) = @_;

 	my $pbs = UniqTmpFile();
	my $pbs_log = $pbs.".log";
	my $label = $name;

	my $text = "#!/bin/bash
#PBS -N $label
#PBS -o $pbs_log
#PBS -j oe
#PBS -l nodes=1:ppn=$K
#PBS -l walltime=20:00:00

dir=`mktemp  -p  $node_dir  -d`
cd \$dir
mkdir bin
mkdir dependencies
cp $bin/run_spliced_alignment.pl bin
cp $bin/gff_from_region_to_contig.pl bin
cp $bin/../dependencies/blosum62.csv dependencies";

	my $spalnPbs = "
cp -r $bin/../dependencies/spaln_table dependencies
cp $bin/../dependencies/spaln dependencies
cp $bin/spaln_to_gff.py bin
cp $bin/../dependencies/spaln_boundary_scorer dependencies
cp $bin/spalnBatch.sh bin

./bin/run_spliced_alignment.pl --nuc $name --prot $db --list $list --cores $K --aligner spaln --min_exon_score $min_exon_score

mv $SPALN_OUT ${name}.gff
cd ..
rm -r \$dir
";

	$text .= $spalnPbs;

	open( OUT, ">$pbs" )||die;
	print OUT $text;
	close OUT;
	chmod  0755, $pbs;


	my $id = `qsub $pbs`;

	chomp $id;
	print STDERR "qsub ID: $id\n" if $debug;

	$jobs{ $id }{'seq'} = $name;
	$jobs{ $id }{'list'} = $list;
	$jobs{ $id }{'out_gff'} = "${name}.gff";
	$jobs{ $id }{'script'} = $pbs;
	$jobs{ $id }{'log'} = $pbs_log;
};
# ------------------------------------------------
sub CreateOutput
{
	print STDERR "creating output\n" if $debug;
	my $OUT_GFF;

	open( $OUT_GFF, ">$SPALN_OUT" ) or die "error on open file: $SPALN_OUT $!\n";

	foreach my $id ( sort keys %jobs )
	{
		print STDERR "[" . localtime() . "] Collecting results for $jobs{$id}{'seq'}\n" if $v;

		open( my $IN_GFF, $jobs{$id}{'out_gff'} ) or die "error on open file: $jobs{$id}{'out_gff'} $!\n";
		while( <$IN_GFF> )
		{
			print $OUT_GFF $_;
		}
		close $IN_GFF;

		# clean tmp files
		if ( -s $jobs{$id}{'log'} )
		{
			print STDERR "warning, check this: $jobs{$id}{'log'}\n";
		}
		else
		{
			unlink $jobs{ $id }{'seq'};
			unlink $jobs{ $id }{'list'};
			unlink $jobs{ $id }{'out_gff'};
			unlink $jobs{ $id }{'script'};
			unlink $jobs{ $id }{'log'};
		}
	}

	close $OUT_GFF;
};
# ------------------------------------------------
sub WaitForIt
{
	my $status = shift;

	print STDERR "waiting for spaln\n" if $debug;

	if ( !%jobs )
	{
		print STDERR "error, no registered jobs\n";
		exit 1;
	}

	while( $status )
	{
		sleep 5;
		$status = 0;

		foreach my $id ( keys %jobs )
		{
			if ( ! -e $jobs{$id}{'log'} )  { $status = 1; }
		}
	}
};
# ------------------------------------------------
sub UniqTmpFile
{
	my ( $fh, $tmp_name ) = tempfile( "pbs_XXXXX" );
	if ( !fileno($fh) ) { die "Can't open temporally  file: $!\n"; }
	close $fh;
	chmod 0755, $fh;
	return $tmp_name;
};
# ------------------------------------------------
sub PreparePairs
{
	my $ref = shift;

	my %h;

	# load list
	my %gene_to_proteins;
	LoadPairs( $list_file, \%gene_to_proteins );

	# load deflines from each files
	foreach my $name ( @{$ref} )
	{
		my @arr = ();
        LoadDeflines( $name, \@arr );

        # save pair info
		SavePairs( $name ."_id", \@arr, \%gene_to_proteins );

		$h{$name} = $name ."_id";
	}

	return %h;
};
# ------------------------------------------------
sub SavePairs
{
	my( $name, $ref_arr, $ref_hash ) = @_;

	print STDERR "[" . localtime() . "] saving id: $name\n" if $v;

	open( my $OUT, ">$name" ) || die "$! on open $name\n";
	foreach my $id ( @{$ref_arr} )
	{
		if ( exists $ref_hash->{$id} )
		{
			print $OUT $ref_hash->{$id};
		}
		else
		{
			print STDERR "warning, no pairs for: $name\n" if $v;
		}
	}
	close $OUT;
};
# ------------------------------------------------
sub LoadDeflines
{
	my( $name, $ref ) = @_;

	print STDERR "[" . localtime() . "] reading defline: $name\n" if $v;

	open( my $IN, $name ) || die "$! on open $name\n";
	while( my $line = <$IN> )
	{
		if( $line !~ /^>/ ) {next;}

		if( $line =~ /^>\s*(\S+)\s*/ )
		{
			push @$ref, $1;
		}
	}
	close $IN;
};
# ------------------------------------------------
sub LoadPairs
{
	my( $name, $ref ) = @_;
	open( my $IN, $name ) || die "$! on open $name\n";
	while( my $line = <$IN> )
	{
		if( $line =~ /^\s*$/ ) {next;}
		if( $line =~ /^\s*#/ ) {next;}

		if( $line =~ /(\S+)\s+(\S+)/ )
		{
			$ref->{ $1 } .= $line;
		}
		else { print STDERR "error, unexpected file format found$ 0: $line\n"; exit 1; }
	}
	close $IN;
}
# ------------------------------------------------
sub PrepareSequences
{
	foreach my $f ( ReadFileNames( abs_path($store) ) ) { unlink $f; }

	chdir $store;
	system(" $bin/../dependencies/probuild  --seq $seq_file  --split_numfile $N  --split_fasta n ") and die("error on probuild");
	chdir $work_dir;

	return ReadFileNames( abs_path($store) );
};
# ------------------------------------------------
sub ReadFileNames
{
	my $dir = shift;
	my @list;

	opendir( DIR, $dir ) or die "error on open directory $0: $dir, $!\n";
	foreach my $file ( grep{ /\S+_\d+$/ } readdir(DIR) )
	{
		$file =  File::Spec->catfile( $dir, $file );
		if ( -f $file )
		{
			$file = abs_path($file);
			push @list, $file;
		}
		else { print STDERR "error, unexpected name found $0: $file\n"; exit 1; }
	}
	closedir DIR;

	my $message = scalar @list ." files in directory";
	print STDERR "$message\n" if $v;

        my @sorted_list =
                map{ $_->[0] }
                sort { $a->[1] <=> $b->[1] }
                map { [$_, $_=~/(\d+)/] }
                        @list;

	return @sorted_list;
};
# ------------------------------------------------
sub CheckBeforeRun
{
	print STDERR "check before run\n" if $debug;

	$bin      = ResolvePath( $bin );
	$work_dir = ResolvePath( $work_dir );
	$seq_file = ResolvePath( $seq_file );
	$list_file = ResolvePath( $list_file );

	if( !$seq_file ) { print STDERR "error, file is missing $0:  option --seq\n"; exit 1; }
	if( !$list_file ) { print STDERR "error, file is missing $0:  option --list\n"; exit 1; }

	if( !$db ) { print STDERR "error, database name is missing $0:  option --db\n"; exit 1; }

	$aligner = lc $aligner;

	if ( !$aligner ) {
		print STDERR "Reuqired option --aligner is missing.\n"; exit 1;
	}

	if ( $aligner ne "spaln" ) {
		if ( $aligner eq "prosplign" ) {
				print STDERR "error, ProSplign is not supported in the current version of this script. A version with ProSplign support ";
				print STDERR "is available in https://github.com/gatech-genemark/ProtHint/releases/tag/v2.4.0.\n";
				exit 1;
			} else {
				print STDERR "error, invalid aligner specified: $aligner. In the current version of this script, only \"Spaln\" is supported.\n";
				exit 1;
			}
	}

	mkdir $store;
	if( ! -e $store ) { print STDERR "error, directory not found: $store\n"; exit 1; }
	$store = ResolvePath($store);
	$cfg->{'d'}->{'store'} = $store;

	if ( !$node_dir ) { print STDERR "error missing tmp directory name on nodes\n"; exit 1; }
};
# ------------------------------------------------
sub ResolvePath
{
	my( $name, $path ) = @_;
	return '' if !$name;
	$name = File::Spec->catfile( $path, $name ) if ( defined $path and $path );
	if( ! -e $name ) { print STDERR "error, file not found $0: $name\n"; exit 1; }
	return abs_path( $name );
};
# ------------------------------------------------
sub ParseCMD
{
	print STDERR "parse cmd\n" if $debug;

	my $cmd = $0;
	foreach my $str (@ARGV) { $cmd .= ( ' '. $str ); }

	my $opt_results = GetOptions
	(
		'seq=s'   => \$seq_file,
		'list=s'  => \$list_file,
		'db=s'    => \$db,
		'tmp=s'   => \$node_dir,
		'N=i'     => \$N,
		'K=i'     => \$K,
		'aligner=s' => \$aligner,
		'verbose' => \$v,
		'debug'   => \$debug,
		'min_exon_score=f' => \$min_exon_score
	);

	if( !$opt_results ) { print STDERR "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print STDERR "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
	$v = 1 if $debug;

	# save information for debug
	$cfg->{'d'}->{'seq_file'}  = $seq_file;
	$cfg->{'d'}->{'list_file'} = $list_file;
	$cfg->{'d'}->{'tmp'}   = $node_dir;
	$cfg->{'d'}->{'bd'}    = $db;
	$cfg->{'d'}->{'N'}     = $N;
	$cfg->{'d'}->{'K'}     = $K;
	$cfg->{'d'}->{'aligner'} = $aligner;
	$cfg->{'d'}->{'v'}     = $v;
	$cfg->{'d'}->{'debug'} = $debug;
	$cfg->{'d'}->{'cmd'}   = $cmd;
	$cfg->{'d'}->{'min_exon_score'} = $min_exon_score;

	print STDERR Dumper($cfg) if $debug;
};
# ------------------------------------------------
sub Usage
{
	my $txt =
"# ----------
Usage: $0 options

run spliced alignment on PBS

  --seq  [s] FASTA file with nucleotide sequences
  --list [s] list with  dna-to-protein mappings pair names
  --aligner [s]  Which spliced alignment tool to use. Current version only supports \"Spaln\"
  --db   [s] proteins are here
  --tmp  [s] name of temporary directory on node
  --N    [i] number of jobs to submit
  --K    [i] number of cores per job (on the same node)
  --v        verbose
  --min_exon_score [f] discard all hints inside/neighboring exons with score lower than minExonScore. Spaln specific option.
  --debug
# -----------
";
	print $txt;
	exit 1;
};
# --------------------------------------------
