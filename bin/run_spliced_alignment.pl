#!/usr/bin/env perl
# ==============================================================
# Tomas Bruna, Alex Lomsadze
# Copyright 2019, Georgia Institute of Technology, USA
# 
# This script takes as input protein and nucleotide sequences,
# runs spliced alignment and scores the output with boundary scorer
# ==============================================================

use strict;
use warnings;

use Getopt::Long qw( GetOptions );
use FindBin qw( $RealBin );
use File::Spec;
use Cwd qw( abs_path cwd );
use Data::Dumper;
use YAML;
use threads;
use Thread::Queue;
use MCE::Mutex;
use File::Temp qw/ tempfile /;

# ------------------------------------------------
my $v = 0;
my $debug = 0;

my $cfg;
my $log;

my $bin = $RealBin;
my $work_dir = cwd;
# ------------------------------------------------
my $nuc_file  = '';
my $prot_file = '';
my $list_file = '';
my $out_file_regions  = '';
my $out_file_asn_regions  = '';
my $out_file_asn_gff_regions  = '';
my $cores = 1;
my $aligner = '';
my $tmp_dir = '';
my $min_exon_score = 25;
# ------------------------------------------------
my $PROSPLIGN_INTRONS_OUT = "scored_introns.gff";
my $PROSPLIGN_OUT = "prosplign.gff";
my $SPALN_OUT = "spaln.gff";
my $BATCH_SIZE = 100;
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

print STDERR "[" . localtime() . "] Starting spliced alignment with $aligner\n" if $v;

if ($aligner eq "spaln") {
	$out_file_regions = "spaln.regions.gff";
} elsif ($aligner eq "prosplign") {
	$out_file_regions = "scored_introns.regions.gff";
	$out_file_asn_regions = "prosplign.regions.asn";
	$out_file_asn_gff_regions = "prosplign.regions.gff";
}

# It is important to create threads before calling ReadSequence as
# threads have a copy of all variables
my $mutex = MCE::Mutex->new;
my $q = Thread::Queue->new();

# Do not create too many temporary files
$q->limit = $cores * 2;
for (1..$cores)
{
	threads->create(\&alignerThread);
}


my %nuc;
my %prot;
my @list;

print STDERR "[" . localtime() . "] Loading alignment pairs into memory\n" if $v;

ReadList( $list_file, \@list, \%nuc, \%prot );

ReadSequence( $nuc_file,  \%nuc);
ReadSequence( $prot_file, \%prot);

if ( $debug )
{
	foreach my $val  (@list)
	{
		if ( ! exists $nuc{ $val->[0] } ) { print STDERR "$val->[0] missing\n"; }
		if ( ! exists $prot{ $val->[1] } ) { print STDERR "$val->[1] missing\n"; }
	}
}

# Reset output files
my $OUT;
open( $OUT, ">$out_file_regions" ) or die( "$!, error on open file $out_file_regions" );
close $OUT;
if ($aligner eq "prosplign") {
	open( $OUT, ">$out_file_asn_regions" ) or die( "$!, error on open file $out_file_asn_regions" );
	close $OUT;
	open( $OUT, ">$out_file_asn_gff_regions" ) or die( "$!, error on open file $out_file_asn_gff_regions" );
	close $OUT;
}

my $tmp_prot_file;
my $tmp_nuc_file;
my $tmp_out_file;

my $nuc_id;
my $prot_id;
my $ref;

my $counter = 0;
my $batchCounter = 1;
my $currentBatchCounter = 0;
my $BATCH;
my $nextPrint = 0;
my $startTime = time();

$ENV{ALN_TAB} = "$bin/../dependencies/spaln_table";

my $pairsCount = (scalar @list);
print STDERR "[" . localtime() . "] Pairs loaded. Number of pairs to align: $pairsCount\n" if $v;
print STDERR "[" . localtime() . "] Starting the alignments\n" if $v;
foreach $ref ( @list )
{

	if ($currentBatchCounter == 0) {
		open($BATCH, ">batch_$batchCounter") or die("$!, error on open file batch_$batchCounter");
	}

	printProgress();

	$counter += 1;
	$currentBatchCounter += 1;

	$nuc_id = $ref->[0];
	$prot_id = $ref->[1];

	$tmp_nuc_file  = "nuc_" . $counter;
	$tmp_prot_file = "prot_" . $counter;
	$tmp_out_file  = "out_" . $counter;

	SaveSingleFasta( $tmp_nuc_file,  $nuc_id, $nuc{$nuc_id} );
	SaveSingleFasta( $tmp_prot_file, $prot_id, $prot{$prot_id} );
	addToBatch($tmp_nuc_file, $tmp_prot_file, $BATCH);

	if ($currentBatchCounter == $BATCH_SIZE) {
		close $BATCH;
		$q->enqueue("batch_$batchCounter");
		$batchCounter += 1;
		$currentBatchCounter = 0;
	}
}

# Enqueue the last remaining batch
if ($currentBatchCounter != 0) {
	close $BATCH;
	$q->enqueue("batch_$batchCounter");
}

# Finish and wait for threads
$q->end();
foreach my $thr (threads->list())
{
	$thr -> join();
}

print STDERR "[" . localtime() . "] $pairsCount/$pairsCount (100%) pairs aligned\n" if $v;
print STDERR "[" . localtime() . "] Alignment of pairs finished\n" if $v;
print STDERR "[" . localtime() . "] Translating coordinates from local pair level to contig level\n" if $v;

# Cerate final output with correct coordinates
if ($aligner eq "spaln") {
	system("$bin/gff_from_region_to_contig.pl --in_gff $out_file_regions --seq $nuc_file --out_gff $SPALN_OUT");
} elsif ($aligner eq "prosplign") {
	system("$bin/gff_from_region_to_contig.pl --in_gff $out_file_regions --seq $nuc_file --out_gff $PROSPLIGN_INTRONS_OUT");
	system("$bin/gff_from_region_to_contig.pl --in_gff $out_file_asn_gff_regions --seq $nuc_file --out_gff $PROSPLIGN_OUT");
	unlink $out_file_asn_gff_regions;
}

unlink $out_file_regions;

print STDERR "[" . localtime() . "] Finished spliced alignment\n" if $v;
exit 0;

#------------------------------------------------
# A single alignment, processed by a single thread
sub alignerThread
{
	while (my $batchFile = $q->dequeue()) {
		if ($aligner eq "spaln") {
			$tmp_out_file = "${batchFile}_out";
			system("$bin/spalnBatch.sh $batchFile $tmp_out_file $min_exon_score");
		} elsif (($aligner eq "prosplign")) {
			alignWithProSplign($tmp_nuc_file, $tmp_prot_file, $tmp_out_file);
			$mutex->lock;
			AppendToFile ($out_file_asn_regions, "${tmp_out_file}_asn");
			AppendToFile ($out_file_asn_gff_regions, "${tmp_out_file}_asn_gff");
			$mutex->unlock;
			unlink "${tmp_out_file}_asn";
			unlink "${tmp_out_file}_asn_gff";
		}

		# Safe write
		$mutex->lock;
		AppendToFile ($out_file_regions, $tmp_out_file);
		$mutex->unlock;

		unlink $tmp_out_file;
		unlink $batchFile;
	}
}

sub alignWithProSplign
{
	my $tmp_nuc_file = shift;
	my $tmp_prot_file = shift;
	my $tmp_out_file = shift;

	# -o     Output file in asn format
	# -eo    Output file with full alignment
	# -nfa   Single nucleotide sequence to read from a FASTA file
	# -pfa   Single protein sequence to read from a FASTA file
	# -nogenbank   Do not use GenBank data loader.
	# -gaps_inf    Show all gaps in info output (if not set, show frameshifts only)
	# -two_stages

	system("$bin/../dependencies/prosplign -o \"${tmp_out_file}_asn\" -eo \"${tmp_out_file}_ali\"  -nfa \"$tmp_nuc_file\" -pfa \"$tmp_prot_file\" -nogenbank -gaps_inf -two_stages");
	# Parse and score introns from alignment file
	system("$bin/../dependencies/prosplign_intron_scorer -i \"${tmp_out_file}_ali\" -o \"$tmp_out_file\" -w 10 -a -s $bin/../dependencies/blosum62.csv");
	# Convert asn output to gff
	system("$bin/asn_to_gff.pl --asn \"${tmp_out_file}_asn\" --out \"${tmp_out_file}_asn_gff\" --exons");

	unlink "${tmp_out_file}_ali";
}

#------------------------------------------------
# Parallel safe, platform independent append
sub AppendToFile	
{
	my $out_file = shift;
	my $tmp_out_file = shift;
	open(my $out_fh, '>>', $out_file) or die "Could not open '$out_file' - $!";

	if (open my $in, '<', "$tmp_out_file") {
		while (my $line = <$in>) {
			print $out_fh $line;
		}
		close $in;
	} else {
		warn "Could not open '$tmp_out_file' for reading\n";
	}

	close $out_fh;
}
#------------------------------------------------
sub UniqTmpFile
{
	my $name = shift;
	$name .= "_XXXXX";
	my ( $fh, $tmp_name ) = tempfile( $name );
	if ( !fileno($fh) ) { die "Can't open temporally  file: $!\n"; }
	close $fh;
	return $tmp_name;
}
#------------------------------------------------
sub SaveSingleFasta
{
	my( $name, $id, $s ) = @_;
	
	open( my $OUT, ">$name" ) or die( "$!, error on open file $name" );
	print $OUT ( ">". $id ."\n");
	print $OUT ( $s  ."\n");
	close $OUT;	
}
#------------------------------------------------
sub addToBatch
{
	my($tmp_nuc_file, $tmp_prot_file, $BATCH) = @_;
	print $BATCH ($tmp_nuc_file . "\t" . $tmp_prot_file  . "\n");
}
#------------------------------------------------
sub ReadList
{
	my( $name, $ref, $h_f, $h_s ) = @_;
	
	open( my $IN, "$name" ) || die "$! on open $name\n";
	while( my $line = <$IN> )
	{
		if( $line =~ /^\s*$/ ) {next;}
		if( $line =~ /^\s*#/ ) {next;}
		
		if( $line =~ /(\S+)\s+(\S+)/ )
		{
			push @{ $ref }, [$1, $2];
			
			$h_f->{$1} = '';
			$h_s->{$2} = '';
		}
		else { print STDERR "error, unexpected file format found$0: $line\n"; exit 1; }
	}
	close $IN;
}
#------------------------------------------------
sub ReadSequence
{
	my ($name, $ref) = @_;
	open( my $IN, "$name" ) || die "$! on open $name\n";

	my $id = "";
	
	my $skip = 0;

	while( my $line = <$IN> )
	{
		if( $line =~ /^>(\S+)\s*/ )
		{
			$id = $1;
			
			$skip = 0;
			if ( ! exists $ref->{$id} ) {
				$skip = 1;
			}
		}
		else
		{
			if( !$id ) { print STDERR "error, fasta record whithout definition line found $0: $line\n"; exit 1; }

			if ( $skip ) {next;}

			$line = uc $line;

			# remove non alphabet
			$line =~ tr/A-Z//dc;

			$ref->{ $id } .= $line;
		}
	}
 
	close $IN;
};
# ------------------------------------------------
sub printProgress
{
	my $permille = int(($counter * 1000) / $pairsCount);
	if ($permille >= $nextPrint) {
		printf STDERR "[" . localtime() . "] Enqueueing pair $counter/$pairsCount (%.1f%%)", $permille / 10 if $v;
		$nextPrint = $permille + 1;
		my $elapsedTime = time() - $startTime;
		if ($counter < $BATCH_SIZE * $cores * 8) {
			print STDERR "\n" if $v;
		} else {
			my $secondsPerPermille = $elapsedTime / $permille;
			my $secondsLeft = int((1000 - $permille)  * $secondsPerPermille) + 1;
			printf STDERR ". Est. time left: %02d:%02d:%02d (hh:mm:ss)\n",
				int($secondsLeft / 3600),
				($secondsLeft / 60) % 60,
				$secondsLeft % 60 if $v;
		}
	}
}
# ------------------------------------------------
sub CheckBeforeRun
{
	print STDERR "check before run\n" if $debug;
		
	$bin      = ResolvePath( $bin );
	$work_dir = ResolvePath( $work_dir );
	
	$nuc_file  = ResolvePath( $nuc_file );
	$prot_file = ResolvePath( $prot_file );
	$list_file = ResolvePath( $list_file );
	
	if( !$nuc_file )  { print STDERR "error, required file name is missing $0:  option --nuc\n"; exit 1; }
	if( !$prot_file ) { print STDERR "error, required file name is missing $0:  option --prot\n"; exit 1; }
	if( !$list_file ) { print STDERR "error, required file name is missing $0:  option --list\n"; exit 1; }
	
	if( $cores < 1 ) { print STDERR "error, out of range prot_ids specified for number of cores $0: $cores\n"; exit 1; }

	$aligner = lc $aligner;

	if ( !$aligner ) {
		print STDERR "Reuqired option --aligner is missing. Please specify the aligner. Valid options are: \"Spaln\", \"ProSplign\".\n"; exit 1;
	} 

	if ( $aligner ne "spaln" && $aligner ne "prosplign" ) 
	{
		print STDERR "error, invalid aligner specified: $aligner. Valid options are: \"Spaln\", \"ProSplign\".\n";
		exit 1;
	}
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
		'nuc=s'     => \$nuc_file,
		'prot=s'    => \$prot_file,
		'list=s'    => \$list_file,
		'cores=i'   => \$cores,
		'aligner=s' => \$aligner,
		'verbose'   => \$v,
		'debug'     => \$debug,
		'min_exon_score=f' => \$min_exon_score
	);

	if( !$opt_results ) { print STDERR "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print STDERR "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
	$v = 1 if $debug;

	# save informaton for debug
	$cfg->{'d'}->{'nuc_file'}  = $nuc_file;
	$cfg->{'d'}->{'prot_file'} = $prot_file;
	$cfg->{'d'}->{'list_file'} = $list_file;
	$cfg->{'d'}->{'cores'}     = $cores;
	$cfg->{'d'}->{'aligner'}   = $aligner;
	$cfg->{'d'}->{'v'}     = $v;
	$cfg->{'d'}->{'debug'} = $debug;
	$cfg->{'d'}->{'cmd'}   = $cmd;
	$cfg->{'d'}->{'min_exon_score'} = $min_exon_score;
	
	print STDERR Dumper($cfg) if $debug;
};
# ------------------------------------------------
sub Usage
{
	print qq(# -------------------
Usage:  $0

Required options:
  --nuc       [name] name of file with nucleotide sequences
  --prot      [name] name of file with protein sequences
  --list      [name] list of nuc to prot mapping
  --aligner   [name] Which spliced alignment tool to use. Valid options are: \"Spaln, ProSplign\"


 Optional parameters:
  --cores            [number] number of threads to use
  --min_exon_score   [number] discard all hints inside/neighboring exons with score lower than minExonScore. Spaln specific option.

Developer options:
  --verbose
  --debug
# -------------------
);
	exit 1;
}
# ================== END sub =====================
