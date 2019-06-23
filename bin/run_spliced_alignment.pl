#!/usr/bin/perl
# ==============================================================
# Tomas Bruna, Alex Lomsadze
# Copyright 2019, Georgia Institute of Technology, USA
# 
# This script takes as input protein and nucleotide sequences 
# and runs Spliced alignment
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
# ------------------------------------------------
my $PROSPLIGN_INTRONS_OUT = "scored_introns.gff";
my $PROSPLIGN_OUT = "prosplign.gff";
my $SPALN_OUT = "spaln.gff";
my $SPALN_MIN_EXON_SCORE = 50;
my $SPALN_MIN_START_SCORE = 100;
my $SPALN_MIN_STOP_SCORE = 100;
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

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
$q->limit = $cores * 4;
for (1..$cores)
{
	threads->create(\&alignerThread);
}


my %nuc;
my %prot;
my @list;

ReadList( $list_file, \@list, \%nuc, \%prot );

print "list size: ". (scalar @list) ."\n" if $v;

ReadSequence( $nuc_file,  \%nuc,  \%nuc );
ReadSequence( $prot_file, \%prot, \%prot );

if ( $debug )
{
	foreach my $val  (@list)
	{
		if ( ! exists $nuc{ $val->[0] } ) { print "$val->[0] missing\n"; }
		if ( ! exists $prot{ $val->[1] } ) { print "$val->[1] missing\n"; }
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


my $key;
my $value;
my $ref;

my $counter = 0;

$ENV{ALN_TAB} = "$bin/../dependencies/spaln_table";

print "looping list\n" if $v;
foreach $ref ( @list )
{

	$key = $ref->[0];
	$value = $ref->[1];

	$counter += 1;

	# Some fasta protein headers may contain a slash
	my $file_value = $value;
	$file_value =~ tr/\//_/;
	# Header may be longer than a file name length limit
	$file_value = substr $file_value, 0, 200;

	print $key ." ". $value ." ". $counter  ."\n" if $v;

	$tmp_nuc_file  = "nuc_". $key ."_". $counter;
	$tmp_prot_file = "prot_". $file_value ."_". $counter;
	$tmp_out_file  = "out_". $key ."_". $file_value ."_". $counter;

	# Spaln does not like dots in file names; orthoDB headers contain dots
	$tmp_nuc_file  =~ tr/\./_/;
	$tmp_prot_file  =~ tr/\./_/;

	if ( $debug )
	{
		print "error on $tmp_nuc_file\n" if( ! -e $tmp_nuc_file );
	}

	SaveSingleFasta( $tmp_nuc_file,  $key,   $nuc{$key} );
	SaveSingleFasta( $tmp_prot_file, $value, $prot{$value} );

	my @files = ($tmp_nuc_file, $tmp_prot_file, $tmp_out_file);
	$q->enqueue(\@files);


}

print "total runs: $counter\n" if $v;

# Finish and wait for threads
$q->end();
foreach my $thr (threads->list())
{
	$thr -> join();
}


# Cerate final output with correct coordinates
if ($aligner eq "spaln") {
	system("$bin/gff_from_region_to_contig.pl --in_gff $out_file_regions --seq $nuc_file --out_gff $SPALN_OUT");
} elsif ($aligner eq "prosplign") {
	system("$bin/gff_from_region_to_contig.pl --in_gff $out_file_regions --seq $nuc_file --out_gff $PROSPLIGN_INTRONS_OUT");
	system("$bin/gff_from_region_to_contig.pl --in_gff $out_file_asn_gff_regions --seq $nuc_file --out_gff $PROSPLIGN_OUT");
	unlink $out_file_asn_gff_regions;
}

unlink $out_file_regions;

exit 0;

#------------------------------------------------
# A single alignment, processed by a single thread
sub alignerThread
{
	while (my $item = $q->dequeue()) {
		
		my $tmp_nuc_file = $$item[0];
		my $tmp_prot_file = $$item[1];
		my $tmp_out_file = $$item[2];

		if ($aligner eq "spaln") {
			alignWithSpaln($tmp_nuc_file, $tmp_prot_file, $tmp_out_file);
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

		unlink $tmp_nuc_file;
		unlink $tmp_prot_file;
		unlink $tmp_out_file;
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
	system("$bin/prosplign_parser/prosplign_parser -i \"${tmp_out_file}_ali\" -o \"$tmp_out_file\" -w 10 -a -s $bin/prosplign_parser/blosum62.csv");
	# Convert asn output to gff
	system("$bin/asn_to_gff.pl --asn \"${tmp_out_file}_asn\" --out \"${tmp_out_file}_asn_gff\" --exons");

	unlink "${tmp_out_file}_ali";
}


sub alignWithSpaln
{
	my $tmp_nuc_file = shift;
	my $tmp_prot_file = shift;
	my $tmp_out_file = shift;

	# Estimate the maximum possible possible length of the alignment, including gaps.
	my $alignmentLength = ((stat "$tmp_nuc_file")[7]) * 2;

	# -Q3    Algorithm runs in the fast heuristic mod
	# -pw    Report result even if alignment score is below threshold value
	# -S1    Dna is in the forward orientation
	# -LS    Smith-Waterman-type local alignment. This option may prune out weakly matched terminal regions.
	# -O1    Output alignment
	# -l     Number of characters per line in alignment

	# Align
	system("$bin/../dependencies/spaln -Q3 -LS -pw -S1 -O1 -l $alignmentLength  \"$tmp_nuc_file\" \"$tmp_prot_file\" 2> /dev/null > ${tmp_out_file}_ali");
	# Parse and score hints
	system("$bin/spaln_parser/spaln_parser -i \"${tmp_out_file}_ali\" -o \"$tmp_out_file\" -w 10 -a -s $bin/spaln_parser/blosum62.csv");

	unlink "${tmp_out_file}_ali";

	# system("$bin/../dependencies/spaln -Q7 -pw -S1 -O4  \"$tmp_nuc_file\" \"$tmp_prot_file\" 2> /dev/null | " .
	#	"$bin/spaln_to_gff.py > \"$tmp_out_file\" --intronScore $SPALN_MIN_EXON_SCORE " .
	#	"--startScore $SPALN_MIN_START_SCORE --stopScore $SPALN_MIN_STOP_SCORE " .
	#	"--gene \"$tmp_nuc_file\" --prot \"$tmp_prot_file\"");
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
			
			if ( defined $h_f )
			{
				$h_f->{$1} = '';
				$h_s->{$2} = '';
			}
		}
		else { print "error, unexpected file format found$0: $line\n"; exit 1; }
	}
	close $IN;
}
#------------------------------------------------
sub ReadSequence
{
	my ($name, $ref, $h_ref ) = @_;
	open( my $IN, "$name" ) || die "$! on open $name\n";

	my $id = "";
	
	my $skip = 0;

	while( my $line = <$IN> )
	{
		if( $line =~ /^>(\S+)\s*/ )
		{
			$id = $1;
			
			$skip = 0;
			if ( defined $h_ref )
			{
				$skip = 1 if ( ! exists $ref->{$id} );
			}
		}
		else
		{
			if( !$id ) { print "error, fasta record whithout definition line found $0: $line\n"; exit 1; }

			if ( defined $h_ref and $skip ) {next;}

			$line = uc $line;

			# remove non alphabet
			$line =~ tr/A-Z//dc;

			$ref->{ $id } .= $line;
		}
	}
 
	close $IN;
};
# ------------------------------------------------
sub CheckBeforeRun
{
	print "check before run\n" if $debug;
		
	$bin      = ResolvePath( $bin );
	$work_dir = ResolvePath( $work_dir );
	
	$nuc_file  = ResolvePath( $nuc_file );
	$prot_file = ResolvePath( $prot_file );
	$list_file = ResolvePath( $list_file );
	
	if( !$nuc_file )  { print "error, required file name is missing $0:  option --nuc\n"; exit 1; }
	if( !$prot_file ) { print "error, required file name is missing $0:  option --prot\n"; exit 1; }
	if( !$list_file ) { print "error, required file name is missing $0:  option --list\n"; exit 1; }
	
	if( $cores < 1 ) { print "error, out of range values specified for number of cores $0: $cores\n"; exit 1; }

	$aligner = lc $aligner;

	if ( !$aligner ) {
		print "Reuqired option --aligner is missing. Please specify the aligner. Valid options are: \"Spaln\", \"ProSplign\".\n"; exit 1; 
	} 

	if ( $aligner ne "spaln" && $aligner ne "prosplign" ) 
	{
		print "error, invalid aligner specified: $aligner. Valid options are: \"Spaln\", \"ProSplign\".\n"; 
		exit 1;
	}
};
# ------------------------------------------------
sub ResolvePath
{
	my( $name, $path ) = @_;
	return '' if !$name;
	$name = File::Spec->catfile( $path, $name ) if ( defined $path and $path );
	if( ! -e $name ) { print "error, file not found $0: $name\n"; exit 1; }
	return abs_path( $name );
};
# ------------------------------------------------
sub ParseCMD
{
	print "parse cmd\n" if $debug;
	
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
	);

	if( !$opt_results ) { print "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
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
	
	print Dumper($cfg) if $debug;
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
  --cores  [number]  number of threads to use

Developer options:
  --verbose
  --debug
# -------------------
);
	exit 1;
}
# ================== END sub =====================
