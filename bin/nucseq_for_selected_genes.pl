#!/usr/bin/env perl
# ==============================================================
# Alex Lomsadze
# Copyright 2019, Georgia Institute of Technology, USA
# 
# This script takes as input nucleotide sequence, gene coordinates and list of gene ID's
# and outputs nucleotide sequence of specified genes with margins
# ==============================================================

use strict;
use warnings;

use Getopt::Long qw( GetOptions );
use FindBin qw( $RealBin );
use File::Spec;
use Cwd qw( abs_path cwd );
use Data::Dumper;
use YAML;

# ------------------------------------------------
my $v = 0;
my $debug = 0;

my $cfg;
my $log;

my $bin = $RealBin;
my $work_dir = cwd;
# ------------------------------------------------
my $seq_file  = '';
my $gene_stat = '';
my $gene_list = '';
my $out_file  = '';

my $margin = 2000;

my $soft_mask = 0;
my $max_gene  = 0;

# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

my %sequence;
ReadSequence( $seq_file, \%sequence );

my $seq_info;
my $gene_info;

( $seq_info, $gene_info ) = YAML::LoadFile( $gene_stat );

my %list;
ReadList( $gene_list, \%list );

foreach my $id (keys %list)
{
	if( !exists $gene_info->{$id} ) { print "warning, gene id not found in gene stat file $0: $id\n"; next; }

	my $strand = $gene_info->{$id}{'strand'};
	my $seq_id = $gene_info->{$id}{'seq_id'};

	my( $L, $R ) = GetPositions( $gene_info->{$id}, $margin, $seq_info->{$seq_id} );

	if ( !$L or !$R )
	{
		$list{$id}{'status'} = 0;
		next;
	}
	
	if ( $max_gene and ( $R - $L + 1 > $max_gene ) )
	{
		$list{$id}{'status'} = 0;
		next;
	}

	$list{$id}{'status'} = 1;

	$list{$id}{'defline'} = ">$id $list{$id}{'counts'} $L $R $strand $seq_id";
	
	if( $strand eq '+' )
	{
		$list{$id}{'s'} = substr( $sequence{$seq_id}, $L - 1, $R - $L + 1 );
	}
	elsif( $strand eq '-' )
	{
		$list{$id}{'s'} = RevComp( substr( $sequence{$seq_id}, $L - 1, $R - $L + 1 ) );
	}
	else {die;}
}

PrintOut( $out_file, \%list );

exit 0;

#------------------------------------------------
sub PrintOut
{
	my( $name, $ref ) = @_;
	
	print "printing output: $name\n" if $v;
	
	open( my $OUT, ">$name" ) or die( "$!, error on open file $name" );
	foreach my $id (keys %{$ref})
	{	
		if ( $ref->{$id}{'status'} )
		{ 
			print $OUT ($ref->{$id}{'defline'} ."\n");
			print $OUT ($ref->{$id}{'s'} ."\n\n");
		}
	}
	close $OUT;
}
#------------------------------------------------
sub RevComp
{
	my $s = shift;
	$s =~ s/T/1/g;
	$s =~ s/C/2/g;
	$s =~ s/A/T/g;
	$s =~ s/G/C/g;
	$s =~ s/1/A/g;
	$s =~ s/2/G/g;
	$s = reverse($s);
	return $s;
};
#------------------------------------------------
sub GetPositions
{
	my( $ref, $m, $max ) = @_;

	my $L = $ref->{'L'};
	my $R = $ref->{'R'};

	if ( !$R or !$L )
	{
		print "error, unexpected format:\n";
		print Dumper($ref);
		exit 1;
	}
		
	if ( $R + $m > $max )
	{
		$R = $max;
	}
	else
	{
		$R = $R + $m;
	}
		
	if ( $L - $m < 1 )
	{
		$L = 1;
	}
	else
	{
		$L = $L - $m;
	}
	
	if ( $L > $R )
	{
		print "error L is more than R:  $L, $R\n";
		print Dumper($ref);
		exit 1;
	}
	
	return ($L, $R);
}
#------------------------------------------------
sub ReadList
{
	my( $name, $ref ) = @_;
	
	print "loading list: $name\n" if $v;
	
	open( my $IN, $name ) || die "$! on open $name\n";	
	while( my $line = <$IN> )
	{
		if( $line =~ /^\s*$/ ) {next;}
		if( $line =~ /^\s*#/ ) {next;}
		
		if( $line =~ /(\S+)\s+(\S+)/ )
		{
			$ref->{$1}{'counts'} += 1;
		}
		else { print "error, unexpected file format found $0: $line\n"; exit 1; }
	}
	close $IN;
	
	print "loading list done\n" if $v;
}
#------------------------------------------------
sub ReadSequence
{
	my ($name, $ref ) = @_;
	
	print "loading sequence: $name\n" if $v;
	
	open( my $IN, $name ) || die "$! on open $name\n";

	my $id = "";

	while( my $line = <$IN> )
	{
		if( $line =~ /^>(\S+)\s*/ )
		{
			$id = $1;
		}
		else
		{
			if( $id eq '' ) { print "error, fasta record whithout definition line found $0: $line\n"; exit 1; }

			if ( $soft_mask )
			{
				$line =~ tr/[atcg]/N/;
			}

			$line = uc $line;

			# remove non alphabet
			$line =~ tr/A-Z//dc;

			# replace allowed nucleic acid code (non A T C G) by N
			$line =~ tr/RYKMSWBDHV/N/;

			# stop if unexpected character
			if( $line =~ m/[^ATCGN]/ )
				{ print "error, unexpected letter foind in sequence $0: $line\n"; exit 1; }

			$ref->{ $id } .= $line;
		}
	}
 
	close $IN;
	
	if( $debug )
	{
		print "id  length\n";
		foreach my $id (keys %{$ref} )
		{
			print ( $id ." ". length($ref->{$id}) ."\n" );
		}
	}
};
# ------------------------------------------------
sub CheckBeforeRun
{
	print "check before run\n" if $debug;
	
	$bin      = ResolvePath( $bin );
	$work_dir = ResolvePath( $work_dir );
	
	$seq_file  = ResolvePath( $seq_file );
	$gene_stat = ResolvePath( $gene_stat );
	$gene_list = ResolvePath( $gene_list );
	
	if( !$seq_file )  { print "error, required file name is missing $0:  option --seq\n"; exit 1; }
	if( !$gene_stat ) { print "error, required file name is missing $0:  option --gene\n"; exit 1; }
	if( !$gene_list ) { print "error, required file name is missing $0:  option --list\n"; exit 1; }
	if( !$out_file )  { print "error, required file name is missing $0:  option --out\n"; exit 1; }
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
		'seq=s'    => \$seq_file,
		'gene=s'   => \$gene_stat,
		'list=s'   => \$gene_list,
		'out=s'    => \$out_file,
		'margin=i' => \$margin,
		'verbose'  => \$v,
		'debug'    => \$debug,
		'soft_mask'=> \$soft_mask,
		'max_gene=i'=> \$max_gene,
	);

	if( !$opt_results ) { print "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
	$v = 1 if $debug;
	
	# save informaton for debug
	$cfg->{'d'}->{'seq_file'}  = $seq_file;
	$cfg->{'d'}->{'gene_stat'} = $gene_stat;
	$cfg->{'d'}->{'gene_list'} = $gene_list;
	$cfg->{'d'}->{'out_file'}  = $out_file;
	$cfg->{'d'}->{'margin'}    = $margin;
	$cfg->{'d'}->{'v'}     = $v;
	$cfg->{'d'}->{'debug'} = $debug;
	$cfg->{'d'}->{'cmd'}   = $cmd;
	$cfg->{'d'}->{'soft_mask'} = $soft_mask;
	$cfg->{'d'}->{'max_gene'}  = $max_gene;
	
	print Dumper($cfg) if $debug;
};
# ------------------------------------------------
sub Usage
{
	print qq(# -------------------
Usage:  $0

Required options:
  --seq   [name] name of file with nucleotide contig sequence
                 file should be in FASTA format with unique ID in definishion line (first word in defline)
  --gene  [name] name of the file with gene coordinates on nucleotite contigs in YAML format
  --list  [name] name of the file with list of gene ID's
  --out   [name] ouput nucleotite sequence into this file in FASTA format

Optional parameters:
  --margin [number] margin around the genes, default: $margin
  --soft_mask       mask lowcase letters
  --max_gene [size] ignore genes longer than ; default is [0] use all; length on the DNA level (including introns)

Developer options:
  --verbose
  --debug
# -------------------
);
	exit 1;
}
# ================== END sub =====================
