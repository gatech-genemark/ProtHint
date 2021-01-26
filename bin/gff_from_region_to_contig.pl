#!/usr/bin/env perl
# ==============================================================
# Alex Lomsadze
# Copyright 2019, Georgia Institute of Technology, USA
# 
# This script takes as input GFF on region level and transfers coordinates to contig level
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
my $in_gff_file  = '';
my $seq_file     = '';
my $out_gff_file = '';
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

my %h;
ParseFastaDefline( $seq_file, \%h );

my @gff;
ParseGFF( $in_gff_file, \%h, \@gff );

PrintToFile( $out_gff_file, \@gff );

exit 0;

# ------------------------------------------------
sub PrintToFile
{
	my( $name, $ref ) = @_;
	
	open( my $OUT, ">$name" ) or die( "$!, error on open file $name" ); 
	foreach my $line ( @{$ref} )
	{
			print $OUT $line;
	} 
	close $OUT;;
}
# ------------------------------------------------
sub ParseGFF
{
	my( $name, $ref, $a ) = @_;

	my %d;
	my $count_in_gff = 0;
	my $new_line;
	my $start;
	my $end;
	
	open( my $IN, $name ) || die "$! on open $name\n";
	while( my $line = <$IN> )
	{
		if( $line =~ /^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t([-+])\t(\S+)\t(.*)/ )
		{
			$d{'id'} = $1;
			$d{'info'} = $2;
			$d{'type'} = $3;
			$d{'start'} = $4;
			$d{'end'} = $5;
			$d{'score'} = $6;
			$d{'strand'} = $7;
			$d{'ph'} = $8;
			$d{'att'} = $9;
			
			if( ! exists $ref->{ $d{'id'} } )
				{ print "error, ID of gene not found in regions file\n"; exit 1; }
				
			if( $ref->{$d{'id'}}{'strand'} eq '+' )
			{
				$start = $d{'start'} + $ref->{$d{'id'}}{'L'} - 1;
				$end   = $d{'end'}   + $ref->{$d{'id'}}{'L'} - 1;
				
				$new_line = $ref->{$d{'id'}}{'ID'} ."\t". $d{'info'} ."\t". $d{'type'} ."\t". $start ."\t". $end ."\t". $d{'score'} ."\t".'+'."\t". $d{'ph'} ."\t". $d{'att'} ." seed_gene_id=". $d{'id'} .";\n";
				
				print $new_line if $debug;
				
				push @{$a}, ($new_line); 
			}
			elsif( $ref->{$d{'id'}}{'strand'} eq '-' )
			{
				$end   = $ref->{$d{'id'}}{'R'} - $d{'start'} + 1;
				$start = $ref->{$d{'id'}}{'R'} - $d{'end'}   + 1;
				
				$new_line = $ref->{$d{'id'}}{'ID'} ."\t". $d{'info'} ."\t". $d{'type'} ."\t". $start ."\t". $end ."\t". $d{'score'} ."\t".'-'."\t". $d{'ph'} ."\t". $d{'att'} ." seed_gene_id=". $d{'id'} .";\n"; 
				
				print $new_line if $debug;
				
				push @{$a}, ($new_line);
			}
			else {die;}
			
			++$count_in_gff;
		}
	}
	close $IN;
	
	print "in: $count_in_gff\n" if $debug;	
}
# ------------------------------------------------
sub ParseFastaDefline
{
	my( $name, $ref ) = @_;

	open( my $IN, $name ) || die "$! on open $name\n";
	while( my $line = <$IN> )
	{
		if( $line =~ /^>(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+([-+])\s+(\S+)/ )
		{
			$ref->{$1}{"L"} = $2;
			$ref->{$1}{"R"} = $3;
			$ref->{$1}{"strand"} = $4;
			$ref->{$1}{"ID"} = $5;
		}
	}
 
	close $IN;	
}
# ------------------------------------------------
sub CheckBeforeRun
{
	print "check before run\n" if $debug;
		
	$bin      = ResolvePath( $bin );
	$work_dir = ResolvePath( $work_dir );
	
	$in_gff_file = ResolvePath( $in_gff_file );
	$seq_file    = ResolvePath( $seq_file );
	
	if( !$in_gff_file )  { print "error, required file name is missing $0:  option --in_gff\n"; exit 1; }
	if( !$seq_file )     { print "error, required file name is missing $0:  option --seq\n"; exit 1; }
	if( !$out_gff_file ) { print "error, required file name is missing $0:  option --out_gff\n"; exit 1; }
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
		'in_gff=s'  => \$in_gff_file,
		'seq=s'     => \$seq_file,
		'out_gff=s' => \$out_gff_file,
		'verbose' => \$v,
		'debug'   => \$debug
	);

	if( !$opt_results ) { print "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
	$v = 1 if $debug;

	# save informaton for debug
	$cfg->{'d'}->{'in_gff_file'}  = $in_gff_file;
	$cfg->{'d'}->{'seq_file'}     = $seq_file;
	$cfg->{'d'}->{'out_gff_file'} = $out_gff_file;
	$cfg->{'d'}->{'v'}     = $v;
	$cfg->{'d'}->{'debug'} = $debug;
	$cfg->{'d'}->{'cmd'}   = $cmd;
	
#	print Dumper($cfg) if $debug;
};
# ------------------------------------------------
sub Usage
{
	print qq(# -------------------
Usage:  $0

Required options:
  --in_gff   [name] name of file with coordinates on the region level
  --seq      [name] name of file with regions of sequence in FASTA format
  --out_gff  [name] output
 
Developer options:
  --verbose
  --debug
# -------------------
);
	exit 1;
}
# ================== END sub =====================
