#!/usr/bin/env perl
# ==============================================================
# Alex Lomsadze, Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
# 
# This script takes as input GFF file and creates a non redundant set
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
my $out_gff_file = '';
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

my %gff;

ParseGFF( $in_gff_file, \%gff );

PrintToFile( $out_gff_file, \%gff );

exit 0;

# ------------------------------------------------
sub PrintToFile
{
	my( $name, $ref ) = @_;
	
	my $line;

	open( my $OUT, ">$name" ) or die( "$!, error on open file $name" ); 
	foreach my $key (keys %$ref) 
	{

		if ($ref->{$key}{'al_score'} == -1) {
			$ref->{$key}{'att'} = ".";
		} else {
			$ref->{$key}{'att'} = "al_score=$ref->{$key}{'al_score'};";
		}

		addAttribute( $ref, $key, 'splice_sites' );
		addAttribute( $ref, $key, 'topProt' );

		$line = $ref->{$key}{'id'}     ."\t";
		$line .= $ref->{$key}{'info'}  ."\t";
		$line .= $ref->{$key}{'type'}  ."\t";
		$line .= $ref->{$key}{'start'} ."\t";
		$line .= $ref->{$key}{'end'}   ."\t";
		$line .= $ref->{$key}{'score'} ."\t";
		$line .= $ref->{$key}{'strand'}."\t";
		$line .= $ref->{$key}{'ph'}    ."\t";
		$line .= $ref->{$key}{'att'}   ."\n";

		print $OUT $line;
	} 
	close $OUT;
}
# ------------------------------------------------
sub addAttribute
{
	my( $ref, $key, $attr ) = @_;
	if ($ref->{$key}{$attr}) {
		if ($ref->{$key}{'att'} eq ".") {
			$ref->{$key}{'att'} = "$attr=$ref->{$key}{$attr};";
		} else {
			$ref->{$key}{'att'} .= " $attr=$ref->{$key}{$attr};";
		}
	}
}
# ------------------------------------------------
sub ParseGFF
{
	my( $name, $ref ) = @_;

	my $count_in_gff = 0;
	my $new_line;

	my $key;
	my $current;
	
	open( my $IN, $name ) || die "$! on open $name\n";
	while( my $line = <$IN> )
	{
		if( $line =~ /^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t([-+])\t(\S+)\t(.*)/ )
		{
			$key = $1 ."_". $3 ."_". $4 ."_". $5 ."_". $7;
			$ref->{$key}{'id'} = $1;
			$ref->{$key}{'info'} = $2;
			$ref->{$key}{'type'} = $3;
			$ref->{$key}{'start'} = $4;
			$ref->{$key}{'end'} = $5;
			$ref->{$key}{'strand'} = $7;
			$ref->{$key}{'ph'} = $8;

			$current = $6;
			my $attribute = $9;

			ParseAlignmentScore($attribute, $ref, $key);

			if ($attribute =~ /splice_sites=([^;]+);.*/) {
				$ref->{$key}{'splice_sites'} = $1;
			}

			if ($attribute =~ /topProt=([^;]+);.*/) {
				$ref->{$key}{'topProt'} = $1;
			}

			if ( $current eq '.' )
			{
				$ref->{$key}{'score'} += 1;
			}
			elsif ( $current =~ /^\d+$/ )
			{
				$ref->{$key}{'score'} += $current;
			}
			else
			{
				print "score $current\n";
			}
			
			++$count_in_gff;
		}
		else
		{
			print $line;
		}
	}
	close $IN;
	
	print "in: $count_in_gff\n" if $debug;	
}
# ------------------------------------------------
sub ParseAlignmentScore
{
	my( $att, $ref, $key) = @_;

	if (!defined $ref->{$key}{'al_score'})
	{
		$ref->{$key}{'al_score'} = -1;
	}

	if( $att =~ /al_score=([^;]+);.*/ )
	{
		if ($1 > $ref->{$key}{'al_score'})
		{
			$ref->{$key}{'al_score'} = $1;
		}
	}
}
# ------------------------------------------------
sub CheckBeforeRun
{
	print "check before run\n" if $debug;
		
	$bin      = ResolvePath( $bin );
	$work_dir = ResolvePath( $work_dir );
	
	$in_gff_file = ResolvePath( $in_gff_file );
	
	if( !$in_gff_file )  { print "error, required file name is missing $0:  option --in_gff\n"; exit 1; }
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
		'out_gff=s' => \$out_gff_file,
		'verbose' => \$v,
		'debug'   => \$debug
	);

	if( !$opt_results ) { print "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
	$v = 1 if $debug;

	# save informaton for debug
	$cfg->{'d'}->{'in_gff_file'}  = $in_gff_file;
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
  --in_gff   [name] input file generated by Spaln boundary scorer
  --out_gff  [name] output
 
Developer options:
  --verbose
  --debug
# -------------------
);
	exit 1;
}
# ================== END sub =====================
