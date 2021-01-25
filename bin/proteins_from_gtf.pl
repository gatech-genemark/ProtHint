#!/usr/bin/env perl
# ==============================================================
# Alex Lomsadze
# Copyright 2019, Georgia Institute of Technology, USA
#
# This script takes as input nucleotide sequence and gene coordinates
# and outputs protein sequences of genes
# ==============================================================

# --
# to do: finish option GFF
# --

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
my $seq_file   = '';
my $annot_file = '';
my $out_file   = '';

my $format = 'GTF';

my $min_protein_length  = 0; 
my $stat_annot_out_file = '';
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

my %genes;

if ( $format eq "GTF" )   { ReadGTF( $annot_file, \%genes ); }
#elsif ($format eq "GFF")  { ReadGFF( $annot_file, \%genes ); }
else { print "error, unexpected annotation format was specified $0: $format\n";  exit 1; }

my %sequence;

ReadSequence( $seq_file, \%sequence );
GetTranscripts( \%genes, \%sequence );

my %gcode;

TranslationTable_gcode_1( \%gcode );
GetProtein( \%genes, \%gcode );
PrintProt( $out_file, \%genes, $min_protein_length );

if( $stat_annot_out_file )
{
	PrintStat( $stat_annot_out_file );
}

exit 0;

#------------------------------------------------
sub PrintStat
{
	my $name = shift;
	
	my %seq_info;
	foreach my $key ( keys %sequence )
	{
		$seq_info{$key} = length( $sequence{$key} );
	}
	
	my %gene_info;
	foreach my $key ( keys %genes )
	{
		if( !exists $genes{$key}{'start'})
		{
			$gene_info{$key}{'start'} = -1;
		}
		else
		{
			$gene_info{$key}{'start'} = $genes{$key}{'start'};
		}
		if ( !exists $genes{$key}{'stop'} )
		{
			$gene_info{$key}{'stop'} = -1;
		}
		else
		{
			$gene_info{$key}{'stop'} = $genes{$key}{'stop'};
		}
		
		$gene_info{$key}{'R'} =  $genes{$key}{'R'};
		$gene_info{$key}{'L'} =  $genes{$key}{'L'};
		
		$gene_info{$key}{'strand'} = $genes{$key}{'strand'};
		$gene_info{$key}{'seq_id'} = $genes{$key}{'seq_id'};
	}
	
	open( my $OUT, ">$name") or die "error on open file $0: $name\n$!\n";
	print $OUT Dump( \%seq_info, \%gene_info );
	close $OUT;		
}
#------------------------------------------------
sub PrintProt
{
	my( $name, $g, $min ) = @_;

	my $count_all = 0;
	my $count_out = 0;

	open( my $OUT, ">$name" ) or die( "$!, error on open file $name" );
	foreach my $gid ( keys %{$g} )
	{
		if ( length( $g->{$gid}{'prot'} ) > $min )
		{
			print $OUT  ">$gid\n";
			print $OUT  $g->{$gid}{'prot'} ."\n\n";
			
			++$count_out;
		}
		
		++$count_all;
	}
	close $OUT;
	
	print "all prot: $count_all\n" if $debug;
	print "out prot: $count_out\n" if $debug;
};
#------------------------------------------------
sub dna_to_aa
{
    my( $s, $ref ) = @_;
    
    my $codon;
    my $AA = '';
    
    for( my $i=0; $i < length($s); $i += 3)
    {
        $codon = substr( $s, $i, 3 );

		if ( exists $ref->{$codon} )
		{
			$AA .= $ref->{$codon}
		}
		else
		{
			$AA .= 'X';
			print "warning, codon was not recognized $0: $codon\n" if $v;
		}
	}
	
	if( $AA =~ /\*$/ )
	{
		chop $AA;
	}
	
	return $AA;
};
#------------------------------------------------
sub GetProtein
{
	my( $g, $code ) = @_;

	foreach my $gid ( keys %{$g} )
	{
		$g->{$gid}{'prot'} = dna_to_aa( $g->{$gid}{'mrna'}, $code );
	}
};
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
# ------------------------------------------------
# in case of incomplete transcripts, trim beginning and end of sequence to full codon
# ------------------------------------------------
sub GetTranscripts
{
	my( $g, $s ) = @_;

	# just to simplify code reading
	my $seq_id;
	my $total_cds_count;
	my $strand;	
	my $L;
	my $R;
	my $ph;
	my $length_mod;

	foreach my $gid ( keys %{$g} )
	{
		$seq_id = $g->{$gid}{'seq_id'};
		$total_cds_count = scalar( @{ $g->{$gid}{"CDS"} } ) ;

		for( my $j = 0; $j < $total_cds_count; ++$j )
		{
			$L = $g->{$gid}{"CDS"}[$j]{'L'};
			$R = $g->{$gid}{"CDS"}[$j]{'R'};
			
			$g->{$gid}{'mrna'} .= substr( $s->{$seq_id}, $L - 1, $R - $L + 1 );
		}

		$strand = $g->{$gid}{'strand'};

		if( $strand eq '+' )
		{
			$ph = $g->{$gid}{"CDS"}[0]{'ph'};
			
			if( ! $ph )
			{
				;
			}
			elsif( $ph == 1 )
			{
				$g->{$gid}{'mrna'} =~ s/^.//;
			}
			elsif( $ph == 2 )
			{
				$g->{$gid}{'mrna'} =~ s/^..//;
			}
			else {die;}
			
			$length_mod = length($g->{$gid}{'mrna'}) % 3;
			
			if( ! $length_mod )
			{
				;
			}
			elsif( $length_mod == 1 )
			{
				chop $g->{$gid}{'mrna'};
			}
			elsif( $length_mod == 2 )
			{
				chop $g->{$gid}{'mrna'};
				chop $g->{$gid}{'mrna'};
			}
			else {die;}
		}
		elsif( $strand eq '-' )
		{
			$ph = $g->{$gid}{"CDS"}[$total_cds_count -1]{'ph'};
			
			if( ! $ph )
			{
				;
			}
			elsif( $ph == 1 )
			{
				chop $g->{$gid}{'mrna'};
			}
			elsif( $ph == 2 )
			{
				chop $g->{$gid}{'mrna'};
				chop $g->{$gid}{'mrna'};
			}
			else {die;}
			
			$length_mod = length($g->{$gid}{'mrna'}) % 3;
			
			if( ! $length_mod )
			{
				;
			}
			elsif( $length_mod == 1 )
			{
				$g->{$gid}{'mrna'} =~ s/^.//;
			}
			elsif( $length_mod == 2 )
			{
				$g->{$gid}{'mrna'} =~ s/^..//;
			}
			else {die;}
			
			$g->{$gid}{'mrna'} = RevComp( $g->{$gid}{'mrna'} );
		}
		else { print "error, strand information is missing $0: $gid\n"; exit 1; }
	}
}
# ------------------------------------------------
# this sub is modification of code from John Besemer (2000 GaTech)
# ------------------------------------------------
sub TranslationTable_gcode_1
{
	my $ref = shift;
    $ref->{"TTT"}="F";  $ref->{"TTC"}="F"; $ref->{"TTA"}="L"; $ref->{"TTG"}="L";
    $ref->{"CTT"}="L";  $ref->{"CTC"}="L"; $ref->{"CTA"}="L"; $ref->{"CTG"}="L";
    $ref->{"ATT"}="I";  $ref->{"ATC"}="I"; $ref->{"ATA"}="I"; $ref->{"ATG"}="M";
    $ref->{"GTT"}="V";  $ref->{"GTC"}="V"; $ref->{"GTA"}="V"; $ref->{"GTG"}="V";
    $ref->{"TCT"}="S";  $ref->{"TCC"}="S"; $ref->{"TCA"}="S"; $ref->{"TCG"}="S";
    $ref->{"CCT"}="P";  $ref->{"CCC"}="P"; $ref->{"CCA"}="P"; $ref->{"CCG"}="P";
    $ref->{"ACT"}="T";  $ref->{"ACC"}="T"; $ref->{"ACA"}="T"; $ref->{"ACG"}="T";
    $ref->{"GCT"}="A";  $ref->{"GCC"}="A"; $ref->{"GCA"}="A"; $ref->{"GCG"}="A";
    $ref->{"TAT"}="Y";  $ref->{"TAC"}="Y"; $ref->{"TAA"}="*"; $ref->{"TAG"}="*";
    $ref->{"CAT"}="H";  $ref->{"CAC"}="H"; $ref->{"CAA"}="Q"; $ref->{"CAG"}="Q";
    $ref->{"AAT"}="N";  $ref->{"AAC"}="N"; $ref->{"AAA"}="K"; $ref->{"AAG"}="K";
    $ref->{"GAT"}="D";  $ref->{"GAC"}="D"; $ref->{"GAA"}="E"; $ref->{"GAG"}="E";
    $ref->{"TGT"}="C";  $ref->{"TGC"}="C"; $ref->{"TGA"}="*"; $ref->{"TGG"}="W";
    $ref->{"CGT"}="R";  $ref->{"CGC"}="R"; $ref->{"CGA"}="R"; $ref->{"CGG"}="R";
    $ref->{"AGT"}="S";  $ref->{"AGC"}="S"; $ref->{"AGA"}="R"; $ref->{"AGG"}="R";
    $ref->{"GGT"}="G";  $ref->{"GGC"}="G"; $ref->{"GGA"}="G"; $ref->{"GGG"}="G";
};
# ------------------------------------------------
sub ReadSequence
{
	my ($name, $ref ) = @_;
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
			if( $id eq '' ) { print "error, fasta record without definition line found $0: $line\n"; exit 1; }

			$line = uc $line;

			# remove non alphabet
			$line =~ tr/A-Z//dc;

      		# replace allowed nucleic acid code (non A T C G) by N
			$line =~ tr/RYKMSWBDHV/N/;

			# stop if unexpected character
			if( $line =~ m/[^ATCGN]/ )
				{ print "error, unexpected letter found in sequence $0: $line\n"; exit 1; }

      		$ref->{ $id } .= $line;
		}
	}
 
	close $IN;
	
	if( $debug )
	{
		print "id length\n";
		foreach my $id (keys %{$ref} )
		{
			print ( $id ." ". length($ref->{$id}) ."\n" );
		}
	}
};
# ------------------------------------------------
sub ReadGTF
{
	my( $name, $ref ) = @_;
	open( my $IN, $name ) or die "$! on open file $name\n";

	my $seq_id;
	my $type;
	my $left;
	my $right;
	my $strand;
	my $phase;
	my $gene_id;
	
	my $count_cds = 0;
	
	while( my $line = <$IN> )
	{

# file format
#                       1       2      3     4      5     6       7        8                        9
#                    seq_id   source type  left   right score   strand   phase                    gene_id

		if( $line =~ /^(.+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t([\+-])\t([012.])\t.*gene_id\s+\"\s*(\S+)\s*\"\s*;/  )
		{
			$seq_id  = $1;			
			$type    = $3;
			$left    = $4;
			$right   = $5;
			$strand  = $7;
			$phase   = $8;
			$gene_id = $9;

			if ( $seq_id =~ /^(\S+)\s*/ )
			{
				$seq_id = $1;
			} 

			if ( exists $ref->{$gene_id} )
			{
				if( $ref->{ $gene_id }{ "strand" } ne $strand )
					{ print "error, mismatch in strand $0 : $line\n"; exit 1; }
					
				if( $ref->{ $gene_id }{"seq_id"} ne $seq_id )
					{ print "error, mismatch in sequence ID $0 : $line\n"; exit 1; }		
			}
			else
			{
				$ref->{ $gene_id }{"strand"} = $strand;
				$ref->{ $gene_id }{"seq_id"} = $seq_id;
			}
	
			if( $type eq "CDS" )
			{	
				push @{ $ref->{ $gene_id }{"CDS"} }, { ph=>$phase, L=>$left, R=>$right };
				
				++$count_cds;
				
				if ( !exists $ref->{ $gene_id }{'L'} )
				{
					$ref->{ $gene_id }{'L'} = $left;
					$ref->{ $gene_id }{'R'} = $right;
				}
				else
				{
					if ( $ref->{ $gene_id }{'L'} > $left )
					{
						$ref->{ $gene_id }{'L'} = $left;
					}

					if ( $ref->{ $gene_id }{'R'} < $right )
					{
						$ref->{ $gene_id }{'R'} = $right;
					}
				}
			}
			elsif( $type eq "start_codon" )
			{
				if ( exists $ref->{$gene_id} and exists $ref->{ $gene_id }{"start"} )
					{ print "error, two start codons found $0 : $line\n"; exit 1; }
				else
				{
					if ( $strand eq '+' )
					{
						$ref->{ $gene_id }{"start"} = $left;
					}
					else
					{
						$ref->{ $gene_id }{"start"} = $right;
					}
				}
			}
			elsif( $type eq "stop_codon" )
			{
				if ( exists $ref->{$gene_id} and exists $ref->{ $gene_id }{"stop"} )
					{ print "error, two stop codons found $0 : $line\n"; exit 1; }
				else
				{
					if ( $strand eq '+' )
					{
						$ref->{ $gene_id }{"stop"} = $right;
					}
					else
					{
						$ref->{ $gene_id }{"stop"} = $left;
					}
				}
			}
			
		}
		else
		{
			if ( $line =~ /^\s+$/ ) { next };
			if ( $line =~ /^#/ )    { next };
			
			print "error, unexpected format found in file $0 $name: $line\n";
			exit 1;
		}
	}

	close $IN;
	
	print "CDS lines found: $count_cds\n" if $debug;
};
# ------------------------------------------------
sub CheckBeforeRun
{
	print "check before run\n" if $debug;
	
	if ( $min_protein_length < 0 )
		{ print "error, out of range values specified for minimum protein length $0: $min_protein_length\n"; exit 1; }
		
	$bin      = ResolvePath( $bin );
	$work_dir = ResolvePath( $work_dir );
	
	$seq_file   = ResolvePath( $seq_file );
	$annot_file = ResolvePath( $annot_file );
	
	if( !$seq_file )   { print "error, required file name is missing $0:  option --seq\n"; exit 1; }
	if( !$annot_file ) { print "error, required file name is missing $0:  option --annot\n"; exit 1; }
	if( !$out_file )   { print "error, required file name is missing $0:  option --out\n"; exit 1; }
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
		'annot=s'  => \$annot_file,
		'out=s'    => \$out_file,
		'format=s' => \$format,
		'min=i'    => \$min_protein_length,
		'stat=s'   => \$stat_annot_out_file,
		'verbose'  => \$v,
		'debug'    => \$debug
	);

	if( !$opt_results ) { print "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
	$v = 1 if $debug;
	
	# save information for debug
	$cfg->{'d'}->{'seq_file'}   = $seq_file;
	$cfg->{'d'}->{'annot_file'} = $annot_file;
	$cfg->{'d'}->{'out_file'}   = $out_file;
	$cfg->{'d'}->{'format'}     = $format;
	$cfg->{'d'}->{'min_protein_length'}  = $min_protein_length;
	$cfg->{'d'}->{'stat_annot_out_file'} = $stat_annot_out_file;
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
  --seq   [name] name of file with nucleotide sequence/s
                 file should be in FASTA format with unique ID in definition line (first word in defline)
  --annot [name] name of file with gene coordinates on nucleotide contigs
                 file should be in GFF or GTF format; each gene should have unique gene ID without white spaces in it
  --out   [name] output protein sequence into this file in FASTA format
                 definition line consist of gene unique ID

Optional parameters:
  --min  [number]  minimum length of the protein to output in AA
  --stat [name]    output statistics of gene/intron/intergenic length into this file
  --format [GTF or GFF] format of file with gene annotation

Developer options:
  --verbose
  --debug
# -------------------
);
	exit 1;
}
# ================== END sub =====================
