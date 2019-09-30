#!/usr/bin/env perl
# ==============================================================
# Alex Lomsadze, Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# This script takes as input alignment info from ProSplign in ASN
# and selects subset of info in GFF format
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
my $asn_file = '';
my $out_file = '';
my $exons = 0;
my $augustus = 0;
my $exonCutoff = 0;
my $closeThreshold = 0;
my $printStats = 0;
my $MAX_AFTER_STOP = 5;
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

my %h;
my %p;
my @arr;

my $txt;
my $last_read = '';
my $endoffile = 0;

open( my $ASN, $asn_file ) or die( "$!, error on open file $asn_file" );
open( my $OUT, ">$out_file" ) or die( "$!, error on open file $out_file" );

my $startLabel="start_codon";
my $stopLabel="stop_codon";
if ($augustus) {
	$startLabel = "start";
	$stopLabel = "stop";
}

while ( !$endoffile and  $txt = ReadFromFile( $ASN )  )
{	
	if( $txt =~ /^\s*Seq-align\s*::=\s*/ ) {;}
	else{die"unexpected format found: $txt\n";}
	
	%h = ();
	%p = ();
	@arr = ();
	
	$h{'score'}{'num_ident'} = 0;
	$h{'score'}{'num_positives'} = 0;
	$h{'score'}{'num_negatives'} = 0;
	$h{'score'}{'product_gap_length'} = 0;
	$h{'score'}{'genomic_gap_length'} = 0;
	$h{'score'}{'align_length'} = 0;
	
	$h{'seq'}{'product-id'} = '';
	$h{'seq'}{'genomic-id'} = '';
	$h{'seq'}{'genomic-strand'} = '';
	$h{'seq'}{'product-length'} = 0;
	
	$h{'exons'} = ();
	$h{'introns'} = ();

	MatchScores();

	my $startFound = 0;
	my $stopFound = 0;
	MatchSeq(\$startFound, \$stopFound);

	FindPoints();

	MatchExons();
	ParseIntrons();

	my $close = "FALSE";

	my $align_length = ($h{'score'}{'align_length'});
	if ($align_length == 0) {
		$align_length = 1
	}

	my $pID = ($h{'score'}{'num_ident'} / $align_length) * 100;
	
	if ($pID >= $closeThreshold) {
		AddTo($OUT, $pID, $startFound, $stopFound);
	}
	

	print Dumper(\%h) if $debug;
}

close $OUT;
close $ASN;

exit 0;

# ------------------------------------------------
sub AddTo
{
	my( $F, $pID, $startFound, $stopFound) = @_;

	my $gff_line;
	$pID = sprintf("%.2f", $pID );

	my $fullProteinAligned = checkFullProteinAligned();

	if ( defined $h{'introns'} )
	{	
		foreach my $i ( 0 .. $#{ $h{'introns'} } )
		{
			$gff_line  = $h{'seq'}{'genomic-id'} ."\tProSplign\tintron\t";
			$gff_line .= $h{'introns'}[$i]{'left'} ."\t";
			$gff_line .= $h{'introns'}[$i]{'right'} ."\t";
			$gff_line .= ".\t+\t.\t";
			if (!$augustus) {
				$gff_line .= ("prot=". $h{'seq'}{'product-id'} .
					"; ssite=". $h{'introns'}[$i]{'ssites'} .
					"; fullProteinAligned=" . $fullProteinAligned .
					"; pID=" . $pID  .
					"; intron_id=". ($i + 1) . ";\n");
				print $F  $gff_line;
			} else {
				$gff_line .= ("grp=" . $h{'seq'}{'genomic-id'} . "_" . $h{'seq'}{'product-id'} .
				";src=P;pri=4;\n");
				print $F  $gff_line;
			}
		} 
	}

	if (defined $h{'exons'} )
	{
		if ($startFound && $h{'exons'}[0]{'prot-start'} == 1) {
			$gff_line  =  $h{'seq'}{'genomic-id'} ."\tProSplign\t$startLabel\t";
			$gff_line .=  $h{'exons'}[0]{'genomic-start'} ."\t";
			$gff_line .=  $h{'exons'}[0]{'genomic-start'} + 2 . "\t";
			$gff_line .=  ".\t+\t0\t";
			if (!$augustus) {
				$gff_line .=  ( "prot=" . $h{'seq'}{'product-id'} .
					"; fullProteinAligned=" . $fullProteinAligned .
					"; pID=" . $pID . ";\n");
				print $F  $gff_line;
			} else {
				$gff_line .= ("grp=" . $h{'seq'}{'genomic-id'} . "_" . $h{'seq'}{'product-id'} .
				";src=P;pri=4;\n");
				print $F  $gff_line;
			}
		}

		my $last = $#{$h{'exons'}};

		if ($stopFound && $h{'exons'}[$last]{'prot-end'} >= $h{'seq'}{'product-length'} - $MAX_AFTER_STOP) {
			$gff_line  =  $h{'seq'}{'genomic-id'} ."\tProSplign\t$stopLabel\t";
			$gff_line .=  $h{'exons'}[$last]{'genomic-end'} + 1 ."\t";
			$gff_line .=  $h{'exons'}[$last]{'genomic-end'} + 3 ."\t";
			$gff_line .=  ".\t+\t0\t";
			if (!$augustus) {
				$gff_line .= ( "prot=" . $h{'seq'}{'product-id'} .
					"; fullProteinAligned=" . $fullProteinAligned .
					"; pID=" . $pID . ";\n");
				print $F  $gff_line;
			} else {
				$gff_line .= ("grp=" . $h{'seq'}{'genomic-id'} . "_" . $h{'seq'}{'product-id'} .
				";src=P;pri=4;\n");
				print $F  $gff_line;
			}
		}

		if ($exons) {
			foreach my $i ( 0 .. $#{ $h{'exons'} } )
			{

				cutExon($h{'exons'}[$i]);

				$gff_line  =  $h{'seq'}{'genomic-id'} ."\tProSplign\t";
				$gff_line  .=  $h{'exons'}[$i]{'label'} ."\t";
				$gff_line .=  $h{'exons'}[$i]{'genomic-start'} ."\t";
				$gff_line .=  $h{'exons'}[$i]{'genomic-end'} ."\t";
				$gff_line .=  ".\t+\t" . $h{'exons'}[$i]{'frame'} . "\t";
				if (!$augustus) {
					$gff_line .=  ("prot=" . $h{'seq'}{'product-id'} .
						"; prot_start=" . $h{'exons'}[$i]{'prot-start'} .
						"; prot_end=" . $h{'exons'}[$i]{'prot-end'} .
						"; partial=" . $h{'exons'}[$i]{'partial'} .
						"; fullProteinAligned=" . $fullProteinAligned .
						"; pID=" . $pID .
						"; frameshift=" . $h{'exons'}[$i]{'frameshifted'} .
						"; exon_id=" . ($i + 1) .";\n");
					print $F  $gff_line;
				} else {
					if ($h{'exons'}[$i]{'frameshifted'} eq 'FALSE') {
						$gff_line .= ("grp=" . $h{'seq'}{'genomic-id'} . "_" . $h{'seq'}{'product-id'} .
						";src=P;pri=4;\n");
						print $F  $gff_line;
					}
				}

			}
		}

		print $F "\n";
	}

}
# ------------------------------------------------
sub checkFullProteinAligned
{
	if (defined $h{'exons'} )
	{

		if ($h{'exons'}[0]{'prot-start'} != 1) {
			return "FALSE";
		}

		if ($h{'exons'}[$#{$h{'exons'}}]{'prot-end'} != $h{'seq'}{'product-length'}) {
			return "FALSE";
		}

		foreach my $i (1 .. $#{$h{'exons'}}){
			if  ($h{'exons'}[$i]{'prot-start'} - $h{'exons'}[$i - 1]{'prot-end'} > 1) {
				return "FALSE";
			}
		}
	}
	return "TRUE";
}
# ------------------------------------------------
sub SaveTo
{
	my( $name ) = @_;
	
	my $gff_line;
	
	open( my $OUT, ">$name" ) or die( "$!, error on open file $name" ); 

	if ( defined $h{'introns'} )
	{	
		foreach my $i ( 0 .. $#{ $h{'introns'} } )
		{
			$gff_line  = $h{'seq'}{'genomic-id'} ."\tProSplign\tIntron\t";
			$gff_line .=  $h{'introns'}[$i]{'left'} ."\t";
			$gff_line .=  $h{'introns'}[$i]{'right'} ."\t";
			$gff_line .= ( ".\t+\t.\tprot=". $h{'seq'}{'product-id'}. "; intron_id=". ($i + 1) .";\n" );
			
			print $OUT $gff_line;
		} 
	}
	
	close $OUT;;
}
# ------------------------------------------------
sub ParseIntrons
{
	if( !defined $h{"exons"} ) {return;}
	
	my $arr_size = scalar @{$h{"exons"}};
	
	if( $arr_size > 1 )
	{
		my $prev = 0;
		my $next = 1;
		my $ref = $h{'exons'};
		
		while( $next < $arr_size )
		{
			if( ($ref->[$prev]{'donor-after-exon'} eq 'GT' ) and ($ref->[$next]{'acceptor-before-exon'} eq 'AG' ) )
			{
				my %intron;
				$intron{'left'}  = $ref->[$prev]{'genomic-end'} + 1;
				$intron{'right'} = $ref->[$next]{'genomic-start'} - 1;
				$intron{'ssites'} = "GT_AG";
				
				push @{$h{'introns'}}, {%intron};
			}
			elsif ( ($ref->[$prev]{'donor-after-exon'} eq 'GC' ) and ($ref->[$next]{'acceptor-before-exon'} eq 'AG' ) )
			{
				my %intron;
				$intron{'left'}  = $ref->[$prev]{'genomic-end'} + 1;
				$intron{'right'} = $ref->[$next]{'genomic-start'} - 1;
				$intron{'ssites'} = "GC_AG";

				push @{$h{'introns'}}, {%intron};
			}
			elsif ( ($ref->[$prev]{'donor-after-exon'} eq 'AT' ) and ($ref->[$next]{'acceptor-before-exon'} eq 'AC' ) )
			{
				my %intron;
				$intron{'left'}  = $ref->[$prev]{'genomic-end'} + 1;
				$intron{'right'} = $ref->[$next]{'genomic-start'} - 1;
				$intron{'ssites'} = "AT_AC";

				push @{$h{'introns'}}, {%intron};
			}
			
			++$prev;
			++$next;
		}
	}
}
# ------------------------------------------------
sub FindEndOf
{
	my $i = shift;
	my $j = 0;
	
	my $current_level = $p{$i}{'level'};
	
	for( my $at = $p{$i}{'at'} + 1; $at < scalar(@arr); ++$at )
	{
		if( $p{$arr[$at]}{'level'} <= $current_level )
		{
			$j = $arr[$at];
			last;
		}
	}
	
	return $j;
}
# ------------------------------------------------
sub cutExon
{
	my $exon = shift;
	my $length = $exon->{'genomic-end'} - $exon->{'genomic-start'} + 1;

	# If the exon is too short after cutoff, report just a codon in the middle of the exon
	# to preserve frame information
	if ($length < 2 * $exonCutoff + 3) {
		$exon->{'genomic-start'} += int($length / 6) * 3;
		$exon->{'genomic-end'} = $exon->{'genomic-start'} + 2;
	} else {
		$exon->{'genomic-start'} += $exonCutoff;
		$exon->{'genomic-end'} -= $exonCutoff;
	}
}
# ------------------------------------------------
sub determineFrame
{
	my ($str, $exon) = @_;

	my $prosplignFrame;
	if ($$str =~ /product-start protpos [^,]+,\sframe (\d+)/) {
		$prosplignFrame = $1;
	} else {
		die;
	}

	if ($prosplignFrame == 1) {
		$exon->{'frame'} = 0;
	} elsif ($prosplignFrame == 2) {
		$exon->{'frame'} = 2;
	} else {
		$exon->{'frame'} = 1;
	}
}
# ------------------------------------------------
sub isPartial
{
	my ($str, $exon) = @_;

	if( $$str =~ /partial (\S+)/ ) {
		$exon->{'partial'} = $1;
	}  else {
		$exon->{'partial'} = 'FALSE';
	}

	$exon->{'label'} = 'CDSpart';

	if ($exon->{'partial'} eq 'FALSE' && $exonCutoff == 0) {
		$exon->{'label'} = 'CDS';
	}
}
sub checkFrameshifts
{
	my ($str, $exon) = @_;

	$exon->{'frameshifted'} = 'FALSE';

	# Gaps in the middle of exon alignment non divisible by 3 cause frameshifts.
	# If the non-divisible gap is at the start or end of an exon, it just means
	# the gap is against a split codon, not frameshift.

	while($$str =~ /, genomic-ins (\d+), diag/g) {
		if ($1 % 3 != 0) {
			$exon->{'frameshifted'} = 'TRUE';
		}
	}

	while($$str =~ /, product-ins (\d+), diag/g) {
		if ($1 % 3 != 0) {
			$exon->{'frameshifted'} = 'TRUE';
		}
	}
}
# ------------------------------------------------
sub MatchOneExon
{
	my $str = shift;
	
	my %exon;

	if( $str =~ /genomic-start (\d+)/ ) { $exon{'genomic-start'} = $1 + 1; }     else{die;}
	if( $str =~ /genomic-end (\d+)/ )   { $exon{'genomic-end'} = $1 + 1; }       else{die;}

	if ($str =~ /product-start protpos \{\samin (\d+)/) {
		$exon{'prot-start'} = $1 + 1;
	} else {
		die;
	}

	if ($str =~ /product-end protpos \{\samin (\d+)/) {
		$exon{'prot-end'} = $1 + 1;
	} else {
		die;
	}

	determineFrame(\$str, \%exon);
	isPartial(\$str, \%exon);
	checkFrameshifts(\$str, \%exon);

	if( $str =~ /donor-after-exon \{ bases \"(\S+)\"/ ) { $exon{'donor-after-exon'} = $1; }  else { $exon{'donor-after-exon'} = ''; }
	if( $str =~ /acceptor-before-exon \{ bases \"(\S+)\"/ ) { $exon{'acceptor-before-exon'} = $1; }  else { $exon{'acceptor-before-exon'} = ''; }
	
	push @{$h{'exons'}}, {%exon};
}
# ------------------------------------------------
sub MatchExons
{
	if( $txt =~ /exons \{/ )
	{
		my $i = $+[0] - 1;
		my $j = FindEndOf($i);
		
		print $i ." ". $j ."\n" if $debug;
		
		# exon section is between $i and $j
		
		my $next = $p{$i}{'at'} + 1;
		
		while( $arr[$next] < $j )
		{
			my $start = $arr[ $next ];
			my $end   = FindEndOf($start);
			
			print "$next $start $end\n" if $debug;

			MatchOneExon( substr( $txt, $start, $end - $start ));
			
			$next = $p{$end}{'at'} + 1;
		}
	}   
}
# ------------------------------------------------
sub FindPoints
{
	my $level = 0;
	my $char;
	my $i;
	
	while( $txt =~ m/[\{\}]/g )
	{
		$i = $-[0];
		$char = substr($txt, $i, 1 );
		
		push @arr, ( $i ); # the same as sorted keys of hash
		$p{$i}{'at'} = scalar( @arr ) - 1;

		$p{$i}{'char'} = $char;
		
		if( $char eq '{' )
		{
			++$level;
			$p{$i}{'level'} = $level;
		}
		elsif( $char eq '}' )
		{
			$p{$i}{'level'} = $level;
			--$level;
		}
#		elsif( $char eq ',' )
#		{
#			$p{$i}{'level'} = $level;
#		}
		else {die;}
		
#		print ( $i ." ". $p{$i}{'char'} ." ". $p{$i}{'level'} ." ". $p{$i}{'at'} ."\n") if $debug;
	}
}
# ------------------------------------------------
sub MatchSeq
{	
	my( $startFound, $stopFound ) = @_;

	if( $txt =~ /product-id local str \"(\S+)\"/ ) { $h{'seq'}{'product-id'} = $1; }     else{die;}
	if(($txt =~ /genomic-id local str \"(\S+)\"/ ) or ( $txt =~ /genomic-id local id (\d+)/ ))
						       { $h{'seq'}{'genomic-id'} = $1; }     else{die "$txt\n";}
	if( $txt =~ /genomic-strand (plus)/ )          { $h{'seq'}{'genomic-strand'} = $1; } else{die;}
	if( $txt =~ /product-length (\d+)/ )           { $h{'seq'}{'product-length'} = $1; } else{die;}


	if( $txt =~ /start-codon-found TRUE/ ) {
		$$startFound = 1;
	}

	if( $txt =~ /stop-codon-found TRUE/ ) {
		$$stopFound = 1;
	}
}
# ------------------------------------------------
sub MatchScores
{	
	if( $txt =~ /num_ident\", value int (\d+)/ )          { $h{'score'}{'num_ident'} = $1; }          else{die;}
	if( $txt =~ /num_positives\", value int (\d+)/ )      { $h{'score'}{'num_positives'} = $1; }      else{die;}
	if( $txt =~ /num_negatives\", value int (\d+)/ )      { $h{'score'}{'num_negatives'} = $1; }      else{die;}
	if( $txt =~ /product_gap_length\", value int (\d+)/ ) { $h{'score'}{'product_gap_length'} = $1; } else{die;}
	if( $txt =~ /genomic_gap_length\", value int (\d+)/ ) { $h{'score'}{'genomic_gap_length'} = $1; } else{die;}
	if( $txt =~ /align_length\", value int (\d+)/ )       { $h{'score'}{'align_length'} = $1; }       else{die;}
}
# ------------------------------------------------
sub ReadFromFile
{
	my( $F ) = @_;
	
	my $str = $last_read;
	my $line_counter = 0;
	
	my $line = '';
	
	while( $line = <$F> )
	{
		++$line_counter;

		if( $line =~ /^\s*$/ ) {next;}
		if( $line =~ /^\s*#/ ) {next;}

		chomp $line;

		if( $line =~ /^\s*Seq-align\s*::=\s*/ )
		{	
			# ini
			if ( !$last_read )
			{
				$last_read = $line;
			}
			else
			{
				$last_read = $line;
				last;
			}
		}

		$str .= ($line ." ");
	}
	
	if ( !defined $line ) { $endoffile = 1; }
	
	print "lines in record: $line_counter\n" if $debug;
	
	$str =~ s/\s+/ /g;
	
	return $str;
}
# ------------------------------------------------
sub ReadFile
{
	my( $name ) = @_;
	
	my $str = '';
	my $line_counter = 0;
	
	open( my $IN, $name ) || die "$! on open $name\n";
	while( my $line = <$IN> )
	{
		++$line_counter;
		
		if( $line =~ /^\s*$/ ) {next;}
		if( $line =~ /^\s*#/ ) {next;}
		
		chomp $line;
		$str .= ($line ." ");
	}
	close $IN;
	
	print "lines in file: $line_counter\n" if $debug;
	
	$str =~ s/\s+/ /g;
	
	return $str;
}
# ------------------------------------------------
sub CheckBeforeRun
{
	print "check before run\n" if $debug;
		
	$bin      = ResolvePath( $bin );
	$work_dir = ResolvePath( $work_dir );
	
	$asn_file  = ResolvePath( $asn_file );

	if ($exonCutoff < 0 || $exonCutoff % 3 != 0) {
		print "error, exonCutoff must be positive and divisible by 3.\n";
		exit 1;
	}

	if( !$asn_file ) { print "error, required file name is missing $0:  option --asn\n"; exit 1; }
	if( !$out_file ) { print "error, required file name is missing $0:  option --out\n"; exit 1; }
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
		'asn=s'        => \$asn_file,
		'out=s'        => \$out_file,
		'exons'        => \$exons,
		'augustus'     => \$augustus,
		'exonCutoff=i'  => \$exonCutoff,
		'close=f'      => \$closeThreshold,
		'stats'        => \$printStats,
		'verbose'      => \$v,
		'debug'        => \$debug
	);

	if( !$opt_results ) { print "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
	$v = 1 if $debug;

	# save informaton for debug
	$cfg->{'d'}->{'asn_file'}  = $asn_file;
	$cfg->{'d'}->{'out_file'}  = $out_file;
	$cfg->{'d'}->{'exons'}  = $exons;
	$cfg->{'d'}->{'augustus'}  = $augustus;
	$cfg->{'d'}->{'exonCutoff'}  = $exonCutoff;
	$cfg->{'d'}->{'close'}  = $closeThreshold;
	$cfg->{'d'}->{'stats'}  = $printStats;
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
  --asn   [name]  name of file with ProSplign output
  --out   [name]  spliced alignment by ProSlign

Optional options:
  --exons             print exons
  --augustus 	      print exons in augustus format
  --close [float]     Print only data from alignments with percent identity higher than this threshold. Default = 0.
  --stats             Print statistics about the alignments
  --exonCutoff [int]  Cut this many bases at the start and end from each exon. The cutoff must be divisible
                      by 3 to preserve the same reading frame. If the exon is shorter than 2 * exonCutoff, a
                      single codon in the middle of the exon is reported to preserve the reading frame information.

Developer options:
  --verbose
  --debug
# -------------------
);
	exit 1;
}
# ================== END sub =====================
