#!/usr/bin/env python
# Author: Tomas Bruna

# TODO:
# * Chained hints for Augustus


import argparse
import os
import sys
import subprocess
import multiprocessing


workDir = ''
binDir = ''
genome = ''
proteins = ''
threads = ''


def main():
    args = parseCmd()
    setEnvironment(args)

    geneMarkGtf = args.geneMarkGtf
    if not geneMarkGtf:
        geneMarkGtf = runGeneMarkES(args.pbs)

    translateSeeds(geneMarkGtf)

    diamondPairs = args.diamondPairs
    if not diamondPairs:
        diamondPairs = runDiamond(args.diamondBin, args.maxProteinsPerSeed, args.evalue)

    prepareSeedSequences(diamondPairs)
    runSpaln(diamondPairs, args.pbs)
    filterSpalnPairs(args.maxSpalnCoverage)

    prepareProSplignPairs(diamondPairs, args.ensureDiamondPairs)
    runProSplign(args.pbs)

    processOutput(args.coverageThreshold, args.alignmentScoreThreshold, args.closeThreshold)

    if cleanup:
        cleanup()


def runGeneMarkES(pbs):
    """Run GeneMark-ES

    Args:
        pbs (boolean): Whether to run on pbs

    Returns:
        string: Path to genemark gff output
    """
    ESDir = workDir + "/GeneMark_ES"
    if not os.path.isdir(ESDir):
        os.mkdir(ESDir)
    os.chdir(ESDir)

    command = ""
    if not pbs:
        command = binDir + "/../dependencies/GeneMarkES/bin/gmes_petap.pl --verbose --cores " + threads + \
                  " --max_intergenic 50000 --ES --seq " + genome + " --soft 1000"
    else:
        command = binDir + "/../dependencies/GeneMarkES/bin/gmes_petap.pl --verbose --pbs \
                  --max_intergenic 50000 --ES --seq " + genome + " --soft 1000"

    subprocess.call(command, shell=True)

    return os.path.abspath("genemark.gtf")


def translateSeeds(geneMarkGtf):
    """Translate GeneMark-ES seeds to proteins

    Args:
        geneMarkGtf (filepath): Path to GeneMark.gtf prediction file
    """
    os.chdir(workDir)

    command = binDir + "/proteins_from_gtf.pl --stat gene_stat.yaml --seq " + \
              genome + " --annot " + geneMarkGtf + " --out gmes_proteins.faa --format GTF"
    subprocess.call(command, shell=True)


def runDiamond(diamondBin, maxProteins, evalue):
    """Run DIAMOND protein search

    Args:
        diamondBin (executable): Path to diamond binary. If none specified,
                                 try to find it in the protein mapping bin folder
        maxProteins (int): Maximum number of protein hits per seed gene.
        evalue (float): Maximum e-value of DIAMOND hits

    Returns:
        string: Path to DIAMOND output
    """
    if not diamondBin:
        diamondBin = binDir + "/../dependencies/diamond"

    if not os.path.isfile(diamondBin):
        sys.stderr.write("error: DIAMOND binary was not found\n")
        sys.exit()

    diamondDir = workDir + "/diamond"
    if not os.path.isdir(diamondDir):
        os.mkdir(diamondDir)
    os.chdir(diamondDir)

    # Make DIAMOND db
    command = diamondBin + " makedb --in " + proteins + " -d diamond_db"
    subprocess.call(command, shell=True)

    # Actual DIAMOND run
    command = diamondBin + " blastp --query ../gmes_proteins.faa --db diamond_db \
              --outfmt 6 qseqid sseqid --out diamond.out --max-target-seqs " + \
              str(maxProteins) + " --max-hsps 1 --threads " + threads + \
              " --evalue " + str(evalue)
    subprocess.call(command, shell=True)

    return os.path.abspath("diamond.out")


def prepareSeedSequences(diamondPairs):
    """Prepare nucleotide sequences for seed genes

    Args:
        diamondPairs (filepath): Path to file with seed gene-protein pairs
    """
    os.chdir(workDir)

    command = binDir + "/nucseq_for_selected_genes.pl --seq " + genome + \
              " --out nuc.fasta --gene gene_stat.yaml --list " + diamondPairs
    subprocess.call(command, shell=True)


def runSpaln(diamondPairs, pbs):
    """Run Spaln spliced alignment

    Args:
        diamondPairs (filePath): Path to file with seed gene-protein pairs to align
        pbs (boolean): Whether to run on pbs
    """
    spalnDir = workDir + "/Spaln"
    if not os.path.isdir(spalnDir):
        os.mkdir(spalnDir)
    os.chdir(spalnDir)

    command = ""
    if not pbs:
        command = binDir + "/run_spliced_alignment.pl --cores " + threads + \
                  " --nuc ../nuc.fasta --list " + diamondPairs + \
                  " --prot " + proteins + " --v --aligner spaln  > loginfo_spaln"
    else:
        command = binDir + "/run_spliced_alignment_pbs.pl --N 120 --K 8 --seq \
                  ../nuc.fasta --list " + diamondPairs + " --db " + \
                  proteins + " --v --aligner spaln  > loginfo_spaln"

    subprocess.call(command, shell=True)


def filterSpalnPairs(maxCoverage):
    """Keep only maxCoverage best proteins supporting each intron/start/stop.
       Kept pairs will be realigned with ProSplign.

    Args:
        maxCoverage (int): number of proteins to keep for each intron/start/stop
    """
    spalnDir = workDir + "/Spaln"
    os.chdir(spalnDir)

    command = binDir + "/select_best_proteins.py spaln.gff " + \
        str(maxCoverage) + " > pairs_filtered.out"
    subprocess.call(command, shell=True)


def prepareProSplignPairs(diamondPairs, k):
    """Prepare alignment pairs for ProSplign by combining filtered Spaln pairs
       and top k pairs from diamond pairs

    Args:
        diamondPairs (filePath): Path to file with seed gene-protein pairs to align
        k (int): Number of best pairs to use from DIAMOND output
    """
    os.chdir(workDir)

    command = binDir + "/combine_alignment_pairs.py  --spalnPairs Spaln/pairs_filtered.out  \
              --diamondPairs " + diamondPairs + " --out pairs_for_prosplign.out --k " + str(k)
    subprocess.call(command, shell=True)


def runProSplign(pbs):
    """Run ProSplign spliced alignment

    Args:
        closeThreshold (int): Percent identity threshold of ProSplign alignment
                              to be treated as an alignment of close homologs.
        pbs (boolean): Whether to run on pbs
    """
    proSplignDir = workDir + "/ProSplign"
    if not os.path.isdir(proSplignDir):
        os.mkdir(proSplignDir)
    os.chdir(proSplignDir)

    command = ""
    if not pbs:
        command = binDir + "/run_spliced_alignment.pl --cores " + threads + \
                  " --nuc ../nuc.fasta --list ../pairs_for_prosplign.out \
                  --prot " + proteins + " --v --aligner prosplign  > loginfo_prosplign"
    else:
        command = binDir + "/run_spliced_alignment_pbs.pl --N 240 --K 8 --seq \
                  ../nuc.fasta --list ../pairs_for_prosplign.out \
                  --db " + proteins + " --v --aligner prosplign  > loginfo_prosplign"

    subprocess.call(command, shell=True)


def processOutput(coverageThreshold, alignmentScoreThreshold, closeThreshold):
    """Prepare the final output
       Flag highly reliable introns from with high coverage and alignment score
       Create chained hints for close alignments

    Args:
        coverageThreshold (int): Minimum coverage for highly reliable introns
        alignmentScoreThreshold (float): Minimum alignment score for highly reliable introns
        closeThreshold (float): Percent identity threshold for close homologs (chained output)
    """
    os.chdir(workDir)

    # All introns
    command = binDir + "/combine_gff_records.pl --in_gff ProSplign/scored_introns.gff --out_gff introns.gff"
    subprocess.call(command, shell=True)

    # Add info about full protein aligned
    command = binDir + "/combine_outputs.py introns.gff ProSplign/prosplign.gff > tmp.gff; mv tmp.gff introns.gff"
    subprocess.call(command, shell=True)

    command = binDir + "/promapp2augustus.py --scoredIntrons introns.gff --al_score " + \
        str(alignmentScoreThreshold) + " --intronCoverage " + str(coverageThreshold) + " > hints.gff"
    subprocess.call(command, shell=True)

    # Start and stop codons
    command = binDir + "/combine_gff_records.pl --in_gff ProSplign/prosplign.gff --out_gff ProSplign/prosplign_combined.gff"
    subprocess.call(command, shell=True)

    # Add info about full protein aligned
    command = binDir + "/combine_outputs.py ProSplign/prosplign_combined.gff ProSplign/prosplign.gff > tmp.gff; \
              mv tmp.gff ProSplign/prosplign_combined.gff"
    subprocess.call(command, shell=True)

    command = binDir + "/promapp2augustus.py --startStops ProSplign/prosplign_combined.gff \
        --startCoverage " + str(coverageThreshold) + " --stopCoverage " + str(coverageThreshold) + " > ProSplign/starts_stops.gff"
    subprocess.call(command, shell=True)

    # Add stops to hints directly
    command = "grep -P \"\t[Ss]top_codon\t\" ProSplign/starts_stops.gff >> hints.gff"
    subprocess.call(command, shell=True)

    # Extra filtering step for starts -- filter by CDS overlaps
    command = "grep -P \"\t[Ss]tart_codon\t\" ProSplign/starts_stops.gff | sort -k1,1 -k4,4n -k5,5n > ProSplign/starts_sorted.gff"
    subprocess.call(command, shell=True)

    command = "grep -P \"\tCDS\t|\tCDS[Pp]art\t\" ProSplign/prosplign_combined.gff | sort -k1,1 -k4,4n -k5,5n > ProSplign/cds_sorted.gff"
    subprocess.call(command, shell=True)

    command = binDir + "/filterStarts.py ProSplign/starts_sorted.gff ProSplign/cds_sorted.gff 4 >> hints.gff"
    subprocess.call(command, shell=True)
    os.remove("ProSplign/cds_sorted.gff")
    os.remove("ProSplign/starts_sorted.gff")

    # Evidence file
    command = "grep \"src=M;\" hints.gff  > evidence.gff"
    subprocess.call(command, shell=True)

    # Chained hints
    # command = binDir + "/asn_to_gff.pl --asn ProSplign/prosplign.gff.regions.asn --out asn_regions.gff \
    #           --exons --exonCutoff 15 --close " + str(closeThreshold) + " --augustus"
    # subprocess.call(command, shell=True)

    # command = binDir + "/gff_from_region_to_contig.pl --in_gff asn_regions.gff --seq nuc.fasta --out chained.gff"
    # subprocess.call(command, shell=True)


def cleanup():
    print()


def setEnvironment(args):
    """Set up and check variables

    Args:
        args (dictionary): Command line arguments

    """
    global workDir, binDir, genome, proteins, threads
    workDir = os.path.abspath(args.workdir)
    binDir = os.path.abspath(os.path.dirname(__file__))

    genome = checkFileAndMakeAbsolute(args.genome)
    proteins = checkFileAndMakeAbsolute(args.proteins)

    if not os.path.isdir(workDir):
        os.mkdir(workDir)

    if args.geneMarkGtf:
        args.geneMarkGtf = checkFileAndMakeAbsolute(args.geneMarkGtf)

    if args.diamondPairs:
        args.diamondPairs = checkFileAndMakeAbsolute(args.diamondPairs)

    if args.threads > 0:
        threads = str(args.threads)
    else:
        threads = str(multiprocessing.cpu_count())


def checkFileAndMakeAbsolute(file):
    """Check if a file with given name exists and make the path to it absolute

    Args:
        file (filepath): Path to the file
    """
    if not os.path.isfile(file):
        sys.stderr.write("error: File \"" + file + "\" was not found.\n")
        sys.exit()
    return os.path.abspath(file)


def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser('promapp.py')

    parser.add_argument('genome', metavar='genome.fasta', type=str,
                        help='Input genomic sequence in FASTA format.')
    parser.add_argument('proteins', metavar='proteins.fasta', type=str,
                        help='Protein database in FASTA format.')

    parser.add_argument('--workdir', type=str, default='.',
                        help='Folder for results and temporary files. If not specified, current directory is used')
    parser.add_argument('--geneMarkGtf', type=str, default='',
                        help='File with GeneMark-ES predictions in gtf format. If this file is provided, GeneMark-ES run is skipped.')
    parser.add_argument('--diamondPairs', type=str, default='',
                        help='File with "seed gene-protein" hits generated by DIAMOND. If this file is provided, DIAMOND search for protein hits is skipped.')
    parser.add_argument('--diamondBin', type=str, default='',
                        help='Path to DIAMOND executable. If not provided, the program looks for the executable in the protein mapping bin folder.')
    parser.add_argument('--maxProteinsPerSeed', type=int, default=250,
                        help='Maximum number of protein hits per seed gene. Lowering this number decreases runtime but may decrease accuracy. \
                        Default is set to 250.')
    parser.add_argument('--evalue', type=float, default=0.001,
                        help='Maximum e-value for DIAMOND alignments hits. Default = 0.001')
    parser.add_argument('--maxSpalnCoverage', type=int, default=10,
                        help='Limit number of proteins supporting each intron and start/stop codon found by Spaln by this number. Lowering this number decreases runtime \
                        but may decrease accuracy. Default is set to 10.')
    parser.add_argument('--ensureDiamondPairs', type=int, default=5,
                        help='Always align (with ProSplign) this many top proteins for each gene in DIAMOND output (no matter the Spaln result). Default = 5.')
    parser.add_argument('--closeThreshold', type=int, default=90,
                        help='Percent identity threshold of ProSplign alignment to be treated as an alignment of close homologs. Full CDS structure \
                        including start and stop codons is reported for alignments with pID higher than this threshold.')
    parser.add_argument('--coverageThreshold', type=int, default=4,
                        help='Flag highly reliable introns with coverage >= coverageThreshold and alignment score >= alignmentScoreThreshold. Default = 4.')
    parser.add_argument('--alignmentScoreThreshold', type=float, default=0.3,
                        help='Flag highly reliable introns with coverage >= coverageThreshold and alignment score >= alignmentScoreThreshold. Default = 0.3.')
    parser.add_argument('--cleanup', default=False, action='store_true',
                        help='Delete temporary files and intermediate results. Cleanup is turned off by default as it is useful to keep these files \
                        for troubleshooting and the intermediate results might be useful on their own.')
    parser.add_argument('--pbs', default=False, action='store_true',
                        help='Run GeneMark-ES, Spaln and ProSplign on pbs.')
    parser.add_argument('--threads', type=int, default=-1,
                        help='Number of threads used by ES, DIAMOND, Spaln and ProSplign. By default, all available threads are used.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
