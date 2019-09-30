#!/usr/bin/env python
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# ProtHint: Pipeline for generating genome wide footprints of homologous
# proteins
# ==============================================================

import argparse
import os
import sys
import subprocess
import multiprocessing
import time
import shutil


VERSION = '2.1.1'
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
    else:
        sys.stderr.write("[" + time.ctime() + "] Skipping GeneMark-ES, using "
                         "the supplied genemark.gtf file instead\n")

    translateSeeds(geneMarkGtf)

    diamondPairs = args.diamondPairs
    if not diamondPairs:
        diamondPairs = runDiamond(args.diamondBin, args.maxProteinsPerSeed, args.evalue)
    else:
        sys.stderr.write("[" + time.ctime() + "] Skipping DIAMOND, using "
                         "the supplied DIAMOND output file instead\n")

    prepareSeedSequences(diamondPairs)
    runSpaln(diamondPairs, args.pbs, args.minExonScore)

    if (not args.ProSplign):
        processSpalnOutput()
    else:
        filterSpalnPairs(args.maxSpalnCoverage)
        prepareProSplignPairs(diamondPairs, args.ensureDiamondPairs)
        runProSplign(args.pbs)
        processProSplignOutput()

    if args.cleanup:
        cleanup()

    sys.stderr.write("[" + time.ctime() + "] ProtHint finished.\n")


def runGeneMarkES(pbs):
    """Run GeneMark-ES

    Args:
        pbs (boolean): Whether to run on pbs

    Returns:
        string: Path to genemark gff output
    """
    sys.stderr.write("[" + time.ctime() + "] Running GeneMark-ES.\n")
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
    sys.stderr.write("[" + time.ctime() + "] GeneMark-ES finished.\n")
    return os.path.abspath("genemark.gtf")


def translateSeeds(geneMarkGtf):
    """Translate GeneMark-ES seeds to proteins

    Args:
        geneMarkGtf (filepath): Path to GeneMark.gtf prediction file
    """
    sys.stderr.write("[" + time.ctime() + "] Translating GeneMark seeds to proteins\n")
    os.chdir(workDir)

    command = binDir + "/proteins_from_gtf.pl --stat gene_stat.yaml --seq " + \
              genome + " --annot " + geneMarkGtf + " --out gmes_proteins.faa --format GTF"
    subprocess.call(command, shell=True)
    sys.stderr.write("[" + time.ctime() + "] Translation of seeds finished\n")


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
    sys.stderr.write("[" + time.ctime() + "] Running DIAMOND\n")
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

    sys.stderr.write("[" + time.ctime() + "] DIAMOND finished\n")
    return os.path.abspath("diamond.out")


def prepareSeedSequences(diamondPairs):
    """Prepare nucleotide sequences for seed genes

    Args:
        diamondPairs (filepath): Path to file with seed gene-protein pairs
    """
    sys.stderr.write("[" + time.ctime() + "] Preparing pairs for alignments\n")
    os.chdir(workDir)

    command = binDir + "/nucseq_for_selected_genes.pl --seq " + genome + \
              " --out nuc.fasta --gene gene_stat.yaml --list " + diamondPairs
    subprocess.call(command, shell=True)
    sys.stderr.write("[" + time.ctime() + "] Preparation of pairs finished\n")


def runSpaln(diamondPairs, pbs, minExonScore):
    """Run Spaln spliced alignment and score the outputs with spaln-boundary-scorer

    Args:
        diamondPairs (filePath): Path to file with seed gene-protein pairs to align
        pbs (boolean): Whether to run on pbs
        minExonScore (float): Discard all hints inside/neighboring exons with
                              score lower than minExonScore
    """
    spalnDir = workDir + "/Spaln"
    if not os.path.isdir(spalnDir):
        os.mkdir(spalnDir)
    os.chdir(spalnDir)

    command = ""
    if not pbs:
        command = binDir + "/run_spliced_alignment.pl --cores " + threads + \
                  " --nuc ../nuc.fasta --list " + diamondPairs + \
                  " --prot " + proteins + " --v --aligner spaln " + \
                  "--min_exon_score " + str(minExonScore)

    else:
        command = binDir + "/run_spliced_alignment_pbs.pl --N 120 --K " + threads + \
                  " --seq ../nuc.fasta --list " + diamondPairs + " --db " + \
                  proteins + " --v --aligner spaln " + \
                  "--min_exon_score " + str(minExonScore)

    subprocess.call(command, shell=True)


def processSpalnOutput():
    """Prepare the final output from Spaln result scored by spaln-boundary-scorer
       Convert the output to GeneMark and Augustus compatible formats
    """
    sys.stderr.write("[" + time.ctime() + "] Processing the output\n")
    os.chdir(workDir)

    processSpalnIntrons()
    processSpalnStops()
    processSpalnStarts()

    # High confidence
    subprocess.call(binDir + "/print_high_confidence.py prothint.gff > evidence.gff", shell=True)

    # Augustus compatible format
    subprocess.call(binDir + "/prothint2augustus.py prothint.gff > prothint_augustus.gff", shell=True)
    subprocess.call(binDir + "/prothint2augustus.py evidence.gff > evidence_augustus.gff", shell=True)

    sys.stderr.write("[" + time.ctime() + "] Output processed\n")


def processSpalnIntrons():
    subprocess.call("grep Intron Spaln/spaln.gff > introns.gff", shell=True)

    subprocess.call(binDir + "/print_high_confidence.py introns.gff "
                    "--intronCoverage 0 --intronAlignment 0.1 "
                    "--addAllSpliceSites > introns_01.gff; rm introns.gff", shell=True)

    subprocess.call(binDir + "/combine_gff_records.pl --in_gff introns_01.gff "
                    "--out_gff prothint.gff; rm introns_01.gff", shell=True)


def processSpalnStops():
    subprocess.call("grep stop_codon Spaln/spaln.gff > stops.gff", shell=True)

    subprocess.call(binDir + "/print_high_confidence.py stops.gff "
                    "--stopCoverage 0 --stopScore 0.01 > stops_01.gff;"
                    "rm stops.gff", shell=True)

    subprocess.call(binDir + "/combine_gff_records.pl --in_gff stops_01.gff "
                    "--out_gff stops_01_combined.gff; cat stops_01_combined.gff >> "
                    "prothint.gff; rm stops_01.gff stops_01_combined.gff", shell=True)


def processSpalnStarts():
    subprocess.call("grep start_codon Spaln/spaln.gff > starts.gff", shell=True)

    subprocess.call(binDir + "/print_high_confidence.py starts.gff "
                    "--startCoverage 0 --startScore 0.01 > starts_01.gff;"
                    "rm starts.gff", shell=True)

    subprocess.call(binDir + "/combine_gff_records.pl --in_gff starts_01.gff "
                    "--out_gff starts_01_combined.gff; rm starts_01.gff", shell=True)

    subprocess.call("sort -k1,1 -k4,4n -k5,5n starts_01_combined.gff > "
                    "starts_01_combined_sorted.gff; rm starts_01_combined.gff", shell=True)

    # Count CDS overlap of starts
    subprocess.call("grep CDS Spaln/spaln.gff > cds.gff", shell=True)

    subprocess.call(binDir + "/combine_gff_records.pl --in_gff cds.gff "
                    "--out_gff cds_combined.gff; rm cds.gff", shell=True)

    # Only count CDS regions which have an upstream support
    # (by start codon or intron) in hints
    subprocess.call(binDir + "/cds_with_upstream_support.py cds_combined.gff "
                    "starts_01_combined_sorted.gff prothint.gff > tmp; \
                    mv tmp cds_combined.gff", shell=True)

    subprocess.call("sort -k1,1 -k4,4n -k5,5n cds_combined.gff > "
                    "cds_combined_sorted.gff; rm cds_combined.gff", shell=True)

    subprocess.call(binDir + "/count_cds_overlaps.py starts_01_combined_sorted.gff "
                    "cds_combined_sorted.gff >> prothint.gff;"
                    "rm starts_01_combined_sorted.gff cds_combined_sorted.gff", shell=True)


def filterSpalnPairs(maxCoverage):
    """Keep only maxCoverage best proteins supporting each intron/start/stop.
       Kept pairs will be realigned with ProSplign.

    Args:
        maxCoverage (int): number of proteins to keep for each intron/start/stop
    """
    sys.stderr.write("[" + time.ctime() + "] Selecting proteins for ProSplign alignment.\n")
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
    """Run ProSplign spliced alignment and score the outputs with prosplign-intron-scorer

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
                  --prot " + proteins + " --v --aligner prosplign"
    else:
        command = binDir + "/run_spliced_alignment_pbs.pl --N 240 --K " + threads + \
                  " --seq ../nuc.fasta --list ../pairs_for_prosplign.out \
                  --db " + proteins + " --v --aligner prosplign"

    subprocess.call(command, shell=True)


def processProSplignOutput():
    """Prepare the final output from full ProSplign output and ProSplign
       introns scored by prosplign-intron-scorer.
       Convert the output to GeneMark and Augustus compatible formats
    """
    sys.stderr.write("[" + time.ctime() + "] Processing the output\n")
    os.chdir(workDir)

    # Collapse all scored introns, add them to the final output
    subprocess.call(binDir + "/combine_gff_records.pl --in_gff \
        ProSplign/scored_introns.gff --out_gff prothint.gff", shell=True)

    # Collapse full ProSplign output
    subprocess.call(binDir + "/combine_gff_records.pl --in_gff ProSplign/prosplign.gff \
                    --out_gff ProSplign/prosplign_combined.gff", shell=True)

    # Add stops to output directly
    subprocess.call("grep -P \"\t[Ss]top_codon\t\" ProSplign/prosplign_combined.gff >> \
                    prothint.gff", shell=True)

    # Count CDS overlaps of starts before adding them to the output file
    subprocess.call("grep -P \"\t[Ss]tart_codon\t\" ProSplign/prosplign_combined.gff | \
                    sort -k1,1 -k4,4n -k5,5n > ProSplign/starts_sorted.gff", shell=True)

    subprocess.call("grep -P \"\tCDS\t|\tCDS[Pp]art\t\" ProSplign/prosplign_combined.gff | \
                    sort -k1,1 -k4,4n -k5,5n > ProSplign/cds_sorted.gff", shell=True)

    subprocess.call(binDir + "/count_cds_overlaps.py ProSplign/starts_sorted.gff \
                    ProSplign/cds_sorted.gff >> prothint.gff", shell=True)
    os.remove("ProSplign/cds_sorted.gff")
    os.remove("ProSplign/starts_sorted.gff")

    # Add info about full protein aligned
    subprocess.call(binDir + "/combine_outputs.py prothint.gff ProSplign/prosplign.gff > \
                    tmp.gff; mv tmp.gff prothint.gff", shell=True)

    # Print high confidence hints
    subprocess.call(binDir + "/print_high_confidence.py --startCoverage 4 "
                    "--startOverlap 3 prothint.gff > evidence.gff", shell=True)

    # Augustus compatible format
    subprocess.call(binDir + "/prothint2augustus.py prothint.gff > prothint_augustus.gff", shell=True)
    subprocess.call(binDir + "/prothint2augustus.py evidence.gff > evidence_augustus.gff", shell=True)

    sys.stderr.write("[" + time.ctime() + "] Output processed\n")


def cleanup():
    """Delete temporary files and intermediate results
    """
    sys.stderr.write("[" + time.ctime() + "] Cleaning up\n")
    os.chdir(workDir)
    try:
        os.remove("gene_stat.yaml")
        os.remove("gmes_proteins.faa")
        os.remove("nuc.fasta")
    except OSError:
        pass

    try:
        os.remove("diamond/diamond_db.dmnd")
    except OSError:
        pass

    try:
        os.remove("Spaln/spaln.gff")
    except OSError:
        pass

    try:
        os.remove("ProSplign/scored_introns.gff")
        os.remove("ProSplign/prosplign.gff")
    except OSError:
        pass

    if os.path.exists("GeneMark_ES"):
        shutil.rmtree("GeneMark_ES/data", ignore_errors=True)
        shutil.rmtree("GeneMark_ES/info", ignore_errors=True)
        shutil.rmtree("GeneMark_ES/run", ignore_errors=True)
        shutil.rmtree("GeneMark_ES/output/data", ignore_errors=True)
        shutil.rmtree("GeneMark_ES/output/gmhmm", ignore_errors=True)


def setEnvironment(args):
    """Set up and check variables

    Args:
        args (dictionary): Command line arguments

    """
    sys.stderr.write("ProtHint Version " + VERSION + "\n")
    sys.stderr.write("Copyright 2019, Georgia Institute of Technology, USA\n\n")
    global workDir, binDir, genome, proteins, threads
    workDir = os.path.abspath(args.workdir)
    binDir = os.path.abspath(os.path.dirname(__file__))

    genome = checkFileAndMakeAbsolute(args.genome)
    proteins = checkFileAndMakeAbsolute(args.proteins)

    if not os.path.isdir(workDir):
        os.mkdir(workDir)

    # Log info about cmd
    callDir = "Called from: " + os.path.abspath(".") + "\n"
    cmd = "Cmd: " + " ".join(sys.argv) + "\n\n"
    sys.stderr.write(callDir)
    sys.stderr.write(cmd)
    with open(workDir + "/cmd.log", "w") as file:
        file.write(callDir)
        file.write(cmd)

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
    parser = argparse.ArgumentParser(description='ProtHint ' + VERSION + ': Pipeline for generating genome wide \
                                     footprints of homologous proteins. The set of high confidence hints \
                                     is generated using default thresholds in print_high_confidence.py \
                                     script. If you wish to use different filtering criteria, re-run\
                                     print_high_confidence.py script with custom thresholds.')

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
    parser.add_argument('--minExonScore', type=float, default=25,
                        help='Discard all hints inside/neighboring exons with score lower than minExonScore. Default = 25')
    parser.add_argument('--cleanup', default=False, action='store_true',
                        help='Delete temporary files and intermediate results. Cleanup is turned off by default as it is useful to keep these files \
                        for troubleshooting and the intermediate results might be useful on their own.')
    parser.add_argument('--pbs', default=False, action='store_true',
                        help='Run GeneMark-ES, Spaln and ProSplign on pbs.')
    parser.add_argument('--threads', type=int, default=-1,
                        help='Number of threads used by ES, DIAMOND, Spaln and ProSplign. By default, all available threads are used.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    parser.add_argument('--ProSplign',  default=False, action='store_true',
                        help='Re-align hints discovered by Spaln with ProSplign and report ProSplign results. This is a legacy option which makes the pipeline \
                        considerably (~20x) slower without a significant effect on the final accuracy. NOT SUPPORTED IN THIS VERSION.')
    parser.add_argument('--maxSpalnCoverage', type=int, default=10,
                        help='For each hint, select maxSpalnCoverage proteins for ProSplign realignment.')
    parser.add_argument('--ensureDiamondPairs', type=int, default=5,
                        help='Always align (with ProSplign) this many top proteins for each gene in DIAMOND output (no matter the Spaln result). Default = 5.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
