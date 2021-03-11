#!/usr/bin/env python3
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
import tempfile


VERSION = '2.6.0'
workDir = ''
binDir = ''
genome = ''
proteins = ''
threads = ''

ProtHintRef = 'https://doi.org/10.1093/nargab/lqaa026'
DIAMONDRef = 'https://doi.org/10.1038/nmeth.3176'
SpalnRef = 'https://doi.org/10.1093/bioinformatics/btn460'


def main():
    args = parseCmd()

    setEnvironment(args)

    processInputProteins(args)

    if args.geneSeeds and args.prevGeneSeeds:
        nextIteration(args)
    else:
        standardRun(args)


def standardRun(args):
    """Execute a standard ProtHint run

    Args:
        args: Command line arguments
    """
    seedGenes = args.geneSeeds
    if not seedGenes:
        seedGenes = runGeneMarkES(args.pbs, args.fungus)
    else:
        sys.stderr.write("[" + time.ctime() + "] Skipping GeneMark-ES, using "
                         "the supplied gene seeds file instead\n")

    translateSeeds(seedGenes)

    diamondPairs = args.diamondPairs
    if not diamondPairs:
        diamondPairs = runDiamond(args.maxProteinsPerSeed, args.evalue)
    else:
        sys.stderr.write("[" + time.ctime() + "] Skipping DIAMOND, using "
                         "the supplied DIAMOND output file instead\n")

    prepareSeedSequences(diamondPairs)

    runSpaln(diamondPairs, args.pbs, args.minExonScore)

    checkOutputs(diamondPairs, seedGenes)
    flagTopProteins(diamondPairs)
    processSpalnOutput()

    if args.cleanup:
        cleanup()

    os.remove(proteins)
    sys.stderr.write("[" + time.ctime() + "] ProtHint finished.\n")


def nextIteration(args):
    """Run a next iteration of ProtHint. ProtHint is only run for gene
    seeds which are new or modified in the --geneSeeds file compared to
    --prevGeneSeeds. Hints for genes which are the same are reused from the
    --prevSpalnGff file but their gene_ids are updated to match the new
    seed file.

    Args:
        args: Command line arguments
    """
    sys.stderr.write("ProtHint is running in the iterative mode.\n")
    prepareDataForNextIteration(args.geneSeeds, args.prevGeneSeeds,
                                args.prevSpalnGff)

    diamondPairs = ""
    if os.path.getsize("uniqueSeeds.gtf") != 0:
        translateSeeds("uniqueSeeds.gtf")
        diamondPairs = runDiamond(args.maxProteinsPerSeed, args.evalue)
        prepareSeedSequences(diamondPairs)
        runSpaln(diamondPairs, args.pbs, args.minExonScore)
        flagTopProteins(diamondPairs)
        # Append subset of hints from the previous iteration to the current result
        os.chdir(workDir)
        with open("Spaln/spaln.gff", "a") as new:
            with open("prevHints.gff", "r") as prev:
                for line in prev:
                    new.write(line)
    else:
        sys.stderr.write("Warning: No unique gene seeds were detected in the " +
                         args.geneSeeds + " input file. ProtHint will only " \
                         "update seed gene IDs of hints to match the IDs in " \
                         "the new seed gene file.\n")
        if not os.path.isdir("Spaln"):
            os.mkdir("Spaln")
        shutil.move("prevHints.gff", "Spaln/spaln.gff")

    processSpalnOutput()

    os.remove(proteins)
    sys.stderr.write("[" + time.ctime() + "] ProtHint finished.\n")


def prepareDataForNextIteration(geneSeeds, prevGeneSeeds, prevSpalnGff):
    """Select gene seeds which are unique in this iteration and hints from
    previous iteration of ProtHint which correspond to seeds which are
    identical in both files (with gene ids updated to match the ids in the
    new seed file.

    Args:
        geneSeeds (filepath): Path to current gene seeds
        prevGeneSeeds (filepath): Path to previous genes seeds
        prevSpalnGff (filepath): Path to previous scored hints
    """
    sys.stderr.write("[" + time.ctime() + "] Selecting a subset of data to run"
                     " in the next iteration\n")

    os.chdir(workDir)

    callScript("print_longest_isoform.py", geneSeeds +
               " > longest_seed_isoforms.gtf")

    callScript("print_longest_isoform.py", prevGeneSeeds +
               " > longest_prevSeed_isoforms.gtf")

    callScript("select_for_next_iteration.py", "--geneSeeds " +
               "longest_seed_isoforms.gtf --prevGeneSeeds " +
               "longest_prevSeed_isoforms.gtf --prevSpalnGff " +
               prevSpalnGff + " --uniqueNewSeedsOut uniqueSeeds.gtf " +
               "--identicalSpalnGffOut prevHints.gff")
    os.remove("longest_prevSeed_isoforms.gtf")
    os.remove("longest_seed_isoforms.gtf")


def runGeneMarkES(pbs, fungus):
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

    pbsFlag = ""
    if pbs:
        pbsFlag = " --pbs"

    fungusFlag = ""
    if fungus:
        fungusFlag = " --fungus"

    callDependency("gmes_petap.pl", "--verbose --cores " + threads + pbsFlag +
                   " --ES --seq " + genome + " --soft auto" + fungusFlag,
                   "GeneMarkES")

    sys.stderr.write("[" + time.ctime() + "] GeneMark-ES finished.\n")
    return os.path.abspath("genemark.gtf")


def translateSeeds(seedGenes):
    """Translate GeneMark-ES seeds to proteins

    Args:
        seedGenes (filepath): Path to file with seed genes
    """
    sys.stderr.write("[" + time.ctime() + "] Translating gene seeds to " +
                     "proteins\n")
    os.chdir(workDir)

    callScript("print_longest_isoform.py", seedGenes +
               " > longest_seed_isoforms.gtf")

    systemCall("grep -P \tCDS\t longest_seed_isoforms.gtf > " +
               "longest_seed_isoforms_cds.gtf")
    os.remove("longest_seed_isoforms.gtf")

    callScript("proteins_from_gtf.pl", "--stat gene_stat.yaml --seq " +
               genome + " --annot longest_seed_isoforms_cds.gtf --out " +
               "seed_proteins.faa --format GTF")
    os.remove("longest_seed_isoforms_cds.gtf")

    sys.stderr.write("[" + time.ctime() + "] Translation of seeds finished\n")


def runDiamond(maxProteins, evalue):
    """Run DIAMOND protein search

    Args:
        maxProteins (int): Maximum number of protein hits per seed gene.
        evalue (float): Maximum e-value of DIAMOND hits

    Returns:
        string: Path to DIAMOND output
    """
    sys.stderr.write("[" + time.ctime() + "] Running DIAMOND\n")

    diamondDir = workDir + "/diamond"
    if not os.path.isdir(diamondDir):
        os.mkdir(diamondDir)
    os.chdir(diamondDir)

    # Make DIAMOND db
    callDependency("diamond", "makedb --in " + proteins + " -d diamond_db " +
                   "--threads " + threads)

    # Actual DIAMOND run
    callDependency("diamond", "blastp --query ../seed_proteins.faa --db " +
                   "diamond_db --outfmt 6 qseqid sseqid --out diamond.out " +
                   "--max-target-seqs " + str(maxProteins) + " --max-hsps 1 " +
                   "--threads " + threads + " --evalue " + str(evalue))

    sys.stderr.write("[" + time.ctime() + "] DIAMOND finished\n")
    return os.path.abspath("diamond.out")


def prepareSeedSequences(diamondPairs):
    """Prepare nucleotide sequences for seed genes

    Args:
        diamondPairs (filepath): Path to file with seed gene-protein pairs
    """
    sys.stderr.write("[" + time.ctime() + "] Preparing pairs for alignments\n")
    os.chdir(workDir)

    callScript("nucseq_for_selected_genes.pl", "--seq " + genome + " --out " +
               "nuc.fasta --gene gene_stat.yaml --list " + diamondPairs)

    sys.stderr.write("[" + time.ctime() + "] Preparation of pairs finished\n")


def runSpaln(diamondPairs, pbs, minExonScore):
    """Run Spaln spliced alignment and score the outputs with spaln-boundary-scorer

    Args:
        diamondPairs (filePath): Path to file with seed gene-protein pairs
                                 to align
        pbs (boolean): Whether to run on pbs
        minExonScore (float): Discard all hints inside/neighboring exons with
                              score lower than minExonScore
    """
    spalnDir = workDir + "/Spaln"
    if not os.path.isdir(spalnDir):
        os.mkdir(spalnDir)
    os.chdir(spalnDir)

    if not pbs:
        callScript("run_spliced_alignment.pl", "--cores " + threads +
                   " --nuc ../nuc.fasta --list " + diamondPairs + " --prot " +
                   proteins + " --v --aligner spaln --min_exon_score " +
                   str(minExonScore))
    else:
        callScript("run_spliced_alignment_pbs.pl", "--N 120 --K " + threads +
                   " --seq ../nuc.fasta --list " + diamondPairs + " --db " +
                   proteins + " --v --aligner spaln --min_exon_score " +
                   str(minExonScore))


def checkOutputs(diamondPairs, seedGenes):
    """Check whether all intermediate outputs were correctly created

    Args:
        diamondPairs (filepath): Path to file with seed gene-protein pairs
    """
    os.chdir(workDir)

    msg = 'This error can be caused by:\n' \
          '    a) The set of input proteins is too small and/or the ' \
          'proteins are too remote.\n' \
          '    b) The gene seeds identified by GeneMark-ES (or the ' \
          'supplied gene seeds in case the option "--geneSeeds" was ' \
          'used) are incorrect. Please try running GeneMark-ES ' \
          'separately to identify errors related to gene seeds (' \
          'https://github.com/gatech-genemark/ProtHint#genemark-es).'

    if os.stat(diamondPairs).st_size == 0:
        sys.exit('error: No homologous proteins were found by DIAMOND (' +
                 diamondPairs + ' is empty).\nThis error can be caused' + msg)

    if os.stat("Spaln/spaln.gff").st_size == 0:
        sys.exit('error: No spliced alignments were created by Spaln (' +
                 workDir + '/Spaln/spaln.gff is empty).\n' + msg)


def flagTopProteins(diamondPairs):
    """Label hints which were mapped from the best DIAMOND target

    Args:
        diamondPairs (filepath): Path to file with seed gene-protein pairs
    """
    sys.stderr.write("[" + time.ctime() + "] Flagging top chains\n")
    os.chdir(workDir)

    callScript("flag_top_proteins.py", "Spaln/spaln.gff " + diamondPairs +
               " > tmp")
    shutil.move("tmp", "Spaln/spaln.gff")


def processSpalnOutput():
    """Prepare the final output from Spaln result scored by spaln-boundary-scorer
       Convert the output to GeneMark and Augustus compatible formats

    """
    sys.stderr.write("[" + time.ctime() + "] Processing the output\n")
    os.chdir(workDir)

    processSpalnIntrons()
    processSpalnStops()
    processSpalnStarts()
    systemCall("sort -k1,1 -k4,4n -k5,5n -o prothint.gff prothint.gff")

    printTopChains()

    # High confidence
    callScript("print_high_confidence.py", "prothint.gff > evidence.gff")

    # Augustus compatible format
    callScript("prothint2augustus.py", "prothint.gff evidence.gff "
               "top_chains.gff prothint_augustus.gff")
    sys.stderr.write("[" + time.ctime() + "] Output processed\n")


def processSpalnIntrons():
    systemCall("grep Intron Spaln/spaln.gff > introns.gff || [[ $? == 1 ]]")

    # Filter out introns with alignment score < 0.1
    callScript("print_high_confidence.py", "introns.gff --intronCoverage 0 " +
               "--intronAlignment 0.1 --addAllSpliceSites > introns_01.gff")
    os.remove("introns.gff")

    callScript("combine_gff_records.pl", "--in_gff introns_01.gff --out_gff " +
               "prothint.gff")
    os.remove("introns_01.gff")


def processSpalnStops():
    systemCall("grep stop_codon Spaln/spaln.gff > stops.gff || [[ $? == 1 ]]")

    # Filter out stops with alignment score < 0.01
    callScript("print_high_confidence.py", "stops.gff --stopCoverage 0 " +
               "--stopAlignment 0.01 > stops_01.gff")
    os.remove("stops.gff")

    callScript("combine_gff_records.pl", "--in_gff stops_01.gff --out_gff " +
               "stops_01_combined.gff")
    systemCall("cat stops_01_combined.gff >> prothint.gff")

    os.remove("stops_01.gff")
    os.remove("stops_01_combined.gff")


def processSpalnStarts():
    systemCall("grep start_codon Spaln/spaln.gff > starts.gff || [[ $? == 1 ]]")

    # Filter out starts with alignment score < 0.01
    callScript("print_high_confidence.py", "starts.gff --startCoverage 0 " +
               "--startAlignment 0.01 > starts_01.gff")
    os.remove("starts.gff")

    callScript("combine_gff_records.pl", "--in_gff starts_01.gff --out_gff " +
               "starts_01_combined.gff")
    os.remove("starts_01.gff")

    # The rest of this function counts CDS overlap of starts

    systemCall("sort -k1,1 -k4,4n -k5,5n starts_01_combined.gff > " +
               "starts_01_combined_sorted.gff")
    os.remove("starts_01_combined.gff")

    systemCall("grep CDS Spaln/spaln.gff > cds.gff || [[ $? == 1 ]]")
    callScript("combine_gff_records.pl", "--in_gff cds.gff --out_gff " +
               "cds_combined.gff")
    os.remove("cds.gff")

    # Only count CDS regions which have an upstream support
    # (by start codon or intron) in hints
    callScript("cds_with_upstream_support.py", "cds_combined.gff " +
               "starts_01_combined_sorted.gff prothint.gff > tmp")
    shutil.move("tmp", "cds_combined.gff")

    systemCall("sort -k1,1 -k4,4n -k5,5n cds_combined.gff > " +
               "cds_combined_sorted.gff")
    os.remove("cds_combined.gff")

    callScript("count_cds_overlaps.py", "starts_01_combined_sorted.gff " +
               "cds_combined_sorted.gff >> prothint.gff")

    os.remove("starts_01_combined_sorted.gff")
    os.remove("cds_combined_sorted.gff")


def printTopChains():
    systemCall("grep topProt=TRUE Spaln/spaln.gff > topProteins.gff " +
               "|| [[ $? == 1 ]]")

    if os.stat("topProteins.gff").st_size == 0:
        sys.exit('error: The "topProt=TRUE" flag is missing in the '
                 'Spaln/spaln.gff output file. This issue can be caused by '
                 'the presence of special characters in the fasta headers of '
                 'input files. Please remove any special characters and '
                 're-run ProtHint. See https://github.com/gatech-genemark/ProtHint#input '
                 'for more details about the input format.')

    callScript("print_high_confidence.py", "topProteins.gff --startCoverage 0 " +
               "--startAlignment 0.01 --stopCoverage 0 --stopAlignment 0.01 " +
               "--intronCoverage 0 --intronAlignment 0.1 --addAllSpliceSites" +
               " > topProteinsFiltered.gff")
    os.remove("topProteins.gff")

    callScript("make_chains.py", "topProteinsFiltered.gff > top_chains.gff")
    os.remove("topProteinsFiltered.gff")


def cleanup():
    """Delete temporary files and intermediate results
    """
    sys.stderr.write("[" + time.ctime() + "] Cleaning up\n")
    os.chdir(workDir)

    try:
        os.remove("gene_stat.yaml")
        os.remove("seed_proteins.faa")
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

    if os.path.exists("GeneMark_ES"):
        shutil.rmtree("GeneMark_ES/data", ignore_errors=True)
        shutil.rmtree("GeneMark_ES/info", ignore_errors=True)
        shutil.rmtree("GeneMark_ES/run", ignore_errors=True)
        shutil.rmtree("GeneMark_ES/output/data", ignore_errors=True)
        shutil.rmtree("GeneMark_ES/output/gmhmm", ignore_errors=True)


def processInputProteins(args):
    """Remove dots from the input file with proteins.
       OrhoDB protein sequences sometimes end with a dot. This format is not
       compatible with DIAMOND.
       Clean fasta headers by removing pipe ("|") characters.
    """
    global proteins
    sys.stderr.write("[" + time.ctime() + "] Pre-processing protein input\n")
    os.chdir(workDir)
    protFile = tempfile.NamedTemporaryFile(delete=False, dir='.', prefix="prot")
    systemCall('sed \"s/\.//\" ' + args.proteins + ' | sed \"s/|/_/g\" > ' +
               protFile.name)
    proteins = checkFileAndMakeAbsolute(protFile.name)


def setEnvironment(args):
    """Set up and check variables

    Args:
        args (dictionary): Command line arguments

    """
    sys.stderr.write("ProtHint Version " + VERSION + "\n")
    sys.stderr.write("Copyright 2019, Georgia Institute of Technology, USA\n\n")
    sys.stderr.write("Please cite\n")
    sys.stderr.write("  - ProtHint: " + ProtHintRef + "\n")
    sys.stderr.write("  - DIAMOND:  " + DIAMONDRef + "\n")
    sys.stderr.write("  - Spaln:    " + SpalnRef + "\n\n")
    global workDir, binDir, genome, threads
    workDir = os.path.abspath(args.workdir)
    binDir = os.path.abspath(os.path.dirname(__file__))

    genome = checkFileAndMakeAbsolute(args.genome)
    args.proteins = checkFileAndMakeAbsolute(args.proteins)

    if args.ProSplign:
        sys.exit("error: ProSplign is not supported in this version of ProtHint. "
                 "For running ProtHint with ProSplign, use ProtHint release v2.4.0: "
                 "https://github.com/gatech-genemark/ProtHint/releases/tag/v2.4.0")

    if args.geneMarkGtf:
        if args.geneSeeds:
            sys.exit("error: please specify either --geneSeeds or\n"
                     "--geneMarkGtf. The arguments are identical,\n"
                     "--geneMarkGtf is supported for backwards compatibility.")
        args.geneSeeds = args.geneMarkGtf

    if args.geneSeeds:
        args.geneSeeds = checkFileAndMakeAbsolute(args.geneSeeds)

    if args.prevGeneSeeds:
        if not args.prevSpalnGff or not args.prevSpalnGff:
            sys.exit("error: --prevSpalnGff and --geneSeeds must be given when\n"
                     "using --prevGeneSeeds")
        args.prevGeneSeeds = checkFileAndMakeAbsolute(args.prevGeneSeeds)

    if args.prevSpalnGff:
        if not args.geneSeeds or not args.prevGeneSeeds:
            sys.exit("error: --geneSeeds and --prevGeneSeeds options must\n"
                     "be given when using --prevSpalnGff option")
        args.prevSpalnGff = checkFileAndMakeAbsolute(args.prevSpalnGff)

    if args.diamondPairs:
        if not args.geneSeeds:
            sys.exit("error: File with gene predictions (--geneSeeds\n"
                     "option) not specified. When using the --diamondPairs\n"
                     "option, a prediction file with seed genes corresponding\n"
                     "to seed genes in the DIAMOND pairs file must be specified.")
        args.diamondPairs = checkFileAndMakeAbsolute(args.diamondPairs)

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


def systemCall(cmd):
    """Make a system call

    Args:
        cmd (string): Command to call
    """
    sys.stderr.flush()
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit('[' + time.ctime() + '] error: ProtHint exited due to an ' +
                 'error in command: ' + cmd)


def callScript(name, args):
    """Call a script located in the ProtHint bin folder

    Args:
        name (string): Name of the script
        args (string): Command line arguments to use in the call
    """
    systemCall(binDir + '/' + name + ' ' + args)


def callDependency(name, args, location=''):
    """Call a dependency. ProtHint first looks for the dependency in the
    dependencies folder. If not found there, it tries to find it in the path.

    Args:
        name (string): Name of the dependency
        args (string): Command line arguments to use in the call
        location (string): Location of the dependency within the dependencies
                           folder
    """
    if location != '':
        location = location + '/'

    if os.path.isfile(binDir + '/../dependencies/' + location + name):
        systemCall(binDir + '/../dependencies/' + location + name + ' ' + args)
    else:
        sys.stderr.write('[' + time.ctime() + '] warning: Could not find ' +
                         name + ' in dependencies/' + location + ' folder.' +
                         ' Attempting to use ' + name + ' in the PATH.\n')
        if shutil.which(name) is not None:
            systemCall(name + ' ' + args)
        else:
            sys.exit('[' + time.ctime() + '] error: Could not find ' + name +
                     ' in the PATH')


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
    parser.add_argument('--geneSeeds', type=str, default='',
                        help='Gene prediction seeds in gtf format. If this file is provided, GeneMark-ES run is skipped.')
    parser.add_argument('--geneMarkGtf', type=str, default='',
                        help='Same as --geneSeeds, for backwards compatibility')
    parser.add_argument('--fungus', default=False, action='store_true',
                        help='Run GeneMark-ES in fungus mode.')
    parser.add_argument('--diamondPairs', type=str, default='',
                        help='File with "seed gene-protein" hits generated by DIAMOND. If this file is provided, DIAMOND search for protein hits is skipped.\
                        The seed genes in this file must correspond to seed genes passed by "--geneSeeds" option. All pairs in the file are used -- option \
                        "--maxProteinsPerSeed" is ignored.')
    parser.add_argument('--maxProteinsPerSeed', type=int, default=25,
                        help='Maximum number of protein hits per seed gene. Increasing this number leads to increased runtime and may improve the\
                        sensitivity of hints. Decreasing has an opposite effect. Default is set to 25.')
    parser.add_argument('--evalue', type=float, default=0.001,
                        help='Maximum e-value for DIAMOND alignments hits. Default = 0.001')
    parser.add_argument('--minExonScore', type=float, default=25,
                        help='Discard all hints inside/neighboring exons with score lower than minExonScore. Default = 25')
    parser.add_argument('--cleanup', default=False, action='store_true',
                        help='Delete temporary files and intermediate results. Cleanup is turned off by default as it is useful to keep these files \
                        for troubleshooting and the intermediate results might be useful on their own.')
    parser.add_argument('--pbs', default=False, action='store_true',
                        help='Run GeneMark-ES and Spaln on pbs.')
    parser.add_argument('--threads', type=int, default=-1,
                        help='Number of threads used by ES, DIAMOND, and Spaln. By default, all available threads are used.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    parser.add_argument('--prevGeneSeeds', type=str, default='',
                        help='File with gene seeds which were used in the previous iteration. Next iteration of ProtHint is only executed for --geneSeeds which \
                        differ from --prevGeneSeeds. --prevSpalnGff is required when this option is used since results from the previous iteration are reused \
                        for seeds which do not differ (Gene ids of such hints are updated to match the new seed genes).')
    parser.add_argument('--prevSpalnGff', type=str, default='',
                        help='Scored hints from previous iteration.')

    parser.add_argument('--ProSplign',  default=False, action='store_true',
                        help=argparse.SUPPRESS)

    return parser.parse_args()


if __name__ == '__main__':
    main()
