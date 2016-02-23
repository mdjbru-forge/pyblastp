### * Description

# Entry points for the command line script

DESCRIPTION = ("Convenient command line tool to run blastp locally. The blastp "
 "database can be generated on the fly from a fasta file, or an existing "
 "database can be used. Multiple cores can be used for the run (in this case "
 "the input file is split and several independent blastp runs are started "
 "instead of using the n_threads option of blastp, resulting in shorter run "
 "times). Note that blastp is for protein-protein blast. The output format is "
 "the tabular format from blastp (\"-outfmt 6\" option from blastp).")

### ** TODO

# Use the parse_known_args methods from the parser to send other options to
# blastp directly.

### * Set up

### ** Import

import os
import sys
import shutil
import argparse
import pyblastp as pyblastp

### * Parser

def makeParser() :

    """Build the parser
    
    Returns:
        ArgumentParser: An argument parser

    """
    parser = argparse.ArgumentParser(description = DESCRIPTION)
    # Query file
    parser.add_argument("query_file", type = str,
                        metavar = "QUERY_FILE",
                        help = "Fasta file containing the query sequences. "
                        "If no subject file is given and no database is "
                        "specified, then the "
                        "query file is also used as the subject file (i.e. "
                        "self-blastp on the query file).")
    # Subject file
    group = parser.add_mutually_exclusive_group()
    group.add_argument("subject_file", type = str,
                        metavar = "SUBJECT_FILE", nargs = "?",
                        help = "Fasta file containing the subject sequences")
    # Blastp task
    parser.add_argument("-t", "--task", choices = ["blastp", "blastp-fast",
                                                  "blastp-short"],
                        default = "blastp",
                        help = "Task to execute (default: blastp)")
    # Eval threshold
    parser.add_argument("-e", "--eval", metavar = "THRESHOLD",
                        type = float, default = 10,
                        help = "Expectation value (E) threshold for saving "
                        "hits (default: 10)")
    # Output file
    parser.add_argument("-o", "--output", metavar = "TABLE", type = str,
                        help = "Output file name (query file name with "
                        "\".out\" suffix if not specified)")
    # Max target seqs
    parser.add_argument("--max_target_seqs", metavar = "INT",
                        type = int, default = 10000,
                        help = "Maximum number of aligned sequences to keep "
                        "(default: 10000)")
    # Database
    group.add_argument("--db", metavar = "DATABASE", type = str,
                       help = "Use this database instead of generating a new "
                       "one from the subject file or the query file. This "
                       "option cannot be used if a subject file is specified.")
    # Parallel run with splitting of the input file
    parser.add_argument("-n", "--n_cores", metavar = "N_CORES", type = int,
                        default = 1, dest = "parallel",
                        help = "Number of cores to use for blastp. The query "
                        "file will be split in multiple fasta files and "
                        "several blastp runs will be started.")
    # Shuffle
    parser.add_argument("-s", "--shuffle", action = "store_true",
                        help = "Shuffle the sequences between the cores "
                        "(useful to make sure that all cores get a similar "
                        "work load in the case the sequences are sorted "
                        "by length in the input file)")
    # Reciprocal blastp
    parser.add_argument("-r", "--reciprocal_blastp", action = "store_true",
                        help = "Perform a reciprocal blastp between QUERY_FILE "
                        " and SUBJECT_FILE. Output to stdout.")
    return parser

### * Main

def main(args = None, stdout = None, stderr = None) :
    """Main entry point

    Args:
        args (namespace): Namespace with script arguments, parse the command 
          line arguments if None
        stdout (file): Writable stdout stream (if None, use `sys.stdout`)
        stderr (file): Writable stderr stream (if None, use `sys.stderr`)

    """
    if args is None :
        parser = makeParser()
        args = parser.parse_args()
    if stdout is None :
        stdout = sys.stdout
    if stderr is None :
        stderr = sys.stderr
    # Test for blast binaries
    # isAvailable
    # Function for one run
    def run(outfmt_option = "6") :
        # Ungap input file
        pyblastp.ungapFasta(args.query_file, args.query_file + ".ungap.tmp")
        inputFile = args.query_file + ".ungap.tmp"
        # Shuffle
        if args.shuffle :
            pyblastp.shuffleFastaFile(inputFile, inputFile + ".shuffled")
            shutil.move(inputFile + ".shuffled", inputFile)
        # Make database
        if args.subject_file is not None :
            db = pyblastp.makeBlastDb(inFile = args.subject_file,
                                      dbtype = "prot")["db"]
        elif args.db is None :
            db = pyblastp.makeBlastDb(inFile = inputFile, dbtype = "prot")["db"]
        else :
            db = args.db
        # Perform blastp
        if args.output is None :
            args.output = args.query_file + ".out"
        run = pyblastp.runBlastp(query = inputFile, db = db,
                                 evalMax = args.eval, task = args.task,
                                 out = args.output,
                                 max_target_seqs = args.max_target_seqs,
                                 cores = args.parallel,
                                 outfmt = outfmt_option)
        # Remove input file and database
        os.remove(inputFile)
        if args.db is None :
            pyblastp.removeBlastpDb(db = db)
    # Run
    if not args.reciprocal_blastp :
        run()
    # Reciprocal blastp
    if args.reciprocal_blastp :
        # Check that a subject file is provided
        assert args.subject_file is not None
        # Prepare the output file names
        qfile = os.path.basename(args.query_file)
        sfile = os.path.basename(args.subject_file)
        outFile1 = qfile + "-vs-" + sfile + ".tsv"
        outFile2 = sfile + "-vs-" + qfile + ".tsv"
        # First run
        args.output = outFile1
        run(outfmt_option = pyblastp.OUTFMT_RECIPROCAL_BLASTP_OPTION)
        # Second run
        args.query_file, args.subject_file = args.subject_file, args.query_file
        args.output = outFile2
        run(outfmt_option = pyblastp.OUTFMT_RECIPROCAL_BLASTP_OPTION)
        # Parse the run results
        out1 = pyblastp.getBestMatch(pyblastp.parseBlastpOutput(outFile1,
                 header = pyblastp.OUTFMT_RECIPROCAL_BLASTP_HEADER))
        out2 = pyblastp.getBestMatch(pyblastp.parseBlastpOutput(outFile2,
                 header = pyblastp.OUTFMT_RECIPROCAL_BLASTP_HEADER))
        bestMatches = pyblastp.reciprocalBestMatch(out1, out2)
        # Write to stdout
        if (len(bestMatches) > 0) :
            header = ["qseqid", "sseqid", "qlen", "slen", "qstart", "qend",
                      "sstart", "send", "length", "pident", "gapopen",
                      "mismatch", "evalue", "bitscore"]
            header = (["A", "B"] + [x + "-AB" for x in header] +
                      [x + "-BA" for x in header])
            stdout.write("\t".join(header) + "\n")
            for l in bestMatches :
                stdout.write("\t".join([l[x].strip() for x in header]) + "\n")
    sys.exit(0)
