### * Description

# Python wrapper for blastp

### * Setup

### ** Import

import os
import sys
import subprocess
import argparse
import random
import math
from Bio import SeqIO
from Bio import AlignIO

### ** Parameters

OUTFMT_HEADER = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

### * Functions

### ** randomTag(n)

def randomTag(n) :
    """Generate a random tag of size n.

    Args:
        n (int): Length of the tag

    Returns:
        str: A tag of length n, which each character in [0-9][a-f]
    """
    o = ""
    allowed = "0123456789abcdef"
    for x in range(n) :
        o += random.choice(allowed)
    return o

### ** makeBlastDb(inFile, dbtype, out = None)

def makeBlastDb(inFile, dbtype, outDb = None) :
    """Prepare a blast database by calling 'makeblastdb'

    Args:
        inFile (str): Input fasta file to prepare the database
        dbtype (str): Database type ('prot' or 'nucl')
        outDb (str): Database file name (if None, a random name is produced)

    Returns:
        dict: Dictionary with
    
        - 'db': database name
        - 'output': output from 'makeblastdb'
        - 'error': stderr from 'makeblastdb'

    """
    command = ["makeblastdb"]
    command += ["-in", inFile]
    command += ["-dbtype", dbtype]
    if outDb is None :
        outDb = "blastp.db." + randomTag(16)
    command += ["-out", outDb]
    p = subprocess.Popen(command, stdout = subprocess.PIPE,
                         stderr = subprocess.PIPE)
    (output, error) = p.communicate()
    p.wait()
    return {"db" : outDb, "output": output, "error": error}

### ** runBlastp(query, db, evalMax, task, out, max_target_seqs, jobs, cores)

def runBlastp(query, db, evalMax, task, out, max_target_seqs = 10000, jobs = 1, cores = 1) :
    """Perform a blastp run. For the parallel version of it, use GNU parallel, 
    as explained here: https://www.biostars.org/p/63816/.

    Args:
        query (str): Name of the query fasta file
        db (str): Name of the subject database (`db` element in the output 
          from :func:`makeBlastDb`)
        evalMax (float): Evalue threshold for saving hits
        task (str): Task to execute (`blastp`, `blastp-fast` or 
          `blastp-short`)
        out (str): Name of the output file. If `None`, use the input file name 
          with '.out' suffix appended.
        jobs (int): Number of cores to use for blastp run. If more than 1, use
          GNU parallel to split the workload across several cores.

    Returns:
        Name of the output file

    """
    assert jobs > 0 and cores > 0
    if out is None :
        out = query + ".out"
    if jobs == 1 and cores == 1:
        command = ["blastp"]
        command += ["-db", db]
        command += ["-query", query]
        command += ["-task", task]
        command += ["-seg", "yes"]
        command += ["-evalue", str(evalMax)]
        command += ["-out", out]
        command += ["-outfmt", "6"]
        command += ["-max_target_seqs", str(max_target_seqs)]
        p = subprocess.Popen(command)
        p.wait()
    elif cores == 1 :
        command = ""
        command += "cat " + query + " | "
        command += "parallel --block 100k --recstart '>' --jobs " + str(jobs)
        command += " --pipe "
        command += "blastp -query - "
        command += "-db " + db + " -task " + task
        command += " -seg yes " + "-evalue " + str(evalMax)
        command += " -outfmt 6 "
        command += "-max_target_seqs " + str(max_target_seqs) + " "
        command += "> " + out
        os.system(command)
    elif jobs == 1:
        # Split the input file into subsets
        inputFiles = splitFastaFile(query, cores)
        # Run blastp for each subset
        blastpProcesses = []
        for inputFile in inputFiles :
            command = ["blastp"]
            command += ["-db", db]
            command += ["-query", inputFile]
            command += ["-task", task]
            command += ["-seg", "yes"]
            command += ["-evalue", str(evalMax)]
            command += ["-out", inputFile + ".out"]
            command += ["-outfmt", "6"]
            command += ["-max_target_seqs", str(max_target_seqs)]
            blastpProcesses.append(subprocess.Popen(command))            
        # Wait until all the processes have stopped
        for process in blastpProcesses :
            process.wait()
        # Pool the results
        out = mergeFiles([x + ".out" for x in inputFiles], out = out)
        # Erase input and output files
        [os.remove(x) for x in inputFiles]
        [os.remove(x + ".out") for x in inputFiles]
    else :
        raise Exception("jobs and cores cannot be both > 1")
    return out


### ** shuffleFastaFile(fastaFile, outFile)

def shuffleFastaFile(fastaFile, outFile) :
    """Shuffle the sequence order in a fasta file (not mixing the sequences 
    names and sequences, keeping the data as is but just shuffling the order 
    in which sequences are presented in the file). Line wrapping is not 
    preserved (biopython does the line wrapping in the output file).

    All the sequences are loaded into memory at the same time, so there should 
    be enough memory to hold them all.

    Args:
        fastaFile (str): Input fasta file
        outFile (str): Output fasta file

    """
    seqs = SeqIO.parse(fastaFile, "fasta")
    seqList = []
    for seq in seqs :
        seqList.append(seq)
    outOrder = range(len(seqList))
    random.shuffle(outOrder)
    with open(outFile, "w") as fo :
        for i in outOrder :
            SeqIO.write(seqList[i], fo, "fasta")
    
### ** splitFastaFile(fastaFile, n)

def splitFastaFile(fastaFile, n) :
    """Split a large fasta file into smaller subsets

    Args:
        fastaFile (str): Input fasta file
        n (int): Number of output files

    Return:
        list of str: Names of the output files

    """
    # Get the number of sequences (total)
    command = "grep -c \">\" " + fastaFile
    p = subprocess.Popen(command, stdout = subprocess.PIPE, shell = True)
    p.wait()
    nSeqs = int(p.communicate()[0])
    # Get the number of sequences per subset
    nSubsets = [math.floor(nSeqs / n)] * (n - 1)
    nSubsets.append(nSeqs - sum(nSubsets))
    # Go through the initial file and split it
    subsetFiles = [fastaFile + "." + randomTag(32) for x in range(n)]
    assert len(nSubsets) == len(subsetFiles)
    with open(fastaFile, "r") as fi :
        l = fi.next()
        for (subsetFile, nSeq) in zip(subsetFiles, nSubsets) :
            try :
                with open(subsetFile, "w") as fo :
                    nWritten = 0
                    while (nWritten < nSeq) :
                        fo.write(l)
                        l = fi.next()
                        if l.startswith(">") :
                            nWritten += 1
            except StopIteration :
                pass
    # Return
    return subsetFiles

### ** mergeFiles(files, out)

def mergeFiles(files, out) :
    """Merge several text files into a single file

    Args:
        files (list of str): List of file names
        out (str): Name of the output file

    Returns:
        str: Name of the output file

    """
    command = "cat "
    command += " ".join(files)
    command += " > " + out
    os.system(command)
    return out

### ** removeBlastpDb(db)

def removeBlastpDb(db) :
    """Delete the files of a blastp database

    Args:
        db (str): Name of the blastp database to delete (`db` element in the
          output from :func:`makeBlastDb`)

    Returns:
        None

    """
    suffixes = ["phr", "pin", "psq"]
    for s in suffixes :
        os.remove(db + "." + s)
    
### ** prepareBlastpABC(results, out = None)

def prepareBlastpABC(results, out = None) :
    """Convert a result file to ABC format

    Args:
        results (str): Name of the blastp result file
        out (str): Name of the output file. If None, append the '.ABC' suffix 
          to the result file

    Returns:
        Name of the output file

    """
    command = "cut -f 1,2,11 "
    command += results
    if out is None :
        out = results + ".ABC"
    command += " > " + out
    os.system(command)
    return out

### ** ungapFasta(inputFile, out)

def ungapFasta(inputFile, out) :
    """Ungap a fasta file (remove all "-" characters except in sequence names)
    Note: input and output files must be different (they will be opened 
    simultaneously).

    Args:
        inputFile (str): Path to the input fasta file
        out (str): Path to the output file
    """
    assert inputFile != out
    with open(inputFile, "r") as fi :
        with open(out, "w") as fo :
            for l in fi :
                if l.strip().startswith(">") :
                    fo.write(l)
                else :
                    fo.write(l.replace("-", ""))

### ** parseBlastpOutput

def parseBlastpOutput(outputFile) :
    """Load the content of a blastp output file. The format is assumed to be
    the one obtained with -outfmt 6

    Args:
        outputFile (str): Path to the blaspt output file to parse

    Returns:
        A dictionary with qseqid as keys and lists of blastp results as values

    """
    out = dict()
    with open(outputFile, "r") as fi :
        for line in fi :
            line = line.rstrip()
            if line != "" :
                line = line.split("\t")
                result = dict(zip(OUTFMT_HEADER, line))
                out[result["qseqid"]] = out.get(result["qseqid"], [])
                out[result["qseqid"]].append(result)
    return out

### ** getBestMatch

def getBestMatch(blastpOutput) :
    """Get the best match for each qseqid from a blastp result. For each 
    qseqid, a list of the sseqid with the best bitscore is returned (there
    can be several sseqid in the list if there are ties).
    

    Args:
        blastpOutput (dict): Output from parseBlastpOutput()

    Returns:
        Dictionary mapping each qseqid to its list of best sseqid matches

    """
    out = dict()
    for (qseqid, result) in blastpOutput.items() :
        sseqids = [(x["sseqid"], float(x["bitscore"])) for x in result]
        sseqids.sort(key = lambda x: x[1], reverse = True)
        bestBitscore = sseqids[0][1]
        sseqids = [x[0] for x in sseqids if x[1] == bestBitscore]
        out[qseqid] = sseqids
    return out

### ** reciprocalBestMatch

def reciprocalBestMatch(bestMatches1, bestMatches2) :
    """Find the reciprocal best matches from two blastp datasets

    Args:
        bestMatches1 (dict): Output from getBestMatch()
        bestMatches2 (dict): Output from getBestMatch()

    Returns:
        A list of best pairs

    """
    out = set()
    for (qseqid, matches) in bestMatches1.items() :
        for match in matches :
            reciprocalMatches = bestMatches2.get(match, [])
            if qseqid in reciprocalMatches :
                out.add(frozenset([qseqid, match]))
    for (qseqid, matches) in bestMatches2.items() :
        for match in matches :
            reciprocalMatches = bestMatches1.get(match, [])
            if qseqid in reciprocalMatches :
                out.add(frozenset([qseqid, match]))
    out = list(out)
    cleanOut = [list(x) for x in out]
    cleanOut = [x if x[0] in bestMatches1.keys() else [x[1], x[0]] for x in cleanOut]
    return cleanOut

