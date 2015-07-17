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
from Bio import AlignIO

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
