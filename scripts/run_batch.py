""" splits de novo clustering analysis across multiple jobs on a LSF cluster.
"""

import os
import subprocess
import time
import random
import argparse
import tempfile

def get_options():
    """ get the command line options
    """
    
    parser = argparse.ArgumentParser(description="Script to batch process de"
        "novo clustering.")
    parser.add_argument("--in", dest="input", required=True, help="Path to"
        "file listing known mutations in genes. See example file in data folder"
        "for format.")
    parser.add_argument("--temp-dir", required=True, help="path to hold intermediate files")
    parser.add_argument("--out", required=True, help="Path to output file.")
    
    args = parser.parse_args()
    
    return args

def count_missense_per_gene(lines):
    """ count the number of missense variants in each gene.
    """
    
    counts = {}
    for x in lines:
        x = x.split("\t")
        gene = x[0]
        consequence = x[3]
        
        if gene not in counts:
            counts[gene] = 0
        
        if consequence != "missense_variant":
            continue
        
        counts[gene] += 1
    
    return counts

def split_denovos(denovo_path, temp_dir):
    """ split de novos from an input file into files, one for each gene
    """
    
    # open the de novos, drop the header line, then sort them (which will be by
    # HGNC symbol, as the first element of each line)
    with open(denovo_path, "r") as handle:
        lines = handle.readlines()
        header = lines.pop(0)
    
    basename = os.path.basename(denovo_path)
    
    # only include genes with 2+ missense SNVs
    counts = count_missense_per_gene(lines)
    counts = dict((k, v) for k, v in counts.items() if v > 1 )
    
    genes = set([])
    for line in sorted(lines):
        gene = line.split("\t")[0]
        
        # open a new output file whenever we hit a different gene
        if gene not in genes and gene in counts:
            genes.add(gene)
            path = os.path.join(temp_dir, "tmp.{}.txt".format(len(genes)))
            output = open(path, "w")
            output.write(header)
        
        if gene in counts:
            output.write(line)
    
    return len(genes)

def is_number(string):
    """ check whether a string can be converted to a number
    
    Args:
        string: value as a string, could be a number
        
    Returns:
        True/False for whether the value can be converted to a number
    """
    
    try:
        number = float(string)
    except ValueError:
        return False
    
    return True

def get_random_string():
    """ make a random string, which we can use for bsub job IDs, so that
    different jobs do not have the same job IDs.
    """
    
    # set up a random string to associate with the run
    hash_string = "%8x" % random.getrandbits(32)
    hash_string = hash_string.strip()
    
    # done't allow the random strings to be equivalent to a number, since
    # the LSF cluster interprets those differently from letter-containing
    # strings
    while is_number(hash_string):
        hash_string = "%8x" % random.getrandbits(32)
        hash_string = hash_string.strip()
    
    return hash_string

def submit_bsub_job(command, job_id=None, dependent_id=None, memory=None, requeue_code=None, logfile=None):
    """ construct a bsub job submission command
    
    Args:
        command: list of strings that forma unix command
        job_id: string for job ID for submission
        dependent_id: job ID, or list of job IDs which the current command needs
            to have finished before the current command will start. Note that
            the list can be empty, in which case there are no dependencies.
        memory: minimum memory requirements (in megabytes)
    
    Returns:
        nothing
    """
    
    if job_id is None:
        job_id = get_random_string()
    
    job = "-J \"{0}\"".format(job_id)
    
    mem = ""
    if memory is not None:
        mem = "-R 'select[mem>{0}] rusage[mem={0}]' -M {0}".format(memory)
    
    requeue = ""
    if requeue_code is not None:
        requeue = "-Q 'EXCLUDE({0})'".format(requeue_code)
    
    dependent = ""
    if dependent_id is not None:
        if type(dependent_id) == list:
            dependent_id = " && ".join(dependent_id)
        dependent = "-w '{0}'".format(dependent_id)
    
    log = "bjob_output.txt"
    if logfile is not None:
        log = logfile
    
    preamble = ["bsub", job, dependent, requeue, "-q", "normal", "-o", log, mem]
    command = ["bash", "-c", "\""] + command + ["\""]
    
    command = " ".join(preamble + command)
    subprocess.call(command, shell=True)

def batch_process(de_novo_path, temp_dir, output_path):
    """ sets up a lsf job array
    """
    
    temp_dir = tempfile.mkdtemp(dir=temp_dir)
    count = split_denovos(de_novo_path, temp_dir)
    
    # set up run parameters
    job_name = "denovonear"
    job_id = "{0}[1-{1}]%20".format(job_name, count)
    
    basename = os.path.basename(de_novo_path)
    infile = os.path.join(temp_dir, "tmp.\$LSB_JOBINDEX\.txt")
    outfile = os.path.join(temp_dir, "tmp.\$LSB_JOBINDEX\.output")
    
    command = ["denovonear", "cluster",
        "--in", infile,
        "--out", outfile]
    submit_bsub_job(command, job_id, memory=3500, requeue_code=134, logfile="clustering.bjob")
    time.sleep(2)
    
    # merge the array output after the array finishes
    merge_id = "{}_merge".format(job_name)
    command = ["head", "-n", "1", os.path.join(temp_dir, "tmp.1.output"), ">", output_path, \
        "; tail", "-q", "-n", "+2", os.path.join(temp_dir, "tmp.*.output"), "|", "sort", ">>", output_path]
    submit_bsub_job(command, merge_id, memory=100, dependent_id=job_id, logfile="clustering.bjob")
    time.sleep(2)
    
    # submit a cleanup job to the cluster
    submit_bsub_job(["rm", "-r", temp_dir], job_id="{}_cleanup".format(job_name), \
        memory=100, dependent_id=merge_id, logfile="clustering.bjob")
    
def main():
    args = get_options()
    
    batch_process(args.input, args.temp_dir, args.out)

if __name__ == '__main__':
    main()
