""" splits de novo clustering analysis across multiple jobs on a LSF cluster.
"""

import os
import subprocess
import time
import random
import argparse

def get_options():
    """ get the command line options
    """
    
    parser = argparse.ArgumentParser(description="Script to batch process de"
        "novo clustering.")
    parser.add_argument("--script", required=True, help="Path to denovonear"
        "clustering script")
    parser.add_argument("--in", dest="input", required=True, help="Path to"
        "file listing known mutations in genes. See example file in data folder"
        "for format.")
    parser.add_argument("--rates", required=True, help="Path to rates file.")
    parser.add_argument("--out", required=True, help="Path to output file.")
    
    args = parser.parse_args()
    
    return args

def split_denovos(denovo_path):
    """ split de novos from an input file into files containing 100 de novos
    """
    
    # open the de novos, drop the header line, then sort them (which will be by
    # HGNC symbol, as the first element of each line)
    f = open(denovo_path, "r")
    lines = f.readlines()
    header = lines.pop(0)
    
    lines = sorted(lines)
    
    prior = ""
    count = 0
    paths = []
    for line in lines:
        count += 1
        tmp = line.split("\t")
        gene = tmp[0]
        
        # once there are more than 100 de novos in the file, and we have cleared
        # any in progress genes, then get ready to start a new file
        if count > 100 and gene != prior:
            count = 1
        
        # open a file handle to write the de novos to
        if count == 1:
            path = "{0}.{1}.txt".format(denovo_path, (len(paths) + 1))
            output = open(path, "w")
            output.write(header)
            paths.append(path)
        
        prior = gene
        output.write(line)
    
    return paths

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

def batch_process_clustering(script, paths, de_novo_path, rates_path, output_path):
    """ sets up a lsf job array
    """
    
    # set up run parameters
    job_name = "denovonear"
    job_id = "{0}[1-{1}]%10".format(job_name, len(paths))
    
    command = ["python", script, \
        "--in", de_novo_path + ".\$LSB_JOBINDEX\.txt", \
        "--out", de_novo_path + ".\$LSB_JOBINDEX\.output.txt", \
        "--rates", rates_path]
    submit_bsub_job(command, job_id, memory=3000)
    time.sleep(2)
    
    # merge the array output after the array finishes
    merge_id = "merge1_" + job_name
    command = ["head", "-n", "1", de_novo_path + ".1.output.txt", ">", output_path, \
        "; tail", "-q", "-n", "+2", de_novo_path + ".*.output.txt", ">>", output_path]
    submit_bsub_job(command, merge_id, dependent_id=job_id)
    time.sleep(2)
    
    # submit a cleanup job to the cluster
    cleanup_id = "cleanup"
    command = ["rm", de_novo_path + ".*.txt"]
    submit_bsub_job(command, cleanup_id, dependent_id=merge_id)
    
def main():
    args = get_options()
    paths = split_denovos(args.input)
    batch_process_clustering(args.script, paths, args.input, args.rates, args.output)

if __name__ == '__main__':
    main()
