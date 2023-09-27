#!/usr/bin/env python3
"""Summarize read counts from kneaddata."""
__author__ = "SX"
__date__ = "20230921"
__version__ = 0.1


from sys import argv, exit
from pathlib import Path
import argparse
import re
import os

import pandas as pd

def parse_args():
    desc = f"{__doc__} Copyright (c) {__author__} {__date__}. Version v{__version__}"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--kneaddata", metavar="sample_stat.log", nargs="+",
            help="kneaddata standard output log file.")
    parser.add_argument("-o", "--output-table", metavar="TSV",
            default="read_processing_summary.txt",
            help="Filename of output table in tsv format [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()



def parse_kneaddata_log(logfiles):
    loglist = []
    for logfile in logfiles:
        logdict = {}
        with open(logfile) as f:
            logdict['sample_name'] = Path(logfile).stem.split("_")[0]
            for line in f:
                if re.search("Initial number of reads",line):
                    temp_raw = int(float(line.split(':')[-1]))
                    if logdict.get('raw_read'):
                        if logdict.get('raw_read') != temp_raw:
                            raise ValueError('Raw fastq had different number of reads.')
                    else:
                        logdict['raw_read'] = temp_raw
                if re.search("Total reads after trimming.*trimmed.1.fastq",line):
                    logdict['after_trimmomatic'] = int(float(line.split(':')[-1]))
                if re.search("Total reads after merging results from multiple databases.*paired_1.fastq",line):
                    logdict['after_bowtie2_host_removal'] = int(float(line.split(':')[-1]))
                
            loglist.append(logdict)
    return loglist
        

if __name__ == "__main__":
    
    args = parse_args()

    knead_sumlist = [os.path.join(args.kneaddata[0], file) for file in os.listdir(args.kneaddata[0]) if re.search('_stat.log',file)]

    data_kneaddata = parse_kneaddata_log(knead_sumlist)

    df = pd.DataFrame(data_kneaddata).set_index("sample_name")

    df.sort_index(inplace = True)

    df.to_csv(args.output_table, sep="\t")



