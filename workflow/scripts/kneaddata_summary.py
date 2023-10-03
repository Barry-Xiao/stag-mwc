#!/usr/bin/env python3
"""Summarize read counts from kneaddata."""


from sys import argv, exit
from pathlib import Path
import argparse
import re
import os

import pandas as pd
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt

def parse_args():
    desc =  "Summarize read counts from kneaddata. Adopted from preprocessing_summary.py by CTMR, Fredrik Boulund"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--kneaddata", metavar="kneaddata/sample.log", nargs="+",
            help="kneaddata standard output log file.")
    parser.add_argument("--kraken2", metavar="host_removal/sample.kraken2.log", nargs="+",
            help="Kraken2 log output.")
    parser.add_argument("-o", "--output-table", metavar="TSV",
            default="read_processing_summary.txt",
            help="Filename of output table in tsv format [%(default)s].")
    parser.add_argument("-p", "--output-plot", metavar="PDF",
            default="",
            help="Filename of output table in PDF format [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()



def parse_kneaddata_log(logfiles):
    for logfile in logfiles:
        logdict = {}
        with open(logfile) as f:
            logdict['Sample'] = Path(logfile).stem.split(".")[0]
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
                
                yield logdict

def parse_kraken2_logs(logfiles):
    for logfile in logfiles:
        with open(logfile) as f:
            sample_name = Path(logfile).stem.split(".")[0]
            for line in f:
                if " unclassified" in line:
                    yield {
                        "Sample": sample_name,
                        "after_kraken2_host_removal": int(line.strip().split()[0]),
                    }


if __name__ == "__main__":
    
    args = parse_args()

    dfs = {
        "kneaddata": pd.DataFrame(),
        "kraken2": pd.DataFrame(),
    }

    data_kneaddata = list(parse_kneaddata_log(args.kneaddata))
    dfs["kneaddata"] = pd.DataFrame(data_kneaddata).set_index("Sample")
    if args.kraken2:
        data_kraken2 = list(parse_kraken2_logs(args.kraken2))
        dfs["kraken2"] = pd.DataFrame(data_kraken2).set_index("Sample")

    df = pd.concat(dfs.values(), axis="columns")

    column_order = [
        "raw_read",
        "after_trimmomatic",
        "after_bowtie2_host_removal",
        "after_kraken2_host_removal",
    ]

    final_columns = [c for c in column_order if c in df.columns]
    df = df[final_columns]

    df.to_csv(args.output_table, sep="\t")

    if args.output_plot:
        fig, ax = plt.subplots(figsize=(6, 5))
        df[final_columns[1:]]\
            .transpose()\
            .plot(kind="line", style=".-", ax=ax)
        ax.set_title("Reads passing through QC and host removal")
        ax.set_xlabel("Stage")
        ax.set_ylabel("Reads")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(0, -0.1))
        fig.savefig(args.output_plot, bbox_inches="tight")


