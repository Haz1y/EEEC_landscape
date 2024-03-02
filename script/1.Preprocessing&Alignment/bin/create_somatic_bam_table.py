#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import logging
import pandas as pd

LOG = logging.getLogger(__name__)

__author__ = ("fan junpeng",)
__version__ = "1.0.0"
__email__ = "jpfan@foxmail.com"


def find_paired_samples(file):

    r = []
    df = pd.read_csv(file, header=0, index_col=0)

    for i, row in df.iterrows():
        if row["sample_type"] == "T":
            normal_index = df[(df["patient"] == row["patient"]) & (df["sample_type"] == "N")].index
            if normal_index.empty:
                r.append(
                    [i, row["sample"], row["sample_type"],
                    row["patient"], row["platform"],
                    row["bam"], row["bai"],
                    "", ""]
                    )
            else:
                r.append(
                    [i, row["sample"], row["sample_type"],
                    row["patient"], row["platform"],
                    row["bam"], row["bai"],
                    df.loc[normal_index[0], "bam"], df.loc[normal_index[0], "bai"]]
                    )
        else:
            pass

    return r

def add_args(parser):

    parser.add_argument("input_csv", help="mapped csv from wgs.")
    parser.add_argument("output_csv", help="output csv for somatic mutaion calling.")

    return parser

def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
filter_paired_mutect_vcf.py - a tool to filter mutect2 vcf from tumor-normal pair

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_args(parser)
    args = parser.parse_args()

    pairs = find_paired_samples(
        args.input_csv
    )


    with open(args.output_csv, "w") as fh:
        fh.write("id,sample,sample_type,patient,platform,tumor_bam,tumor_bai,normal_bam,normal_bai\n")
        for i in pairs:
            fh.write(",".join(i)+"\n")


if __name__ == "__main__":
    main()
