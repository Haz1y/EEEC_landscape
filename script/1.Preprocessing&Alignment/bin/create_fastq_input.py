#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import logging
import pandas as pd

LOG = logging.getLogger(__name__)

__author__ = ("fan junpen",)
__version__ = "1.0.0"
__email__ = "jpfan@foxmail.com"


def merge_input(file):

    r = []
    run = pd.read_excel(file, sheet_name='run')
    sample = pd.read_excel(file, sheet_name='sample')
    r = run.merge(sample, how='left', on='sample').loc[:,\
        ["run","sample","sample_type","patient","platform","lane",\
        "forward_file_name","reverse_file_name"]].rename(columns={"run": "id"})
    r["sample_type"] = r["sample_type"].str.replace("癌旁","N").replace("原发","T").replace("转移","T")

    return r

def add_args(parser):

    parser.add_argument("input_xlsx", help="run and sample xlsx.")
    parser.add_argument("output_csv", help="output csv for fastq channel")

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
create_fastq_input.py - a tool to prepare input of fastq channel from xlsx

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_args(parser)
    args = parser.parse_args()

    fastq_csv = merge_input(
                args.input_xlsx
    )
    fastq_csv.to_csv(args.output_csv,index=False)





if __name__ == "__main__":
    main()
