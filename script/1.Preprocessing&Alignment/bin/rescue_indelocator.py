#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import gzip
import argparse
import logging
import datetime


import vcf

LOG = logging.getLogger(__name__)

__author__ = ("hu zhe",)
__version__ = "1.0.0"
__email__ = "u201810332@hust.edu.cn"


def rescue_indelocator(input_vcf, min_tumor_dp, min_normal_dp, min_tumor_af,  max_normal_af, tumor_ID):
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    vcf_W1 = vcf.Writer(open(tumor_ID+".indelocator.somatic_only.vcf",'w'), vcf_reader)
    for record in vcf_reader:
        if record.INFO.has_key('SOMATIC'):
            vcf_W1.write_record(record)
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    vcf_W2 = vcf.Writer(open(tumor_ID+".indelocator.rescue_only.vcf",'w'),vcf_reader)
    for record in vcf_reader:
        if ((not(record.INFO.has_key('SOMATIC'))) and(float(record.INFO['T_AC'][1]) / float(record.INFO['T_DP']) >= min_tumor_af) and (float(record.INFO['N_AC'][1]) / float(record.INFO['N_DP']) <= max_normal_af ) and (record.INFO['N_DP'] >= min_normal_dp) and (record.INFO['T_DP'] >= min_tumor_dp) ):
            vcf_W2.write_record(record)
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    vcf_W3 = vcf.Writer(open(tumor_ID+".indelocator.total.vcf",'w'),vcf_reader)
    for record in vcf_reader:
        if (record.INFO.has_key('SOMATIC')) or ((float(record.INFO['T_AC'][1]) / float(record.INFO['T_DP']) >= min_tumor_af) and (float(record.INFO['N_AC'][1]) / float(record.INFO['N_DP']) <= max_normal_af) and (record.INFO['N_DP'] >= min_normal_dp) and (record.INFO['T_DP'] >= min_tumor_dp) ) :
            vcf_W3.write_record(record)
    return 0


def add_args(parser):

    parser.add_argument("input_vcf", help="vcf")
    parser.add_argument("--min_tumor_dp", type=int, default=50, help="minimum depth of position in tumor (50).")
    parser.add_argument("--min_normal_dp", type=int, default=50, help="minimum depth of position in normal (50).")
    parser.add_argument("--min_tumor_af", type=float, default=0.2, help="minimum variant allele fraction in tumor (0.20).")
    parser.add_argument("--max_normal_af", type=float, default=0.05, help="max variant allele fraction in normal (0.05).")
    parser.add_argument("--tumor_ID", required=True, help="prefix of output vcf file")

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
rescue_indelocator.py - a tool to rescue indel  from indelocator vcf according to 2016 NG paper

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_args(parser)
    args = parser.parse_args()

    rescue_indelocator(
        args.input_vcf, args.min_tumor_dp, args.min_normal_dp,
        args.min_tumor_af,  args.max_normal_af, args.tumor_ID
    )


if __name__ == "__main__":
    main()
