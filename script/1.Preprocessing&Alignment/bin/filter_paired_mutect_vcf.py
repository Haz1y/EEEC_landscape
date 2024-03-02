#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import gzip
import argparse
import logging
import datetime


from pysam import VariantFile

LOG = logging.getLogger(__name__)

__author__ = ("fan junpeng",)
__version__ = "1.0.0"
__email__ = "jpfan@foxmail.com"


def filter_mutect2_vcf(vcf, filter_name, min_tumor_dp, min_normal_dp, min_tumor_ad, min_tumor_af,  max_normal_af, out):

    vcf_in  = VariantFile(vcf, "r")

    if out.endswith("gz"):
        vcf_out = gzip.open(out, "wt")
    else:
        vcf_out = open(out, "w")

    version = __version__
    now = datetime.datetime.now()
    date_time = now.__format__("%B %d, %Y %I:%M:%S %p CST")
    vcf_in.header.add_line(
        f'##CommandLine=<ID=fiter_paired_mutect_vcf,CommandLine="fiter_paired_mutect_vcf.py --filter {filter_name} --min_tumor_dp {min_tumor_dp} \
        --min_normal_dp {min_normal_dp} --min_tumor_ad {min_tumor_ad} --min_tumor_af {min_tumor_af} --max_normal_af {max_normal_af}",Version="{version}",Date="{date_time}"'
        )
    header = str(vcf_in.header).split("\n")
    _h = header[-2].split("\t")
    header = header[:-2] + ["\t".join(_h[:-2] + [_h[-1]])]
    vcf_out.write("\n".join(header)+"\n")
    # get tumor and normal
    tumor_name = ""
    normal_name = ""
    for x in vcf_in.header.records:
        if x.key == "tumor_sample":
            tumor_name = x.value

        if x.key == "normal_sample":
            normal_name = x.value
    if not (tumor_name and normal_name):
        LOG.error("tumor %r and normal %r is not paired!")
        sys.exit(1)

    for rec in vcf_in:

        tumor = rec.samples[tumor_name]
        normal = rec.samples[normal_name]

        if tumor["DP"] >= min_tumor_dp and normal["DP"] >= min_normal_dp and filter_name in rec.filter:
            pass_alts = []
            for n, alt in enumerate(rec.alts):
                tumor_pass = normal_pass = False
                if tumor["AD"][n+1] >= min_tumor_ad and tumor["AF"][n] >= min_tumor_af:
                    tumor_pass =True
                if normal["AF"][n] < max_normal_af:
                    normal_pass = True

                if tumor_pass and normal_pass:
                    pass_alts.append(n)
            if pass_alts:
                if rec.id is None:
                    _id = "."
                else:
                    _id = rec.id
                _filter = ";".join(rec.filter)
                _alt = ",".join([rec.alts[i] for i in pass_alts])

                _info = ""
                for k, v in rec.info.items():
                    if isinstance(v, tuple):
                        _v = []
                        for i in v:
                            if isinstance(i, float):
                                _v.append("%g" % i)
                            else:
                                _v.append(str(i))
                        _v = ",".join(_v)
                    else:
                        _v = v
                    _info += "%s=%s;" %(k, _v)
                _format = ":".join(rec.format.keys())

                _sample = []
                for sample in rec.samples:
                    if sample == normal_name:
                        continue
                    _s = ""
                    for k in rec.format.keys():
                        _v = []
                        if k == "AD":
                            _v = [str(rec.samples[sample][k][i+1]) for i in ([-1] + pass_alts)]
                        elif k == "AF":
                            _v = ["%g" % rec.samples[sample][k][i] for i in pass_alts]
                        elif k == "DP":
                            _v = [str(rec.samples[sample][k])]
                        elif k == "GT":
                            _v = ["|".join(map(str,rec.samples[sample][k]))]
                        elif k == "PGT":
                            _v = ["".join(map(str,rec.samples[sample][k]))]
                        elif k == "PS":
                            _v = [str(rec.samples[sample][k])]
                        elif k == "PID":
                            _v = ["".join(map(str,rec.samples[sample][k]))]
                        else:
                            _v = map(str, rec.samples[sample][k])
                        _s += ",".join(_v)+":"
                    _sample.append(_s[:-1])
                _sample = "\t".join(_sample)
                vcf_out.write(f"{rec.chrom}\t{rec.pos}\t{_id}\t{rec.ref}\t{_alt}\t.\t{_filter}\t{_info}\t{_format}\t{_sample}\n")

    vcf_out.close()

    return 0


def add_args(parser):

    parser.add_argument("vcf", help="vcf")
    parser.add_argument("--filter", default="PASS", help="allowed filter (PASS).")
    parser.add_argument("--min_tumor_dp", type=int, default=8, help="minimum depth of position in tumor (8).")
    parser.add_argument("--min_normal_dp", type=int, default=8, help="minimum depth of position in normal (8).")
    parser.add_argument("--min_tumor_ad", type=int, default=5, help="minimum variant allele depth in tumor (5).")
    parser.add_argument("--min_tumor_af", type=float, default=0.05, help="minimum variant allele fraction in tumor (0.05).")
    parser.add_argument("--max_normal_af", type=float, default=0.02, help="max variant allele fraction in normal (0.02).")
    parser.add_argument("-o", "--output", required=True, help="output vcf file")

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

    filter_mutect2_vcf(
        args.vcf, args.filter, args.min_tumor_dp, args.min_normal_dp,
        args.min_tumor_ad, args.min_tumor_af,  args.max_normal_af, args.output
    )


if __name__ == "__main__":
    main()
