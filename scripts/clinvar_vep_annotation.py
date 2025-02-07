#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path as P
def getargs():
    parser = argparse.ArgumentParser(description="Annotate ClinVar VCF with VEP")
    parser.add_argument("-i", "--input", type=str, help="Input VCF file", required=True)
    parser.add_argument("-v", "--vep", type=P, help="VEP binary", required=True)
    parser.add_argument("-f", "--fasta", type=str, help="Fasta ref file", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output vep annotation file", required=True)
    parser.add_argument("-t", "--thread", type=int, help="fork of VEP", default=1)
    parser.add_argument("--grch37", action='store_true', help="use GRCh37 assembly")
    args = parser.parse_args()
    args.thread = args.thread if args.thread > 0 else 1
    return args

def main(args):
    vep_args = [
        args.vep.expanduser(),
        "--vcf",
        "--offline",
        "--force_overwrite",
        "--assembly", "GRCh38" if not args.grch37 else "GRCh37",
        "--symbol",
        "--numbers",
        "--domains",
        "--hgvs",
        "--hgvsg",
        "--merged",
        "--fork", str(args.thread),
        "--no_stats",
        "--quiet",
        "-i", args.input,
        "-o", args.output,
        "--fasta", args.fasta,
    ]
    ret = subprocess.Popen(args=vep_args)
    returncode = ret.wait()
    if returncode:
        raise Exception(f"VEP annotation failed, return code: {returncode}")

if __name__ == "__main__":
    args = getargs()
    main(args)