import json
from pathlib import Path
import argparse


def get_args():
    parser = argparse.ArgumentParser("Making VEP config json file")
    parser.add_argument(
        "-b", "--base_dir",
        default="",
        help='Base directory for database files (default: "", this will be replaced with ~/holmes_database/ in holmes runtime)',
    )
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument("--k_genome", default="", help="Path to K genome database file")
    parser.add_argument("--clinvar", default="", help="Path to ClinVar database file")
    parser.add_argument("--dvd", default="", help="Path to DVD database file")
    parser.add_argument("--coverage", default="", help="Path to coverage database file")
    parser.add_argument("--gnom", default="", help="Path to gnomAD database file")
    parser.add_argument(
        "--gene_info", default="", help="Path to gene info database file"
    )
    parser.add_argument("--fasta", default="", help="Path to FASTA database file")
    parser.add_argument("--gtf", default="", help="Path to GTF database file")
    parser.add_argument("--uniprot", default="", help="Path to UniProt database file")
    return parser.parse_args()


def main(args):
    if args.base_dir != "":
        norm_base_dir = str(Path(args.base_dir).expanduser().resolve())
    else:
        norm_base_dir = ""
    config = {
        "base_dir": norm_base_dir,
        "db": {
            "k_genome": args.k_genome,
            "clinvar": args.clinvar,
            "dvd": args.dvd,
            "coverage": args.coverage,
            "gnom": args.gnom,
            "gene_info": args.gene_info,
            "fasta": args.fasta,
            "gtf": args.gtf,
            "uniprot": args.uniprot,
        },
    }

    with open(args.output, "w") as f:
        json.dump(config, f, indent=4)


if __name__ == "__main__":
    args = get_args()
    main(args)
