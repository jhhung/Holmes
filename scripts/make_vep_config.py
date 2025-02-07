import json
from pathlib import Path
import argparse

def get_args():
    parser = argparse.ArgumentParser("Making Holmes database config json file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output config file path")
    parser.add_argument("--vep_repo_dir", type=str, default="", help="VEP repository directory")
    parser.add_argument("--vep_data_dir", type=str, default="", help="VEP data directory")
    parser.add_argument("--assembly", type=str, default="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", help="Genome assembly fa gz file")
    parser.add_argument("--MaxEntScan_file", type=str, default="fordownload",
                        help="MaxEntScan model & script dir, can be downloaded and extracted from http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz")
    parser.add_argument("--REVEL_file", type=str, default="new_tabbed_revel_grch38.tsv.gz",
                        help="REVEL score file, can be downloaded and made by instruction from https://github.com/Ensembl/VEP_plugins/blob/release/110/REVEL.pm")
    parser.add_argument("--CADD_file", nargs='+', 
                        default=["GRCh38_whole_genome_SNVs.tsv.gz", "GRCh38_gnomad.genomes.r3.0.indel.tsv.gz"],
                        help="CADD score file, can be more than one file")
    args = parser.parse_args()
    return args

def main(args):
    config = {
        "options": [
            "--vcf",
            "--offline",
            "--force_overwrite", 
            "--variant_class",
            "--symbol",
            "--sift b",
            "--polyphen b",
            "--numbers",
            "--domains",
            "--canonical",
            "--hgvs",
            "--hgvsg",
            "--mane",
            "--mane_select",
            "--merged",
            "--no_stats",
            "--quiet"
        ],
        "plugins": {
            "MaxEntScan": [
                args.MaxEntScan_file
            ],
            "pLI": [],
            "REVEL": [
                args.REVEL_file
            ],
            "CADD": args.CADD_file,
            "NMD": []
        },
        "assembly": args.assembly,
        "vep_repo_dir": args.vep_repo_dir,
        "vep_data_dir": args.vep_data_dir
    }

    with open(args.output, "w") as f:
        json.dump(config, f, indent=4)

if __name__ == "__main__":
    args = get_args()
    main(args)
