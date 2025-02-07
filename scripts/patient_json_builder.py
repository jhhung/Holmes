import argparse
import json

def parse_args():
    parser = argparse.ArgumentParser(description="Build patient JSON file")
    parser.add_argument("patient_vcf", help="patient VCF file")
    parser.add_argument("--dad_vcf", help="dad VCF file")
    parser.add_argument("--mom_vcf", help="mom VCF file")
    parser.add_argument("--is_sick", action="store_true", help="is test person sick")
    parser.add_argument("--is_denovo", action="store_true", help="is test person denovo")
    parser.add_argument("--sex", action="store_true", help="is test person a male")
    parser.add_argument("--is_dad_sick", action="store_false", help="dad is healthy")
    parser.add_argument("--is_mom_sick", action="store_false", help="mom is healthy")
    parser.add_argument("--healthy_family", help="healthy family list file")
    parser.add_argument("--sick_family", help="sick family list file")
    parser.add_argument("--output", default="/tmp/patient.json", help="output file path")
    return parser.parse_args()

def build_json(args):
    patient = {
        "patient_vcf": args.patient_vcf,
        "dad_vcf": args.dad_vcf,
        "mom_vcf": args.mom_vcf,
        "is_dad_sick": args.is_dad_sick,
        "is_mom_sick": args.is_mom_sick,
        "sex": args.sex,
        "is_sick": args.is_sick,
        "is_denovo": args.is_denovo
    }

    if args.healthy_family:
        with open(args.healthy_family) as f:
            healthy = [line.strip() for line in f]
        patient["healthy_family"] = healthy

    if args.sick_family:
        with open(args.sick_family) as f:
            sick = [line.strip() for line in f]
        patient["sick_family"] = sick

    with open(args.output, "w") as f:
        json.dump(patient, f, indent=4)

if __name__ == "__main__":
    args = parse_args()
    build_json(args)