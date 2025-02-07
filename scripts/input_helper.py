import argparse
import pathlib as pl
from pathlib import Path
import json

def get_args():
    parser = argparse.ArgumentParser("helper")

    subcom = parser.add_subparsers(dest='subcommand')
    parser_multi = subcom.add_parser('multiple')
    parser_multi.add_argument('-i', required=True, help='input dir')
    parser_multi.add_argument('-t', required=True, help='template json')
    parser_multi.add_argument('-o', default='/dev/stdout', help='output json or dir for separated sample json')
    parser_multi.add_argument('-s', action='store_true', help='separate sample to different json?')
    parser_multi.add_argument('-g', action='store_true', help='is gz?')
    parser_multi.add_argument('-l', action='store_true', help='input is a file list')

    parser_single = subcom.add_parser('single')
    parser_single.add_argument('-i', required=True, help='input vcf')
    parser_single.add_argument('-t', required=True, help='template json')
    parser_single.add_argument('-o', default='/dev/stdout', help='output json')

    args = parser.parse_args()
    return args

def load_template(template_path):
    template_path = Path(template_path)
    assert template_path.exists(), f"Template file {template_path} not exists"
    template = json.load(open(template_path, 'r'))
    return template

def abs_path_str(path: Path):
    return str(Path(path).resolve())

def multi(args):
    input_dir = Path(args.i)
    input_js = load_template(args.t)
    patient_base = input_js['patient'][0]
    if not args.s:
        patients = []
        file_generator = (
            input_dir.glob(f"*.vcf{'.gz' if args.g else ''}") 
                if not args.l else 
            open(input_dir, 'r').read().strip().split()
        )
        for vcf in file_generator:
            patient = dict(patient_base)
            patient['patient_vcf'] = abs_path_str(vcf)
            patients.append(patient)

        input_js['patient'] = patients
        json.dump(input_js, open(args.o, 'w'), indent='\t')
    else:
        args.o = Path(args.o)
        assert args.o != '/dev/stdout' and args.o.is_dir(), "Need output dir in separated sample json mode"
        for vcf in input_dir.glob(f"*.vcf{'.gz' if args.g else ''}"):
            patient = dict(patient_base)
            patient['patient_vcf'] = abs_path_str(vcf)
            
            sample_name = vcf.stem
            copy_input_js = dict(input_js)
            copy_input_js['patient'] = [patient]
            json.dump(copy_input_js, open(args.o / f"input_{sample_name}.json", 'w'), indent='\t')

def single(args):
    input_js = load_template(args.t)
    patient_base = input_js['patient'][0]
    patient = dict(patient_base)
    patient['patient_vcf'] = abs_path_str(args.i)
    input_js['patient'] = [patient]
    json.dump(input_js, open(args.o, 'w'), indent='\t')

if __name__ == '__main__':
    args = get_args()
    if args.subcommand == 'multiple':
        multi(args)
    elif args.subcommand == 'single':
        single(args)
    else:
        raise NotImplementedError