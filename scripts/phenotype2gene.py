import pandas as pd
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Phenotype to gene')
    parser.add_argument('-i', '--input', type=str, required=True,
        help='Input file, shoud at least contain a column named "phenotype", which contains HPO terms separated by ", "')
    parser.add_argument('-o', '--output', type=str, required=True,
        help='Output file, will be input file but added a column named "genes" with the genes associated with the phenotypes')
    parser.add_argument('-p', '--hpo', type=str, required=True,
        help='HPO file, should be the file `phenotype_to_genes.txt` downloaded from https://hpo.jax.org/app/data/annotations')
    args = parser.parse_args()
    return args

def main(args):
    # header is: hpo_id	hpo_name	ncbi_gene_id	gene_symbol	disease_id
    hpo = pd.read_csv(args.hpo, sep='\t', engine='pyarrow', dtype_backend='pyarrow')

    cases = pd.read_csv(args.input)

    case_genes = []
    for i, row in cases.iterrows():
        phenotypes = row['phenotype'].split(', ')
        genes = []
        for p in phenotypes:
            gene = hpo[hpo['hpo_name'] == p]['gene_symbol'].values
            if len(gene) > 0:
                genes.extend(gene)
            else:
                print('No gene found for {}'.format(p))
        case_genes.append(genes)

    cases['genes'] = [';'.join(set(sorted(g))) for g in case_genes]
    cases.to_csv(args.output, index=False)

if __name__ == '__main__':
    args = get_args()
    main(args)