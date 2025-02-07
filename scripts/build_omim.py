import pandas as pd
import typer
from pathlib import Path
from collections import defaultdict

PATTERNS = {
    "Autosomal dominant",
    "Autosomal recessive",
    "Pseudoautosomal dominant",
    "Pseudoautosomal recessive",
    "Digenic dominant",
    "Digenic recessive",
    "Isolated cases",
    "Inherited chromosomal imbalance",
    "Mitochondrial",
    "Multifactorial",
    "Somatic mosaicism",
    "Somatic mutation",
    "X-linked",
    "X-linked dominant",
    "X-linked recessive",
    "Y-linked",
}
UNKNOWN = "Unknown"
SEP_INHE = '$'
SEP_PHENO = '|'
SEP_PATTERN = '/'

def search_patterns(pheno: str):
    pats: list[str] = []
    for pat in PATTERNS:
        if pat in pheno:
            pats.append(pat)

    if len(pats) == 0:
        ret = UNKNOWN
    elif len(pats) > 1:
        print(f"WARNING: '{pheno}' has multiple patterns: {pats}")
        ret = SEP_PATTERN.join(pats)
    else:
        ret = pats[0]
    return ret

def get_omim_table(omim: Path):
    from io import StringIO
    omim_lines = omim.read_text().splitlines()
    omim_lines = filter(lambda x: (not x.startswith('#')) or x.startswith('# Chromosome\t'), omim_lines)
    omim_text = '\n'.join(omim_lines)
    df = pd.read_csv(
        StringIO(omim_text),
        sep="\t",
        dtype={
            '# Chromosome': str,
            'Genomic Position Start': int,
            'Genomic Position End': int,
        }
    )
    return df

def main(
    omim: Path
):
    df = get_omim_table(omim)
    print(list(df))
    # cols: '# Chromosome', 'Genomic Position Start', 'Genomic Position End', 
    # 'Cyto Location', 'Computed Cyto Location', 'Mim Number', 'Gene Symbols', 
    # 'Gene Name', 'Approved Symbol', 'Entrez Gene ID', 'Ensembl Gene ID', 
    # 'Comments', 'Phenotypes', 'Mouse Gene Symbol/ID'
    df['Phenotypes'] = df['Phenotypes'].fillna('')
    df['Gene Symbols'] = df['Gene Symbols'].fillna('')
    df['Approved Symbol'] = df['Approved Symbol'].fillna('')

    db = defaultdict(list)

    def filter_empty(l):
        return filter(lambda x: x != '', l)
    
    def get_coordinates(row):
        ret = [row[k] for k in ['# Chromosome', 'Genomic Position Start', 'Genomic Position End']]
        assert type(ret[0]) is str, "Chromosome is not str"
        assert type(ret[1]) is int, "Start is not int"
        assert type(ret[2]) is int, "End is not int"
        return ret

    def get_phenotypes(row):
        for special_token in [SEP_PHENO, SEP_INHE]:
            assert special_token not in row['Phenotypes'], f"Phenotypes str contains special token '{special_token}'"
        phenotypes = row['Phenotypes'].split('; ')
        phenotypes = filter_empty(phenotypes)
        return list(set(phenotypes))

    def get_all_symbols(row):
        symbols = row['Gene Symbols'].split(', ')
        symbols.append(row['Approved Symbol'])
        symbols = filter_empty(symbols)
        return list(set(symbols))

    for i, row in df.iterrows():
        coors = get_coordinates(row)
        all_symbols = get_all_symbols(row)
        phenotypes = get_phenotypes(row)
        phenotypes_inhe_patts = [(p, search_patterns(p)) for p in phenotypes]

        for symbol in all_symbols:
            key = (*coors, symbol)
            if key in db:
                print(f"WARNING: {key} has multiple row in OMIM")
            db[key].extend(phenotypes_inhe_patts)

    for key in db:
        if len(db[key]) == 0:
            print(f"WARNING: {key} has no phenotypes")
            db[key].append(('None', UNKNOWN))
    
    chr, s, e, gene_symbols = zip(*db.keys())
    phenotype_inhes = [SEP_PHENO.join([f"{t1}{SEP_INHE}{t2}" for t1, t2 in ps]) for ps in db.values()]

    df = pd.DataFrame({
        'chr': chr,
        'start': s,
        'end': e,
        'gene_symbol': gene_symbols,
        'phenotypes': phenotype_inhes
    })
    df.to_csv("omim.csv", index=False)

if __name__ == "__main__":
    typer.run(main)