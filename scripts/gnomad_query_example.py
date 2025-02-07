import requests
import argparse

gnomad_endpoint = "https://gnomad.broadinstitute.org/api"
gnomad_query_str = """
{
    variant(variantId: "%s", dataset: gnomad_r4){
        variant_id
        joint {
            ac
            an
            faf95 {
                popmax
                popmax_population
            }
        }
    }
}
"""

def query_gnomad(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
):
    chrom = chrom.replace("chr", "")
    variant = f"{chrom}-{pos}-{ref}-{alt}"
    res = requests.post(gnomad_endpoint, json={'query': gnomad_query_str % variant})
    if not res.ok:
        print(f"Failed to query gnomad for {variant}, response: {res.text}")
        return None
    data = res.json()
    variant_data = data['data']['variant']
    if variant_data is None:
        print(f"Empty variant entry for {variant}, response: {res.text}")
        return None
    return variant_data
    

def main(
    args
):
    res = query_gnomad(args.chrom, args.pos, args.ref, args.alt)
    print(res)

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("chrom", type=str)
    argparser.add_argument("pos", type=int)
    argparser.add_argument("ref", type=str)
    argparser.add_argument("alt", type=str)
    args = argparser.parse_args()
    main(args)