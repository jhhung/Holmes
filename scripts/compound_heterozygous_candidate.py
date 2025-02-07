import pandas as pd
import sys
from pathlib import Path
import argparse

SCORES_TABLE = {"155": 1.0, "154": 1.5, "156": 1.0}

AF_NOT_HIGH_RULES = [
    "161",
    "160",
    "167",
    "164",  # somewhat high
    "101",
    "135",  # pathogenic range
]


def modify_score(g_all: pd.DataFrame):  # for groupby apply'
    if g_all.shape[0] == 0:
        return g_all

    gene_sym = g_all["gene_symbol"].iloc[0]
    g = g_all[g_all["genotype"] == "hete"]
    
    # filter AR gene
    g = g[g['inheritance_pattern'].str.split(':').str[0] == 'AR']

    # should have atleast 2 hete variant in the same gene
    if g["HGVSg"].unique().shape[0] < 2:
        return g_all

    # TODO: This is for rule#156, only variants that have af not above 'High' can use 156
    g_rules = g["rule"].str.split(";")
    af_not_high_filter = g_rules.transform(
        lambda rs: any(map(lambda r: r.split(":")[0] in AF_NOT_HIGH_RULES, rs))
    )
    g_not_high = g[af_not_high_filter]
    print(g_all.index)
    print(g_not_high.index)
    # input()

    max_score_in_group = g["score"].max()
    if max_score_in_group > 3.9:  # pathogenic / likely pathogenic
        use_rule = ["154"]
    elif max_score_in_group > -3.9:  # vus
        use_rule = ["155"]
    else:
        print(
            f"No variant in this gene ({gene_sym}) has high enough score to meet the criteria"
        )
        return g_all

    # TODO: If variant has phased info, rule#156 can be further aaplied to in trans variant
    #   in cis variant will be excluded, too

    print(f"Used rules:")
    add_score = 0.0
    for r in use_rule:
        additional_score = SCORES_TABLE[r]
        add_score += additional_score
        print(f"\tEV{r.rjust(4, '0')}, additional score: {additional_score}")

    print_col = ["chr", "pos", "ref", "alt", "score", "rule"]
    print(f"Compound Heterozygous rules can apply to these variant:")
    print(g_not_high[print_col])
    print("-" * 80)

    print(g_not_high.index)

    def extend_rule(rs):
        return ";".join(rs + list(map(lambda x: f"{x}:E", use_rule)))

    g_all.loc[g_not_high.index, "rule"] = (
        g_not_high["rule"].str.split(";").transform(extend_rule)
    )
    g_all.loc[g_not_high.index, "score"] = g_not_high["score"] + add_score
    new_cons = []
    for s in g_all.loc[g_not_high.index, "score"]:
        """
        if      ( score >= 5. )  consequence_idx = 0; // pathogenic
          else if ( score >= 4. )  consequence_idx = 1; // likely pathogenic
          else if ( score <= -5. ) consequence_idx = 2; // benign
          else if ( score <= -4. ) consequence_idx = 3; // likely benign
          else                     consequence_idx = 4; // uncertain significance
        """
        if s >= 4.9:
            new_cons.append("pathogenic")
        elif s >= 3.9:
            new_cons.append("likely_pathogenic")
        elif s <= -4.9:
            new_cons.append("benign")
        elif s <= -3.9:
            new_cons.append("likely_benign")
        else:
            new_cons.append("uncertain_significance")
    g_all.loc[g_not_high.index, "consequence"] = new_cons

    return g_all


def analyzer(args):
    res_df_file = Path(args.file)
    output = args.output
    if output is None:
        output = res_df_file.with_suffix(".compound_hete")
    else:
        output = Path(output)

    df = pd.read_csv(res_df_file, sep="\t")
    # FIXME: try to reserve the original df to output?
    # mod_df = df.dropna(subset=["gene_symbol"], axis=0)
    # ...
    # df.loc[mod_df.index, :] = mod_df
    df = df.dropna(subset=["gene_symbol"], axis=0)
    

    df = df.groupby("gene_symbol").apply(modify_score)

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compound_Heterozygous")
    parser.add_argument("-f", "--file", required=True, help="The sherloc output file")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="if the output are not specified, the result will be output with suffix '.compound_hete'",
    )
    args = parser.parse_args()

    analyzer(args)
