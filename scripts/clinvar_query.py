import argparse
import ftplib
import subprocess

NCBI_HOST = "ftp.ncbi.nlm.nih.gov"
CLINVAR_ENTRY = "pub/clinvar"


def check_executable(executable: str):
    # check exec are installed
    try:
        subprocess.run(
            [executable, "--version"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError:
        print(f"{executable} are not installed")
        exit(1)


def get_args():
    parser = argparse.ArgumentParser(description="Query ClinVar", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-i",
        "--id",
        type=str,
        required=True,
        help="query ID, either Variation ID or AlleleID",
    )
    parser.add_argument(
        "-a",
        "--allele_id",
        action="store_true",
        default=False,
        help="A swith option. Turn this on to search Allele ID instead of Variation ID",
    )
    parser.add_argument(
        "-d", "--date", type=str, default=None, help="Search certain file date. example 20230710"
    )
    parser.add_argument(
        "-g",
        "--grch",
        type=int,
        choices=[37, 38],
        default=38,
        help="Search GRCh version",
    )
    parser.add_argument(
        "--bcftools", type=str, default="bcftools", help="bcftools binary path"
    )
    parser.add_argument(
        "-r", "--region_hint", type=str, default=None, help="Region hint for bcftools"
    )
    parser.add_argument(
        "-s",
        "--search",
        action="store_true",
        default=False,
        help="bruce force search all files until found",
    )
    parser.add_argument(
        "-l",
        "--limit",
        type=int,
        default=10,
        help="search limit, only works with --search option",
    )
    return parser.parse_args()


def run_bcftools(arg: argparse.Namespace, query_expr: str, file: str):
    bcftools_args = [arg.bcftools, "view", "-H", "-i", query_expr]
    if arg.region_hint is not None:
        bcftools_args.extend(["-r", arg.region_hint])
    bcftools_args.append(file)
    print(f"Query result of {file}")
    complete = subprocess.run(
        bcftools_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if complete.returncode != 0:
        print("bcftools error")
        print(complete.stderr.decode("utf-8"))
        return False
    out = complete.stdout.decode("utf-8")
    if out:
        print(out)
        return True
    else:
        print("not found")
        return False


def main(args):
    # check exec are installed
    check_executable(args.bcftools)

    ftp = ftplib.FTP(host=NCBI_HOST, user="anonymous", passwd="", acct="")
    weekly_dir = f"{CLINVAR_ENTRY}/vcf_GRCh{args.grch}/weekly"
    ftp.cwd(weekly_dir)

    files: list[str] = []
    ftp.retrlines("NLST", files.append)
    files = list(filter(lambda f: f.endswith("vcf.gz") and "papu" not in f, files))
    print(f"ClinVar available weekly files count: {len(files)}")

    if args.date:
        target_file = f"clinvar_{args.date}.vcf.gz"
    else:
        target_file = "clinvar.vcf.gz"

    if target_file not in files:
        print(f"ClinVar vcf file `{target_file}` not found")
        exit(1)

    # run query
    query_expr = f'ALLELEID={args.id}' if args.allele_id else f"\"ID='{args.id}'\""

    if args.search:
        search_space = files.copy()
        search_space.remove("clinvar.vcf.gz")
        search_space.sort(reverse=True)

        count = 0
        for f in search_space[: args.limit]:
            count += 1
            file = f"https://{NCBI_HOST}/{weekly_dir}/{f}"
            if run_bcftools(args, query_expr, file):
                return

        print(f"Search limit reached, {count} files searched")
    else:
        file = f"https://{NCBI_HOST}/{weekly_dir}/{target_file}"
        found = run_bcftools(args, query_expr, file)
        if found and args.date is None: # found the variant and use the newest file
            if args.allele_id:
                print(f"URL: https://www.ncbi.nlm.nih.gov/clinvar?term={args.id}[AlleleID]")
            else:
                print(f"URL: https://www.ncbi.nlm.nih.gov/clinvar/variation/{args.id}")


if __name__ == "__main__":
    main(get_args())
