from pathlib import Path
import subprocess
import argparse

def get_args():
    parser = argparse.ArgumentParser("Making VEP config json file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_dir", required=True, help="Input database dir to compress")
    parser.add_argument("-o", "--output_dir", required=True, help="Output file path")
    parser.add_argument("-a", "--archive_compressor", required=True, help="archive_compressor binary path")
    parser.add_argument("-b", "--backend", choices=['gzip', 'zstd'], default='gzip', help="compression backend")
    return parser.parse_args()

def main(args):
    arc_comp = Path(args.archive_compressor).expanduser().resolve()

    input_dir = Path(args.input_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()

    if not input_dir.exists():
        raise FileNotFoundError("Input directory does not exist")

    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    for arc in input_dir.rglob("*.arc"):
        if arc.is_file():
            rel = arc.relative_to(input_dir)
            output_file = output_dir / rel
            output_file.parent.mkdir(parents=True, exist_ok=True)
            print(f"Compressing {arc} to {output_file}")
            complete = subprocess.run([
                arc_comp,
                "-i", str(arc),
                "-o", str(output_file),
                "-b", args.backend
            ])
            if complete.returncode != 0:
                raise RuntimeError(f"Failed to compress {arc} to {output_file}")

if __name__ == "__main__":
    args = get_args()
    main(args)
