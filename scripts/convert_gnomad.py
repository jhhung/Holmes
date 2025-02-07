import argparse
import pathlib
import subprocess
import shutil

def convert_one_file(binary, old_file_path, new_file_path):
    cmd = [
        binary,
        "--convert",
        "--gnom", old_file_path,
        "--output", new_file_path
    ]
    return subprocess.run(cmd)

def get_all_file_list(old_gnom_path, new_gnom_path):
    file_list_of_files = []
    file_txt_to_be_copy = []
    with open(old_gnom_path / 'gnomAD_list.txt', 'r') as f:
        file_list_of_files.extend(f.read().strip().split())
    file_txt_to_be_copy.append((
        old_gnom_path / 'gnomAD_list.txt',
        new_gnom_path / 'gnomAD_list.txt'
    ))
    
    file_list = []
    for flf in file_list_of_files:
        file_txt_to_be_copy.append((
            old_gnom_path / flf,
            new_gnom_path / flf
        ))
        with open(old_gnom_path / flf, 'r') as f:
            file_list.extend(f.read().strip().split())
    return file_list, file_txt_to_be_copy
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='convert_gnomad')
    parser.add_argument("-g", "--gnom_builder_bin",
        help="The gnom_archive_builder binary path",
        required=True)
    parser.add_argument("-i", "--input_database_dir",
        help="The directory that contains '/gnomAD_list.txt' file. "
             "This dir is also the 'database' entry in holmes input json file + '/gnomAD'. "
             "ex. /path/to/<database>/gnomAD",
        required=True)
    parser.add_argument("-o", "--output_database_dir",
        help="This dir should be a mirror to the '--input-database-dir' "
             "All new gnomAD arc file will be place in the same relative path.",
        required=True)
    args = parser.parse_args()

    binary = pathlib.Path(args.gnom_builder_bin)
    in_gnom_dir = pathlib.Path(args.input_database_dir)
    out_gnom_dir = pathlib.Path(args.output_database_dir)

    assert(binary.is_file())
    assert(in_gnom_dir.is_dir())
    assert(out_gnom_dir.is_dir())

    all_arc_files_postfix, file_txt_to_be_copy = get_all_file_list(in_gnom_dir, out_gnom_dir)

    for arc in all_arc_files_postfix:
        print(f"Processing {arc}...")
        out_arc = out_gnom_dir / pathlib.Path(arc)
        
        # create the dir of output arc file
        out_arc.parent.mkdir(parents=True, exist_ok=True)
        try:
            convert_one_file(binary, in_gnom_dir / arc, out_arc).check_returncode()
        except:
            print(f"Error when processing {in_gnom_dir / arc} and {out_arc}!!")
            exit(1)
    print("All old arc files are converted.")
    print("Copying meta file ('gnomAD_list.txt' etc)...")
    for old_f, new_f in file_txt_to_be_copy:
        shutil.copyfile(pathlib.Path(old_f).resolve(), pathlib.Path(new_f).resolve())
    print("Copying coverage.arc ...")
    out_cov = pathlib.Path(out_gnom_dir / "coverage" / "coverage.arc")
    out_cov.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(
        pathlib.Path(in_gnom_dir / "coverage" / "coverage.arc").resolve(),
        out_cov.resolve())
    print("All done.")
