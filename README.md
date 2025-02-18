# Holmes

Repository of Holmes: High-Throughput Genetic Variant Interpretation with an Automated Sherloc Framework

## Informations

Holmes is a robust and efficient variant interpretation platform based on Sherloc framework.

Holmes leverages Sherloc-guided criteria, augmented with novel disease-specific modules, to classify genomic variants efficiently.

## Requirement

- [Ubuntu 22.04](https://ubuntu.com/download) OS or higher
- GNU [g++-12](https://gcc.gnu.org/gcc-12/) or higher
- [CMake 3.16.0](https://cmake.org/download/) or higher to build the C++ Core Module
- [Boost 1.74.0](https://www.boost.org/users/history/version_1_74_0.html) to build the C++ Core Module
- [Ensembl VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) release/107 to run C++ Core Module
- [python 3.7](https://www.python.org/downloads/) or higher to run Compound Heterozygosity Module
- [libcurl4-openssl-dev](https://github.com/curl/curl) is required to enable htslib to open remote files

### Clone this project

```bash
# `--recurse-submodules` is needed to clone all submodules
git clone --recurse-submodules https://github.com/jhhung/Holmes.git
```

### Install VEP

see [doc for VEP](docs/VEP.md)

## Usage

### Build the C++ Core Module

```bash
# config
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc-12 -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 -S ./ -B ./build/

# build
cmake --build ./build --config Release
```

### Build Holmes database

see [doc for Database](docs/Database.md)

### Install python env for Compound Heterozygosity Module

```bash
# create venv (optional)
python3 -m venv .venv
source .venv/bin/activate

# install requirements
pip install -r requirements.txt
```

### Run Holmes

#### Make input configs

DB config:

```bash
python scripts/make_db_config.py --help
```

VEP config:

```bash
python scripts/make_vep_config.py --help
```

Input config:

```bash
python scripts/input_helper.py --help
```

#### Run Core Module

```bash
# run raw executable
./bin/sherloc --help

# or run with a wrapper script
## run with no args to see usage
bash scripts/run_one.sh 
```
