# Database

## Database Building Workflow

1. Download / copy source files
2. Annotate ClinVar
3. Build Holmes Database
4. (Optional) Build Holmes VEP cache

## Database source files

```bash
# this scripts is for GRCh38 version
## The GRCh37 part is commented out and will not be executed
bash scripts/download_dbsource.sh <output_dir>
```

**Important Note:** You have to request the access to the `genemap2.txt` file from OMIM your self and place it under `<output_dir>/NoVersion/`

## ClinVar annotation

After downloading the above files, you can use them to perform VEP annotation on the ClinVar VCF.

This step will generate clinvar_annotated.vcf, which is the file that should be fed into the next step.

```bash
python scripts/clinvar_vep_annotation.py --help

usage: clinvar_vep_annotation.py [-h] -i INPUT -v VEP -f FASTA -o OUTPUT [-t THREAD] [--grch37]

Annotate ClinVar VCF with VEP

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input VCF file
  -v VEP, --vep VEP     VEP binary
  -f FASTA, --fasta FASTA
                        Fasta ref file
  -o OUTPUT, --output OUTPUT
                        Output vep annotation file
  -t THREAD, --thread THREAD
                        fork of VEP
  --grch37              use GRCh37 assembly
```

## Database builder

After completing the above steps, you can start building the database for Holmes.

```bash
./bin/database_builder -h
Allowed options:
  -h [ --help ]                show help message
  -k [ --k ] arg               1000 Genomes Project file
  -c [ --clinvar ] arg         VEP annotated Clinvar file
  -v [ --coverage ] arg        Coverage file
  --omim arg                   Gene info: OMIM genemap file
  --ncbigene arg               Gene info: NCBI gene file
  -d [ --dvd ] arg             DVD database file
  -u [ --uniprot ] arg         UniProt database file
  -f [ --fasta ] arg           Reference Fasta file
  --ensembl_gtf arg            Ensembl GTF file
  --refseq_gtf arg             Refseq GTF file
  -g [ --gnom_file ] arg       gnomAD url list file (each chr url/path per 
                               line)
  -o [ --output ] arg (=/tmp/) Output directory
  -t [ --thread ] arg (=4)     Thread num for parallel building GnomAD
```

## VEP Cache builder

If you want to use the VEP cache to speed up VEP annotation, you can use `vep_cache_builder`.

```bash
./bin/vep_cache_builder 
Allowed options:
  -h [ --help ]                         show help message
  -i [ --input ] arg                    input vep annotation file
  -o [ --output ] arg (=/tmp/vep_cache/)
                                        output vep cache archive dirname, this 
                                        dir will end up containing all the 
                                        `chr.arc` file
  -c [ --config ] arg                   vep config file for running vep. If not
                                        specified, will assume the `--input` 
                                        file is already annotated.
  -t [ --thread_num ] arg (=1)          Thread num (for VEP)
  --grch37                              Use grch37 coordinate
```

The input VCF can be any VCF file, such as sites with 1000 Genomes AF > 0.05.

If no config file is specified, it is assumed that the input VCF has already undergone VEP annotation.

If a VEP config is provided in the config file, VEP will be run according to the specified parameters to annotate the input VCF.

It is recommended that the config file be the same as the one used with the `--vep_config` option when running Holmes.
