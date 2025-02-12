# VEP installation steps

## Clone repo and install dependencies

```bash
# clone repo
git clone https://github.com/Ensembl/ensembl-vep.git -b release/107

# install dependencies
sudo apt update && sudo apt install \
    libarchive-zip-perl \
    libdbi-perl \
    libdbd-mysql-perl \
    libset-intervaltree-perl \
    libjson-c-dev \
    libperlio-gzip-perl \
    libmodule-build-perl \
    libwww-perl \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    zip \
    curl \
    tabix \
    liblist-moreutils-perl
```

## Install VEP

```bash
cd ensembl-vep

# install VEP software
perl INSTALL.pl --NO_TEST --PLUGINS all -n -a pa

# download and bgzip Human genome from Ensembl
wget ftp://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bgzip Homo_sapiens.GRCh38.dna.primary_assembly.fa

# download and extract indexed merged VEP cache (contains both Ensembl and RefSeq annotation)
cd $HOME/.vep
curl -O ftp://ftp.ensembl.org/pub/release-107/variation/indexed_vep_cache/homo_sapiens_merged_vep_107_GRCh38.tar.gz
tar xzf homo_sapiens_merged_vep_107_GRCh38.tar.gz

# download data for plugins
mkdir $HOME/.vep/vep_plugin_data
cd $HOME/.vep/vep_plugin_data

## MaxEntScan
wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz
gunzip fordownload.tar.gz
tar -xvf fordownload.tar

## pLI
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/107/pLI_values.txt > $HOME/.vep/Plugins/pLI_values.txt

## REVEL
wget 'https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip'
### Proprocess revel file for VEP plugin (steps are from VEP Plugin github)
unzip revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv
zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz
```
