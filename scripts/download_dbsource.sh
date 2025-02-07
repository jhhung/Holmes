###### script for downloading nearly all GRCh38 dbsources ######
set -ex

holmes_base_dir=$(dirname $0)/..
holmes_base_dir=$(realpath $holmes_base_dir)
default_output_dir=$holmes_base_dir/data/db
output_dir=${1:-$default_output_dir}

dir_no_ver=$output_dir/NoVersion
dir37=$output_dir/GRCh37
dir38=$output_dir/GRCh38

mkdir -p $dir_no_ver $dir37 $dir38

############################################
######### No version difference db #########
############################################
cd $dir_no_ver

# need download
ncbigene=https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
uniprot=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz

### You shoud request access and download OMIM genemap2.txt under $output_dir/NoVersion/

## get data ##
wget $ncbigene
wget $uniprot
gunzip *.gz

# #############################
# ######### GRCh37 db #########
# #############################
# cd $dir37

# # need download
# # NOTE: 1kg grch37 vcf contains ',' in ALT, need a normalization step to split them
# kg_37=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
# function norm_1kg() {
#     vcf=$1
#     out=$2
#     ref=$3

#     # Normalize
#     bcftools norm \
#         --atomize \
#         --fasta-ref $ref \
#         --check-ref w \
#         --multiallelics - \
#         --output-type z \
#         --output $out \
#         $vcf
# }
# clinvar_37=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar_20240107.vcf.gz
# ensembl_37=https://ftp.ensembl.org/pub/grch37/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
# refseq_37=https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz
# fasta_37=http://ftp.ensembl.org/pub/grch37/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
# dvd_37=https://deafnessvariationdatabase.org/public/releases/v9/DVD.r9.2021-01-04.download.vcf.bgz

# # save urls in file
# coverage_37_url=https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1/coverage/genomes/gnomad.genomes.coverage.summary.tsv.bgz
# function get_gnomad_37_urls() {
#     for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do
#         echo https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.$chr.vcf.bgz
#     done
#     # v2 genome has no chrY, so we use exome
#     echo https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.Y.vcf.bgz
# }

# ## get data
# wget $clinvar_37
# wget $ensembl_37
# wget $refseq_37
# wget $fasta_37
# wget $kg_37
# wget $dvd_37
# gunzip *.gtf.gz # unzip gtf
# norm_1kg \
#     ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz \
#     ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.normed.vcf.gz \
#     Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

# echo $coverage_37_url > coverage_37_url.txt
# for url in $(get_gnomad_37_urls); do
#     echo $url >> gnomad_37_urls.txt
# done

#############################
######### GRCh38 db #########
#############################
cd $dir38

# need download
kg_38=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz
clinvar_38=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20240107.vcf.gz
ensembl_38=https://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz
refseq_38=https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

# save urls in file
coverage_38_url=https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz
function get_gnomad_38_urls() {
    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
        echo https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr$chr.vcf.bgz
    done
}

## get data
wget $kg_38
wget $clinvar_38
wget $ensembl_38
wget $refseq_38
gunzip *.gtf.gz # unzip gtf

echo $coverage_38_url > coverage_38_url.txt
for url in $(get_gnomad_38_urls); do
    echo $url >> gnomad_38_urls.txt
done

### Following is GRCh38 building workflow example:
# 1. download (this script)

# 2. annotate clinvar
# > python ./scripts/clinvar_vep_annotation.py \
#       -i $dir38/clinvar_20240107.vcf.gz \
#       -v <vep executable> \
#       -f $dir38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
#       -o $dir38/clinvar_20240107_annotated.vcf \
#       -t <num threads> \

# 3. build holmes db
# # --fasta is deprecated, don't need to specify
# > ./bin/database_builder \
#     -k $dir38/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz \
#     --clinvar $dir38/clinvar_20240107_annotated.vcf \
#     --coverage $dir38/coverage_38_url.txt \
#     --omim $dir_no_ver/genemap2.txt \
#     --ncbigene $dir_no_ver/Homo_sapiens.gene_info \
#     --dvd $dir38/DVD_fixed_hg38.vcf \
#     --uniprot $dir_no_ver/uniprot_sprot.xml \
#     --ensembl_gtf $dir38/Homo_sapiens.GRCh38.107.gtf.gz \
#     --refseq_gtf $dir38/GCF_000001405.40_GRCh38.p14_genomic.gtf \
#     --gnom_file $dir38/gnomad_38_urls.txt \
#     -o <output database dir> \
#     -t <num threads>