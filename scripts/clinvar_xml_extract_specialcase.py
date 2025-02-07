import sys
import typer
import xml.etree.ElementTree as ET
from pathlib import Path
from dataclasses import dataclass
from tqdm import tqdm
from enum import Enum
import json

class GenomeVersion(str , Enum):
    GRCh37 = 'GRCh37'
    GRCh38 = 'GRCh38'

@dataclass
class VCFRecord:
    chrom: str
    pos: str
    ref: str
    alt: str

    def __str__(self):
        return f"{self.chrom}\t{self.pos}\t{self.ref}\t{self.alt}\tPopulation"

@dataclass
class ClinvarRecord:
    variation_id: str
    allele_id: str
    vcf_record: VCFRecord | None = None

    @classmethod
    def from_xml(cls, root: ET.Element, genome_version: GenomeVersion):
        results: list[cls] = []
        required_attribs_map = {
            'VariationID': 'variation_id',
            'AlleleID': 'allele_id',
        }
        record_cols = ['Chr', 'positionVCF', 'referenceAlleleVCF', 'alternateAlleleVCF']
        for simple_allele in root.iter('SimpleAllele'):
            
            if not all(attr in simple_allele.attrib for attr in required_attribs_map):
                # print(f"Skipping SimpleAllele: attr={simple_allele.attrib}")
                continue

            var_id_kwargs = {v: simple_allele.attrib[k] for k, v in required_attribs_map.items()}
            
            vcf_records = []
            for location in simple_allele.iter('Location'):
                for seq_loc in location.iter('SequenceLocation'):
                    if 'Assembly' in seq_loc.attrib and seq_loc.attrib['Assembly'] == genome_version.value:
                        # check if all columns are present
                        if all(col in seq_loc.attrib for col in record_cols):
                            vcf_records.append(VCFRecord(*[seq_loc.attrib[col] for col in record_cols]))

            assert len(vcf_records) <= 1, f"Multiple VCF records found for {var_id_kwargs}"

            results.append(cls(
                vcf_record=vcf_records[0] if len(vcf_records) > 0 else None,
                **var_id_kwargs
            ))
        return results
    
def xml_streamer(file: Path):
    with open(file, "r") as f:
        line: str
        while line := f.readline():
            if line.startswith('<VariationArchive'):
                lines = [line]
                while l := f.readline():
                    lines.append(l)
                    if l.startswith('</VariationArchive>'):
                        break
                yield "".join(lines)


app = typer.Typer(pretty_exceptions_enable=False)

@app.command()
def filter(
    input_xml: Path,
    keyword: str = 'outweigh',
):
    pass_count = 0
    for rep in tqdm(xml_streamer(input_xml), desc=f'filtering xml records by keyword: {keyword}'):
        if keyword in rep:
            pass_count += 1
            print(rep, end='')

    print(f"Total records containing keyword '{keyword}': {pass_count}", file=sys.stderr)


@app.command()
def convert(
    input_xml: Path,
    output_tsv: Path,
    genome_version: GenomeVersion = GenomeVersion.GRCh38
):
    # header of special case tsv
    with open(output_tsv, "w") as f:
        f.write("#CHROM\tPOS\tREF\tALT\tAttribute\n")
        vcf_records = []
        for rep in tqdm(xml_streamer(input_xml), desc=f'parsing xml records'):
            root = ET.fromstring(rep)
            records = ClinvarRecord.from_xml(root, genome_version)
            if len(records) == 0:
                print(f"No records found in {rep}")
            elif len(records) > 1:
                print(f"Multiple records found in {rep}")
            else:
                record = records[0]
                if record.vcf_record is None:
                    print(f"Variation ID {record.variation_id} does not has vcf record")
                else:
                    vcf_records.append(record.vcf_record)
        vcf_records.sort(key=lambda x: (int(x.chrom.replace('chr', '')), int(x.pos)))
        for record in vcf_records:
            f.write(f"{record}\n")
if __name__ == "__main__":
    app()