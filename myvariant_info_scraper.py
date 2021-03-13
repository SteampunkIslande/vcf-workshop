import requests
import sys
from datetime import datetime

inf = sys.stdin # inf is the input for this module, eiter a file or stdin

error_message = []

def variant_fields_generator():
    """
    Returns the next non-empty ID field
    """
    while True:
        try:
            _data = next(inf) # Read next line of input
        except StopIteration:
            return
        chrom,pos,_id,ref,alt = _data.split()
        if _id != "." and _id != "":
            yield { "CHROM":chrom,
                    "POS":pos,
                    "ID":_id,
                    "REF":ref,
                    "ALT":alt} # Only return line if _data is not empty
            

def dump_variant_annotations(variants :list):
    log = open(f"log - {datetime.strftime(datetime.now(),'%d-%m-%Y - %H:%M:%S')}","w+")
    variant_ids = ",".join([v["ID"] for v in variants]) # Retrieve a list of variant's ids, in the form of a comma-separated string
    resp = requests.post("http://myvariant.info/v1/variant",{
                                                            "fields":"gnomad_exome",
                                                            "ids":variant_ids,
                                                            "assembly":"hg19"})

    assert None not in variants, "Error, when calling dump_variant_annotations, make sure that variant list has no None"
    assert len(variants)<=1000
    requested_variants_count = len(variants)
    annotated_variants_count = 0
    if resp.status_code == 200:
        try:
            data = resp.json()
            for variant_annot,variant in zip(data,variants):
                try:
                    af = variant_annot["gnomad_exome"]["af"]["af"]

                    chrom,pos,_id,ref,alt = variant.values()

                    print(f"{chrom}\t{pos}\t{ref}\t{alt}\t{af}")
                    annotated_variants_count += 1
                except KeyError:
                    error_message.append(f"Cannot find allele frequency for variant {variant}")
                    continue

        except ValueError:
            error_message.append(f"No valid data found for variants {variant_ids} !")
    log.write(f"###LOG : {annotated_variants_count} out of {len(variants)} variants requests were successful")



def main():
    """
    This script takes a list of variants (CHROM, POS, ID, REF, ALT in this order) from stdin, and queries for each non-empty variant ID its gnomAD_exome information.
    The output (dumped to stdout) is a TSV that can be used by bcftools to annotate a VCF (see snakefile for more info)
    """
    all_variants = variant_fields_generator()
    variants_chunk = [None] * 1000
    idx = 0
    for variant in all_variants:
        variants_chunk[idx%1000] = variant
        idx+=1
        if idx%1000 == 0:# Every 1000 variants
            dump_variant_annotations(variants_chunk) # Retrieve variant's annotation and dumps it to stdout

            variants_chunk = [None] * 1000 # Reset variants_chunk

    if idx % 1000 != 0:# This is one way to tell that we exited the loop without requesting variant data
        dump_variant_annotations([v for v in variants_chunk if v is not None]) # Request variant data, making sure that there is no None variant


if __name__ == "__main__":
    exit(main())