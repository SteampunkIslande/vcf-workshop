rule get_gnomadex_annotations:
   input:
      "{sample}.ann.vcf"
   output:
      "gnomadex.tsv"
   shell:
      "cat {input} | grep -v '#' | cut -f1-5 | python 'myvariant_info_scraper.py' > gnomadex.tsv"

rule compress_and_index_annotation_file:
   input:
      "gnomadex.tsv"
   output:
      "gnomadex.tsv.gz"
   shell:
      "cat {input} | bgzip > {output} && tabix -s1 -b2 -e2 {output}"

rule annotate_file_with_gnomadex_annotations:
   input:
      "gnomadex.tsv.gz",
      "gnomadex_annot",
      "{sample}.ann.vcf"
   output:
      "{sample}.ann.gnomadex.vcf"
   shell:
      "{input[0]} -c 'CHROM,POS,REF,ALT,GNOMADEX' -h {input[1]} 'input{2}' > {output}"

rule all:
   input:
      "{sample}.ann.gnomadex.vcf"
