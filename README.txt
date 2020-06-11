gencode.v22.annotation.gtf is from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz (main annotation file at https://www.gencodegenes.org/human/release_22.html)
gencode.gene.info.v22.tsv is from https://api.gdc.cancer.gov/data/b011ee3e-14d8-4a97-aed4-e0b10f6bbe82 (GDC.h38 GENCODE TSV file at https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files)
These have the exact same, unique Ensembl IDs (60,483 of them) based on:
  "tail -n +2 gencode.gene.info.v22.tsv | awk '{print $1}' | sort -u"
  "tail -n +6 gencode.v22.annotation-just_genes.gtf | awk '{split($10,arr,"\""); print arr[2]}' | sort -u"
