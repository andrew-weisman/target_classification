#!/bin/bash

filename="/data/BIDS-HPC/private/projects/dmi/data/tree/Public/ALL/mRNA-seq/Phase1/L3/expression/BCCA/HS0825.gene.quantification.txt"
weblink="https://target-data.nci.nih.gov/Public/ALL/mRNA-seq/Phase1/L3/expression/BCCA/HS0825.gene.quantification.txt"

# awk -v filename=$filename -v weblink=$weblink '{
#     if (NR==1)
#         printf("%s\t%s\t%s\t%s\t%s\n", "gene-pretty", "gene-ugly", toupper($4), "filename", "weblink")
#     else {
#         split($1, arr, "|")
#         printf("%s\t%s\t%s\t%s\t%s\n", toupper(arr[1]), toupper(arr[2]), $4, filename, weblink)
#     }
# }' $filename

awk '{
    if (NR==1)
        printf("%s\t%s\t%s\n", "gene-pretty", "gene-ugly", toupper($4))
    else {
        split($1, arr, "|")
        printf("%s\t%s\t%s\n", toupper(arr[1]), toupper(arr[2]), $4)
    }
}' $filename

filenames="/data/HS9999.gene.quantification.txt /data/HS0825.gene.quantification.txt"
weblinks="https://target-data.nci.nih.gov/HS9999.gene.quantification.txt https://target-data.nci.nih.gov/HS0825.gene.quantification.txt"
