### To edit before putting into workflow summary above

#### Create a gene lookup table generation (for transforming between HGNC symbols and Ensembl IDs)

1. Go [here](https://biomart.genenames.org)
1. Click on the [Gene mart](https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart)
1. Deselect "Status" and "Approved name" and select "Ensembl gene ID" (URL then changes to [this](https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart&attributes=hgnc_gene__hgnc_gene_id_1010%2Chgnc_gene__approved_symbol_1010%2Chgnc_gene__ensembl_gene__ensembl_gene_id_104))
1. Select "Go" at the bottom
1. URL doesn't seem to change; it is still [this](https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart&attributes=hgnc_gene__hgnc_gene_id_1010%2Chgnc_gene__approved_symbol_1010%2Chgnc_gene__ensembl_gene__ensembl_gene_id_104)
1. Click on "Download data", which downloads a file called `results.txt` containing tab-separated data
1. Copy this local file to `${project_dir}data/gene_lookup_table.txt`

#### Prepare text files by making their contents uppercase

```bash
tail -n +2 "${project_dir}data/gene_lookup_table.txt" | awk '{printf("%s\t|%s|\n", $1, toupper($2))}' > "${project_dir}data/lookup_symbol_uppercase.txt"
tail -n +2 "${project_dir}data/gene_lookup_table.txt" | awk '{printf("%s\t|%s|\n", $1, toupper($3))}' > "${project_dir}data/lookup_id_uppercase.txt"
```

#### Create a file in the data directory called whether_genes_are_known.txt containing three columns for each unique best gene name in the entire dataset: (1) the number of matches of symbols from the lookup table, (2) the number of matches of the IDs from the lookup table, (3) the unique best gene name itself

```bash
show_whether_genes_are_known "${project_dir}data/" "$working_dir"
```

#### For each of the four types of gene names identified in show_whether_genes_are_known() above (i.e., how they do in the HGNC lookup table), do as best as possible using the HGNC lookup table to get the Ensembl IDs, creating a full-size lookup table with blank third columns if really no Ensembl ID is currently known

```bash
create_partial_global_lookup_table "${project_dir}data/"
```

#### Write to file the gene names that need to identified using REST

```bash
awk 'NF==2{print $2}' "${project_dir}data/partial_global_lookup_table.txt" | sort -u > "${project_dir}data/to_lookup.txt"
```

#### Use the main.ipynb Jupyter notebook

* Run Steps 1 through 3: Perform a POST lookup on the unknown symbol names and add them to a Python version of the global lookup table
