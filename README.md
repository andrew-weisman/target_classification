# Workflow summary for TARGET classification

## Preliminary notes

* See the highly-commented libraries `target_class_lib.sh` and `target_class_lib.py` (same directory as this README) for more details on the following commands that reference these libraries; the document you're reading is just an overall summary of the commands.
* Additional documentation can be found at [this Collaborate site](https://collaborate.nci.nih.gov/x/kJHrDg).

## Workflow summary

Clone and set up this repository (that which holds this README file) on Biowulf:

```bash
mkdir /data/BIDS-HPC/private/projects/dmi
cd !!:1
git clone git@github.com:andrew-weisman/target_classification.git checkout
mkdir data # non-version-controlled files will be stored here
```

Get a compute node (e.g., `sinteractive --mem=20g`) and assign some variables:

```bash
project_dir="$(pwd)/"
working_dir="/home/weismanal/notebook/2020-04-08/scraping_target_site/"
target_website="https://target-data.nci.nih.gov/"
```

Load the Bash library:

```bash
source "${project_dir}checkout/target_class_lib.sh"
```

Crawl through the [TARGET website](https://target-data.nci.nih.gov)'s file tree and create empty file placeholders in the data directory:

```bash
# Note: This takes about 67 minutes to complete
get_file_placeholders "$target_website" "Controlled/ Public/" "${project_dir}data/" |& tee "${working_dir}get_file_placeholders_out_and_err.txt"
```

From the resulting file list, namely `all_files_in_tree.txt`, we need to determine which files we'd like to use in our dataset and then place the appropriate command in the `get_datafile_list()` function in the `target_class_lib.sh` library (the discussion resulting in our current set of 2261 datafiles is located at the [Collaborate page](https://collaborate.nci.nih.gov/x/kJHrDg)).

Download the files containing the data we'd like to include in our dataset:

```bash
download_files "$target_website" "${project_dir}data/" |& tee "${working_dir}download_files_out_and_err.txt"
```

Obtain the unique headers for the downloaded files:

```bash
get_unique_headers "${project_dir}data/"
```

This creates the file `${project_dir}data/unique_headers.txt` containing the unique headers, which likely requires manual cleanup, as explained in the `get_unique_headers()` function. This includes manually creating a new file `${project_dir}data/unique_basenames.txt` in order to help fully differentiate between the different file formats.

Using these two files, we can now create an index file `${project_dir}data/format_mapping.txt` containing the format numbers (corresponding to the line numbers in `${project_dir}data/unique_headers.txt` and `${project_dir}data/unique_basenames.txt`) for each file in our dataset:

```bash
get_format_mapping "${project_dir}data/"
```

Print out a text report containing useful properties of the files corresponding to each format:

```bash
# Note: This takes about 1-2 minutes to complete
ensure_format_consistency "${project_dir}data/"
```

This file should be reviewed in order to ensure that each file format does indeed differentiate the filetypes sufficiently. For example, it was using this file that we learned that we needed to create the `${project_dir}data/unique_basenames.txt` file in order to supplement `${project_dir}data/unique_headers.txt` for fully differentiating between file formats.

### To edit before putting into workflow summary above

#### Create a file called all_best_gene_names.txt in the data directory that writes out the "best" gene name in each line of each datafile

```bash
get_best_gene_names_from_all_files "${project_dir}data/"
```

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
awk '{print $1}' "${project_dir}data/all_best_gene_names.txt" | sort -u | awk '{printf("|%s|\n", toupper($1))}' > "${project_dir}data/unique_best_gene_names_uppercase.txt"
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

## Next up

* Use `${project_dir}data/uniformity_check.txt` and Mark's most recent email to write blocks in the `extract_data()` function (one of which is already there) to process the files of each format into TSV files that we can later read into a Pandas dataframe using the `load_tsv_files()` function in the `target_class_lib.py` Python library.
* In doing this, I need to pick up with hopefully hearing back from the Ensembl team about how to perform a POST-like operation using their xrefs method, and then implement it on the remaining unknown gene names in the Python version of the global lookup table... after that I will have done the best that I can and can start extracting the data from the datafiles using extract_data() or a Python version of that (bottom line was I needed to try my best to obtain an Ensembl ID for each of the 128,610 unique gene names!)
* Don't forget to run the overall processes by Mark, in particular, see my partial Collaborate page to him
