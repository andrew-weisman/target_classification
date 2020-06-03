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
get_file_placeholders "$target_website" "Controlled/ Public/" "${project_dir}data/" |& tee "${working_dir}get_file_placeholders_out_and_err.txt" # creates ${project_dir}data/all_files_in_tree.txt
```

From the resulting file list, `all_files_in_tree.txt`, we need to determine which files we'd like to use in our dataset and then place the appropriate command in the `get_datafile_list()` function in the `target_class_lib.sh` library (the discussion resulting in our current set of 2261 datafiles is located at the [Collaborate page](https://collaborate.nci.nih.gov/x/kJHrDg)).

Download the files containing the data we'd like to include in our dataset:

```bash
download_files "$target_website" "${project_dir}data/" |& tee "${working_dir}download_files_out_and_err.txt"
```

Obtain the unique headers for the downloaded files:

```bash
get_unique_headers "${project_dir}data/" # creates ${project_dir}data/unique_headers.txt

# This step also requires you create ${project_dir}data/unique_basenames.txt as well, as described below and in the function comments
```

This creates the file `${project_dir}data/unique_headers.txt` containing the unique headers, which likely requires manual cleanup, as explained in the `get_unique_headers()` function. This includes manually creating a new file `${project_dir}data/unique_basenames.txt` in order to help fully differentiate between the different file formats.

Using these two files, we can now create an index file `${project_dir}data/format_mapping.txt` containing the format numbers (corresponding to the line numbers in `${project_dir}data/unique_headers.txt` and `${project_dir}data/unique_basenames.txt`) for each file in our dataset:

```bash
get_format_mapping "${project_dir}data/" # creates ${project_dir}data/format_mapping.txt
```

Print out a text report in the data directory called `uniformity_check.txt` containing useful properties of the files corresponding to each format:

```bash
# Note: This takes about 1-2 minutes to complete
ensure_format_consistency "${project_dir}data/" # creates ${project_dir}data/uniformity_check.txt
```

This file should be reviewed in order to ensure that each file format does indeed differentiate the filetypes sufficiently. For example, it was using this file that we learned that we needed to create the `${project_dir}data/unique_basenames.txt` file in order to supplement `${project_dir}data/unique_headers.txt` for fully differentiating between file formats.

Using the report contained in `${project_dir}data/uniformity_check.txt`, write a function in `${project_dir}checkout/target_class_lib.sh` called `get_best_gene_name_from_file()` that extracts the best gene name for each line in each type of file, and subsequently generate a master list, `${project_dir}data/all_best_gene_names.txt`, of the "best" gene names contained in every line of every file in the dataset:

```bash
# This takes about 3 minutes
get_best_gene_names_from_all_files "${project_dir}data/" # creates ${project_dir}data/all_best_gene_names.txt
```

Get the unique gene names in all these files:

```bash
# This takes about 8 minutes
awk '{print $1}' "${project_dir}data/all_best_gene_names.txt" | sort -u | awk '{print $1}' | sort -u > "${project_dir}data/unique_best_gene_names.txt" # 128,612 of these
```

Split this file into Ensembl IDs and other gene names:

```bash
grep -E "^ENSG[0-9]{11}$" "${project_dir}data/unique_best_gene_names.txt" > "${project_dir}data/unique_ensembl_ids.txt" # 73,615 of these
grep -E -v "^ENSG[0-9]{11}$" "${project_dir}data/unique_best_gene_names.txt" > "${project_dir}data/unique_other_names.txt" # 54,997 of these
```

Create a gene lookup table generation using HGNC's Biomart:

 1. Go [here](https://biomart.genenames.org)
 1. Click on the [Gene mart](https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart)
 1. Deselect "Status" and "Approved name" and select "Ensembl gene ID" (URL then changes to [this](https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart&attributes=hgnc_gene__hgnc_gene_id_1010%2Chgnc_gene__approved_symbol_1010%2Chgnc_gene__ensembl_gene__ensembl_gene_id_104))
 1. Select "Go" at the bottom
 1. URL doesn't seem to change; it is still [this](https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart&attributes=hgnc_gene__hgnc_gene_id_1010%2Chgnc_gene__approved_symbol_1010%2Chgnc_gene__ensembl_gene__ensembl_gene_id_104)
 1. Click on "Download data", which downloads a file called `results.txt` containing tab-separated data
 1. Copy this local file to `${project_dir}data/gene_lookup_table.txt`

Extract the datafiles into consistently-formatted TSV files to be subsequently read into Pandas dataframes in Python (this also writes the project and file metadata to a JSON file):

```bash
# This takes about 5 minutes
# This creates ${project_dir}data/tsv_files and all the .tsv files within it, as well as ${project_dir}data/metadata.json
extract_data "https://target-data.nci.nih.gov/" "${project_dir}data/" |& tee "${working_dir}extract_data_out_and_err.txt"
```

## Next up

* Perhaps read the TSV files into a Pandas dataframe using the `load_tsv_files()` function in the `target_class_lib.py` Python library... probably read them all into their own Pandas dataframes first
* Don’t forget to address names with slashes in them, see e.g. the comment in target_class_lib.py; basically, we can just search the output for WARNING perhaps
* Look for duplicate data and probably delete the ones that are subsets of the other, e.g.
* Use other names to help identify other synonyms... go through get_uniform_version_of_file() function and search for "name"
* At the end of gene name resolution, confirm all IDs are in the current Ensembl version
* Explore the data we're given, e.g., in format #8, the status is LOW DATA, should we include these?
* Ask Patee about being able to trust that all datafiles in Public/ are tumor data of the type specified simply by the directory name
