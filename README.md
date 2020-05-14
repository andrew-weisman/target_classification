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
download_files "$target_website" "${project_dir}data/" "${working_dir}all_files_in_tree.txt" |& tee "${working_dir}download_files_out_and_err.txt"
```

Obtain the unique headers for the downloaded files:

```bash
get_unique_headers "${project_dir}data/" "${working_dir}all_files_in_tree.txt"
```

This creates the file `${project_dir}data/unique_headers.txt` containing the unique headers, which likely requires manual cleanup, as explained in the `get_unique_headers()` function. This includes manually creating a new file `${project_dir}data/unique_basenames.txt` in order to help fully differentiate between the different file formats.

Using these two files, we can now create an index file `${project_dir}data/format_mapping.txt` containing the format numbers (corresponding to the line numbers in `${project_dir}data/unique_headers.txt` and `${project_dir}data/unique_basenames.txt`) for each file in our dataset:

```bash
get_format_mapping "${project_dir}data/" "${working_dir}all_files_in_tree.txt"
```

Print out a text report containing useful properties of the files corresponding to each format:

```bash
ensure_format_consistency "${project_dir}data/" > "${working_dir}uniformity_check.txt"
```

This file should be reviewed in order to ensure that each file format does indeed differentiate the filetypes sufficiently. For example, it was using this file that we learned that we needed to create the `${project_dir}data/unique_basenames.txt` file in order to supplement `${project_dir}data/unique_headers.txt` for fully differentiating between file formats.

## Next up

* Use `${working_dir}uniformity_check.txt` and Mark's most recent email to write blocks in the `extract_data()` function (one of which is already there) to process the files of each format into TSV files that we can later read into a Pandas dataframe using the `load_tsv_files()` function in the `target_class_lib.py` Python library.
