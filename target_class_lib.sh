#!/bin/bash

# Note documentation for this library is at https://collaborate.nci.nih.gov/x/kJHrDg

# Run this script like, e.g.:
#
#   (source /data/BIDS-HPC/private/projects/dmi/checkout/target_class_lib.sh; get_file_placeholders "https://target-data.nci.nih.gov/" "Controlled/ Public/" "/data/BIDS-HPC/private/projects/dmi/data/" |& tee -a all_dirs.txt)
#   (source /data/BIDS-HPC/private/projects/dmi/checkout/target_class_lib.sh; get_file_placeholders "https://target-data.nci.nih.gov/Public/RT/mRNA-seq/" "L3/ METADATA/" "/data/BIDS-HPC/private/projects/dmi/data/" |& tee -a all_dirs.txt)
#   (source /data/BIDS-HPC/private/projects/dmi/checkout/target_class_lib.sh; download_files "https://target-data.nci.nih.gov/" "/data/BIDS-HPC/private/projects/dmi/data/" "/home/weismanal/notebook/2020-04-08/scraping_target_site/all_files_in_tree.txt" |& tee -a downloads.txt)
#

# Main recursive function to process all the directories and files at the input URL; probably don't call this function directly
crawl() {

    # Define function parameters
    local url_orig=$1
    local directories=$2
    datadir=$3

    # Constants
    index_file="index.html"
    tmp_file="tmp.txt"

    # For each directory in the directory list...
    for dir in $directories; do

        # Since we loop over all the directories in the input directory, we need to reset the URL each time since we redefine it for simplicity later
        url=$url_orig

        # Create the directory $dir if it doesn't already exist
        if [ ! -d "$dir" ]; then
            mkdir "$dir"
        else
            echo "WARNING: Directory '$dir' already exists"
        fi

        # Enter the directory $dir
        cd "$dir" || (echo "ERROR: Cannot cd into directory '$dir'"; exit 1)
        echo "In directory: $(pwd | awk -v FS="$datadir" '{print $2}')"

        # Create the file index.html in the directory $dir
        url="$url$dir"
        if [ ! -f "$index_file" ]; then
            wget "$url" 2>> "$datadir/wget_err-crawl.txt" || (echo "ERROR: Problem requesting file '$index_file' from URL '$url'"; exit 2)
        else
            echo "WARNING: File '$dir$index_file' from URL '$url' already exists"
        fi

        # Output to tmp.txt the relevant parts of index.html
        IFS_old=$IFS
        IFS=$'\n'
        awk 'BEGIN{ doprint=0 } { if ($0~"indexbreakrow") { if (doprint) doprint=0; else doprint=1 } else { if (doprint && !($0~"PARENTDIR")) print } }' "$index_file" > $tmp_file || (echo "ERROR: Problem getting the relevant chunk of file '$index_file'"; exit 3)
        IFS=$IFS_old

        # Get lists of the directories, files, and file sizes at the current web address
        directories=$(grep "theme/icons/folder.png" $tmp_file | awk -v FS="indexcolname" '{print $2}' | awk -v FS="href=\"" '{print $2}' | awk -v FS="\">" '{print $1}')
        files=$(grep -v "theme/icons/folder.png" $tmp_file | awk -v FS="indexcolname" '{print $2}' | awk -v FS="href=\"" '{print $2}' | awk -v FS="\">" '{print $1}')
        file_sizes=$(grep -v "theme/icons/folder.png" $tmp_file | awk -v FS="indexcolsize" '{print $2}' | awk -v FS="\">" '{print $2}' | awk -v FS="<" '{print $1}')

        # Create empty files, appended with the appropriate file sizes, for all the files at the current URL
        ifile=0
        for file in $files; do
            ifile=$((ifile+1))
            file_size=$(echo $file_sizes | awk -v ifile=$ifile '{print $ifile}') # probably not double-quotting $file_sizes on purpose here
            touch "${file}____$file_size"
        done

        # Recurse into the directories at the current URL
        crawl "$url" "$directories" "$datadir"

        # Delete the temporary file and go back to the parent directory
        rm -rf $tmp_file
        cd ..

    done

}


# Main function to call directly in order to create the file placeholders
get_file_placeholders() {

    # Sample calls:
    #   get_file_placeholders "https://target-data.nci.nih.gov/" "Controlled/ Public/" "/data/BIDS-HPC/private/projects/dmi/data/"
    #   get_file_placeholders "https://target-data.nci.nih.gov/Controlled/AML/mRNA-seq/L3/structural/" "BCCA/ NCI-Meerzaman/" "/data/BIDS-HPC/private/projects/dmi/data/"
    #   get_file_placeholders "https://target-data.nci.nih.gov/Public/RT/mRNA-seq/" "L3/ METADATA/" "/data/BIDS-HPC/private/projects/dmi/data/"

    # Parameters
    base_url=$1
    top_directories=$2
    datadir=$3

    # Constant
    working_dir=$(pwd)

    # Append "tree" to the data directory, create that directory, and enter it
    datadir="${datadir}tree"
    if [ ! -d "$datadir" ]; then
        mkdir "$datadir"
        cd "$datadir" || (echo "ERROR: Cannot cd into directory '$datadir'"; exit 5)
    else
        echo "ERROR: '$datadir' already exists"
        exit 4
    fi

    # Perform the crawl
    echo "Start date/time on $(hostname): $(date)"
    echo "Inputs:"
    echo "  base_url: $base_url"
    echo "  top_directories: $top_directories"
    echo "  datadir: $datadir"
    echo "Working directory: $working_dir"
    crawl "$base_url" "$top_directories" "$datadir"
    echo "End date/time on $(hostname): $(date)"

    # Create two other lists of the files/dirs in the tree we just created
    mv "$datadir/wget_err.txt" "$working_dir"
    tree > "$working_dir/all_files_and_dirs.txt"
    find . -type f > "$working_dir/all_files_in_tree.txt"

    # Go back into the original working directory
    cd "$working_dir" || (echo "ERROR: Cannot cd into directory '$working_dir'"; exit 6)

    # Print out the likely next command we'd want to run
    echo "Likely next command to run:"
    echo "  download_files \"$base_url\" \"$datadir\" \"$working_dir/all_files_in_tree.txt\""

}


# This function should not in general be called directly; it generates the list of datafiles from the file index
get_datafile_list() {
    file_index=$1
    grep "/L3/" "$file_index" | awk -v FS="____" '{print $1}' | grep -E ".*\.gene\.quantification\.txt$|.*\.gene\.fpkm\.txt$|chla\.org_NBL\.HumanExon\.Level-3\.BER\..*_gene\..*\.txt$|.*\.expression\.txt$" | grep -iv "/archive" | awk '{sub(/^\.\//,""); printf("%s\n",$0)}' | sort -u
}


# Main function to call directly in order to download the appropriate files from the TARGET site
download_files() {

    # Sample call:
    #   download_files "https://target-data.nci.nih.gov/" "/data/BIDS-HPC/private/projects/dmi/data/" "/home/weismanal/notebook/2020-04-08/scraping_target_site/all_files_in_tree.txt"

    # Parameters
    base_url=$1
    datadir=$2
    file_index=$3

    # Append "tree" to the data directory
    datadir="${datadir}tree/"

    # Get the list of files to download from the full list of indexed files
    datafiles=$(get_datafile_list "$file_index")
    #datafiles=$(grep "/L3/" "$file_index" | awk -v FS="____" '{print $1}' | grep -E ".*\.gene\.quantification\.txt$|.*\.gene\.fpkm\.txt$|chla\.org_NBL\.HumanExon\.Level-3\.BER\..*_gene\..*\.txt$|.*\.expression\.txt$" | grep -iv "/archive" | awk '{sub(/^\.\//,""); printf("%s\n",$0)}' | sort -u)

    # For each file to download...
    for datafile in $datafiles; do

        # Determine the local filename to download, the weblink where it's located, and the current placeholder file (commented out)
        filename="$datadir${datafile}"
        weblink="$base_url$datafile"
        #ls "${filename}____"* # this will get the placeholder file

        # If the downloaded file doesn't already exist, download it; otherwise, print a warning that it isn't being downloaded
        if [ ! -f "$filename" ]; then
            echo "Downloading file '$filename'"
            wget "$weblink" --output-document="$filename" 2>> "wget_err-download.txt"
        else
            echo "WARNING: File '$filename' already exists; skipping download"
        fi
    done

}


# Main function to call directly in order to extract the datafiles into TSV files to be subsequently read into Pandas dataframes in Python
extract_data() {

    # Sample call:
    #   extract_data "https://target-data.nci.nih.gov/" "/data/BIDS-HPC/private/projects/dmi/data/" "/home/weismanal/notebook/2020-04-08/scraping_target_site/all_files_in_tree.txt"

    # Parameters
    base_url=$1
    datadir=$2
    file_index=$3

    # Create a directory to hold the TSV files
    tsv_dir="${datadir}tsv_files/"
    if [ ! -d "$tsv_dir" ]; then
        mkdir "$tsv_dir"
    else
        echo "ERROR: TSV file directory '$tsv_dir' already exists"
        exit 8
    fi

    # Set the JSON filenames
    filedata_json="${datadir}filedata.json"
    metadata_json="${datadir}metadata.json"

    # Append "tree" to the data directory
    datadir="${datadir}tree/"

    # Get the list of files that should have been downloaded using the full list of indexed files
    datafiles=$(get_datafile_list "$file_index")
    #datafiles=$(grep "/L3/" "$file_index" | awk -v FS="____" '{print $1}' | grep -E ".*\.gene\.quantification\.txt$|.*\.gene\.fpkm\.txt$|chla\.org_NBL\.HumanExon\.Level-3\.BER\..*_gene\..*\.txt$|.*\.expression\.txt$" | grep -iv "/archive" | awk '{sub(/^\.\//,""); printf("%s\n",$0)}' | sort -u)
    ndatafiles=$(echo "$datafiles" | awk '{print NF}')

    # For each datafile...
    idatafile=0
    for datafile in $datafiles; do

        # Determine the local filename and the weblink where it's theoretically located
        filename="$datadir${datafile}"
        weblink="$base_url$datafile"

        # If the downloaded file doesn't exist, exit with an error; otherwise, extract its data
        if [ ! -f "$filename" ]; then
            echo "ERROR: Datafile '$filename' does not exist"
            exit 7
        else

            # Output the current file we're processing
            echo "Processing $filename..."

            # Define the name of the TSV file we're going to create
            tsv_file="${tsv_dir}tsv_file_$(printf "%07i" $idatafile)"
            


            # For now, assuming a uniform type of file, process the file accordingly, saving the result in the current TSV file (next up: create such a block for each file format!!)
            awk '{
                if (NR==1)
                    printf("%s\t%s\t%s\n", "gene-pretty", "gene-ugly", toupper($4))
                else {
                    split($1, arr, "|")
                    printf("%s\t%s\t%s\n", toupper(arr[1]), toupper(arr[2]), $4)
                }
            }' "$filename" > "$tsv_file"



            # Save the data corresponding to each file so we can later create a JSON file of all the file data
            if [ "x$ndatafiles" != "x1" ]; then
                if [ "x$idatafile" == "x0" ]; then
                    filenames="['$filename'"
                    weblinks="['$weblink'"
                    idatafiles="[$idatafile"
                    tsv_files="['$tsv_file'"
                elif [ "x$idatafile" == "x$((ndatafiles-1))" ]; then
                    filenames="$filenames, '$filename']"
                    weblinks="$weblinks, '$weblink']"
                    idatafiles="$idatafiles, $idatafile]"
                    tsv_files="$tsv_files, '$tsv_file']"
                else
                    filenames="$filenames, '$filename'"
                    weblinks="$weblinks, '$weblink'"
                    idatafiles="$idatafiles, $idatafile"
                    tsv_files="$tsv_files, '$tsv_file'"
                fi
            else
                filenames="['$filename']"
                weblinks="['$weblink']"
                idatafiles="[$idatafile]"
                tsv_files="['$tsv_file']"
            fi

        fi

        # Increate the datafile index
        idatafile=$((idatafile+1))
        
    done

    # Save the file data to a JSON file
    echo "
    filedata = {
        'filenames': $filenames,
        'weblinks': $weblinks,
        'idatafiles': $idatafiles,
        'tsv_files': $tsv_files
    }
    " > "$filedata_json"

    # Save the metadata to a JSON file
    echo "
    metadata = {
        'base_url': '$base_url',
        'datadir': '$datadir',
        'file_index': '$file_index',
        'working_dir': '$(pwd)',
        'ndatafiles': $ndatafiles
    }
    " > "$metadata_json"

}
