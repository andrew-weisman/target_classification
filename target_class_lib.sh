#!/bin/bash

#
# See README.md (same directory) for instructions on using this library
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

    # For each file to download...
    for datafile in $(get_datafile_list "$file_index"); do

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


# Function used for printing the unique headers in the datafiles; this determines the formats of the datafiles
get_unique_headers() {

    # Sample call:
    #   get_unique_headers "/data/BIDS-HPC/private/projects/dmi/data/" "/home/weismanal/notebook/2020-04-08/scraping_target_site/all_files_in_tree.txt"
    #   THEN GO THROUGH PART (2) OF THIS FUNCTION!!

    # Parameters
    datadir=$1
    file_index=$2

    # (1) Print the unique headers to file
    for datafile in $(get_datafile_list "$file_index"); do
        filename="${datadir}tree/${datafile}"
        head -n 1 "$filename"
    done | sort -u > "$datadir/unique_headers.txt"

    # (2) MANUALLY delete (all but one of) the ones whose headers include file-specific information, e.g., the patient ID, e.g., "probeset_id	TARGET-30-PAAPFA"
    # Then change those headers to a general regexp, e.g., "probeset_id	TARGET-30-PAAPFA" --> "probeset_id	TARGET-.+"
    # Finally, add a "^" to the beginning of each line in unique_headers.txt and add a "$" to the end of each line

    # (3) Since from the ensure_format_consistency() function we know that format #4 is not consistent, duplicate this row.
    # Then, since we know that defining the filename will make the formats unique, create another file called unique_basenames.txt with the same number of rows as unique_headers, and in the two rows corresponding to the non-distinguishing header formats, distinguish the filenames using regular expressions.

}


# Convert a Bash string to a Python list
bash_str_to_python_list() {
    x=$1
    echo "[$(echo $x | awk '{gsub(" ",","); print}')]" # purposefully not quoting $x here in case it contains newlines!
}


# Determine the file format number given a filename and the files containing the unique headers and unique basenames
determine_file_format() {

    # Parameters
    filename=$1
    unique_headers=$2
    unique_basenames=$3

    # Determine the 1 or more possible headers this file fits and store the possible number of headers
    header_ids=$(awk -v header="$(head -n 1 "$filename")" 'header~$0{print NR}' "$unique_headers")
    nheader_ids=$(echo "$header_ids" | wc -l)

    # If the number of headers is 1 then we have determined the format so we can just return what we already know
    if [ "$nheader_ids" -eq 1 ]; then
        echo "$header_ids"

    # Otherwise, we haven't fully narrowed down the format, so we need to use the filename (basename) as well
    else

        # Determine the 1 or more possible basenames this file fits
        basename_ids=$(awk -v bn="$(basename "$filename")" 'bn~$0{print NR}' "$unique_basenames")

        # Get the Python lists of the possible formats as determined by the headers and the basenames
        header_list=$(bash_str_to_python_list "$header_ids")
        basename_list=$(bash_str_to_python_list "$basename_ids")

        # Return the intersection of these two lists, which should result in a single number (if the steps were appropriately followed in get_unique_headers())
        python -c "print(list(set($header_list) & set($basename_list))[0])"

    fi

}


# Run this function to generate a master file format list
get_format_mapping() {

    # Sample call:
    #   get_format_mapping "/data/BIDS-HPC/private/projects/dmi/data/" "/home/weismanal/notebook/2020-04-08/scraping_target_site/all_files_in_tree.txt"

    # Parameters
    datadir=$1
    file_index=$2

    # Print a file containing the format number next to the filename
    for datafile in $(get_datafile_list "$file_index"); do
        filename="${datadir}tree/${datafile}"
        format_num=$(determine_file_format "$filename" "$datadir/unique_headers.txt" "$datadir/unique_basenames.txt")
        echo "$format_num" "$filename"
    done > "$datadir/format_mapping.txt"

}


# Print out a report of features of the files for each unique file format
ensure_format_consistency() {

    # Sample call:
    #   ensure_format_consistency "/data/BIDS-HPC/private/projects/dmi/data/" > uniformity_check.txt

    # Note: From the results here we see that format #4 results in non-uniform formats, that's why we also had to implement splitting based on the basename (and potentially later, on the full filepath)

    # Parameter
    datadir=$1

    # Variables
    unique_headers="${datadir}unique_headers.txt"
    mapping_file="${datadir}format_mapping.txt"

    # Determine the number of unique file formats
    nformats=$(wc -l "$unique_headers" | awk '{print $1}')

    # For each unique file format...
    nfiles_holder=""
    for format_num in $(seq 1 "$nformats"); do

        # Get the list of filenames of the current format, their corresponding basenames, and the number of files of this format
        filenames=$(awk -v format_num="$format_num" '$1==format_num{print $2}' "$mapping_file")
        basenames=$(echo "$filenames" | awk '{len=split($1,arr,"/"); print arr[len]}' | sort -u)
        nfiles=$(echo "$basenames" | wc -l)
        nfiles_holder="$nfiles_holder,$nfiles"
 
        # Output the current file format being analyzed
        echo -e "\n---- On format number $format_num ----\n"

        # Output the sorted list of basenames in order to observe their patterns
        echo -e "Basenames:\n"
        echo -e "$basenames\n"

        # Output the header line describing the current format
        echo -e "Header line:\n"
        awk -v format_num="$format_num" 'NR==format_num{print}' "$unique_headers"

        # Output the first non-header line of each file of the current format (along with the filename in order to help narrow down differences)
        echo -e "\nFirst non-header line of each file:\n"
        for filename in $filenames; do
            second_line=$(head -n 2 "$filename" | tail -n 1)
            echo -e "$second_line\t$filename"
        done

        # Output the number of files of the current format
        echo -e "\nNumber of files: $nfiles\n"

        # Confirm a single number of unique fields in the datafiles of the current format
        echo -e "Unique numbers of fields:\n"
        for filename in $filenames; do
            awk '{print NF}' "$filename" | sort -u
        done | sort -u
 
    done

    # Output the number of files in each format and ensure their sum equals, currently, 2261
    echo -e "\n\n-------------------"
    echo -e "Overall quantities:\n"
    python -c "
nfiles_holder=[${nfiles_holder:1:${#nfiles_holder}}]
print('Number of files in each format: {}'.format(nfiles_holder))
print('Total number of files: {}'.format(sum(nfiles_holder)))
    "
    
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



            #### Don't forget to save this to the JSON file as well!!!!
            format_num=$(determine_file_format "$filename" "$datadir/unique_headers.txt" "$datadir/unique_basenames.txt")
            


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
