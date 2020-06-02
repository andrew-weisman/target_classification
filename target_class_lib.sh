#!/bin/bash

#
# See README.md (same directory) for instructions on using this library
#
# Example of project_dir: /data/BIDS-HPC/private/projects/dmi/
# Examples of working_dir:
#                           /home/weismanal/notebook/2020-04-08/scraping_target_site/
#                           /home/weismanal/notebook/2020-05-18/dmi/
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
    #   get_file_placeholders "https://target-data.nci.nih.gov/" "Controlled/ Public/" "${project_dir}data/"
    #   get_file_placeholders "https://target-data.nci.nih.gov/Controlled/AML/mRNA-seq/L3/structural/" "BCCA/ NCI-Meerzaman/" "${project_dir}data/"
    #   get_file_placeholders "https://target-data.nci.nih.gov/Public/RT/mRNA-seq/" "L3/ METADATA/" "${project_dir}data/"

    # Parameters
    base_url=$1
    top_directories=$2
    datadir=$3

    # Check for existence of all_files_in_tree.txt
    if [ -f "${datadir}all_files_in_tree.txt" ]; then
        echo "ERROR: File ${datadir}all_files_in_tree.txt already exists"
        exit
    fi

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
    #find . -type f > "$working_dir/all_files_in_tree.txt"
    find . -type f > "${datadir}/../all_files_in_tree.txt"

    # Go back into the original working directory
    cd "$working_dir" || (echo "ERROR: Cannot cd into directory '$working_dir'"; exit 6)

    # Print out the likely next command we'd want to run
    echo "Likely next command to run:"
    #echo "  download_files \"$base_url\" \"$datadir\" \"$working_dir/all_files_in_tree.txt\""
    echo "  download_files \"$base_url\" \"$datadir\" \"${datadir}/../all_files_in_tree.txt\""

}


# This function should not in general be called directly; it generates the list of datafiles from the file index
get_datafile_list() {
    file_index=$1
    grep "/L3/" "$file_index" | awk -v FS="____" '{print $1}' | grep -E ".*\.gene\.quantification\.txt$|.*\.gene\.fpkm\.txt$|chla\.org_NBL\.HumanExon\.Level-3\.BER\..*_gene\..*\.txt$|.*\.expression\.txt$" | grep -iv "/archive" | awk '{sub(/^\.\//,""); printf("%s\n",$0)}' | sort -u
}


# Main function to call directly in order to download the appropriate files from the TARGET site
download_files() {

    # Sample call:
    #   download_files "https://target-data.nci.nih.gov/" "${project_dir}data/"

    # Parameters
    base_url=$1
    datadir=$2

    # Constant
    file_index="${datadir}all_files_in_tree.txt"

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
    #   get_unique_headers "${project_dir}data/"
    #   THEN GO THROUGH PART (2) OF THIS FUNCTION!!

    # Parameters
    datadir=$1

    # Ensure we're not overwriting the output file
    if [ -f "$datadir/unique_headers.txt" ]; then
        echo "ERROR: File $datadir/unique_headers.txt already exists"
        exit
    fi

    # Constant
    file_index="${datadir}all_files_in_tree.txt"

    # (1) Print the unique headers to file
    for datafile in $(get_datafile_list "$file_index"); do
        filename="${datadir}tree/${datafile}"
        head -n 1 "$filename"
    done | sort -u > "$datadir/unique_headers.txt"

    # (2) MANUALLY delete (all but one of) the ones whose headers (in unique_headers.txt) include file-specific information, e.g., the patient ID, e.g., "probeset_id	TARGET-30-PAAPFA"
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
    #   get_format_mapping "${project_dir}data/"

    # Parameters
    datadir=$1

    # Ensure we don't overwrite the output file
    if [ -f "$datadir/format_mapping.txt" ]; then
        echo "ERROR: File $datadir/format_mapping.txt already exists"
        exit
    fi

    # Constant
    file_index="${datadir}all_files_in_tree.txt"

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
    #   ensure_format_consistency "${project_dir}data/"

    # Note: From the results here we see that format #4 results in non-uniform formats, that's why we also had to implement splitting based on the basename (and potentially later, on the full filepath)

    # Parameter
    datadir=$1

    # Check for existence of uniformity_check.txt
    if [ -f "${datadir}uniformity_check.txt" ]; then
        echo "ERROR: File ${datadir}uniformity_check.txt already exists"
        exit
    fi

    # Variables
    unique_headers="${datadir}unique_headers.txt"
    mapping_file="${datadir}format_mapping.txt"

    # Determine the number of unique file formats
    nformats=$(wc -l "$unique_headers" | awk '{print $1}')

    # For each unique file format...
    nfiles_holder=""
    {
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
    } > "${datadir}uniformity_check.txt"
    
}


# Extract the "best" gene name from a line in a file of each type
get_best_gene_name_from_file() {
    format_num=$1
    filename=$2
    case "$format_num" in
        1)
            awk -v format_num="$format_num" -v filename="$filename" '{
                if (NR!=1) {
                    printf("%s\t%s\t%s\n", toupper($1), format_num, filename)
                }
            }' "$filename"
            ;;
        2)
            awk -v format_num="$format_num" -v filename="$filename" '{
                if (NR!=1) {
                    printf("%s\t%s\t%s\n", toupper($1), format_num, filename)
                }
            }' "$filename"
            ;;
        3)
            awk -v format_num="$format_num" -v filename="$filename" '{
                if (NR!=1) {
                    printf("%s\t%s\t%s\n", toupper($1), format_num, filename)
                }
            }' "$filename"
            ;;
        4)
            awk -v format_num="$format_num" -v filename="$filename" '{
                if (NR!=1) {
                    split($1, arr, "|")
                    printf("%s\t%s\t%s\n", toupper(arr[2]), format_num, filename)
                }
            }' "$filename"
            ;;
        5)
            awk -v format_num="$format_num" -v filename="$filename" '{
                if (NR!=1) {
                    printf("%s\t%s\t%s\n", toupper($1), format_num, filename)
                }
            }' "$filename"
            ;;
        6)
            awk -v format_num="$format_num" -v filename="$filename" '{
                if (NR!=1) {
                    split($1, arr, ".")
                    printf("%s\t%s\t%s\n", toupper(arr[1]), format_num, filename)
                }
            }' "$filename"
            ;;
        7)
            awk -v format_num="$format_num" -v filename="$filename" '{
                if (NR!=1) {
                    printf("%s\t%s\t%s\n", toupper($1), format_num, filename)
                }
            }' "$filename"
            ;;
        8)
            awk -v format_num="$format_num" -v filename="$filename" '{
                if (NR!=1) {
                    printf("%s\t%s\t%s\n", toupper($1), format_num, filename)
                }
            }' "$filename"
            ;;
        *)
            echo "ERROR: Unknown format number ($format_num)"
            exit
            ;;
    esac
}


# Create a file called all_best_gene_names.txt in the data directory that writes out the "best" gene name in each line of each datafile
get_best_gene_names_from_all_files() {

    # Sample call:
    #   get_best_gene_names_from_all_files "${project_dir}data/"

    # Parameters
    datadir=$1

    # Ensure we're not overwriting an existing file
    if [ -f "${datadir}all_best_gene_names.txt" ]; then
        echo "ERROR: File ${datadir}all_best_gene_names.txt already exists"
        exit
    fi

    # Constant
    file_index="${datadir}all_files_in_tree.txt"

    # Get the list of files that should have been downloaded using the full list of indexed files
    datafiles=$(get_datafile_list "$file_index")

    # For each datafile...
    {
        for datafile in $datafiles; do

            # Determine the local filename
            filename="${datadir}tree/${datafile}"

            # If the downloaded file doesn't exist, exit with an error; otherwise, extract its data
            if [ ! -f "$filename" ]; then
                echo "ERROR: Datafile '$filename' does not exist"
                exit 7
            else

                # Get the format number for the current file
                format_num=$(determine_file_format "$filename" "$datadir/unique_headers.txt" "$datadir/unique_basenames.txt")
                
                # Get the best gene name from the current file
                get_best_gene_name_from_file "$format_num" "$filename"

            fi

        done
    } > "${datadir}all_best_gene_names.txt"

}


# Given a raw data file, process it using awk to make it a more uniform format prior to reading into Pandas
get_uniform_version_of_file() {
    format_num=$1
    filename=$2
    case "$format_num" in
        1) # done
            awk '{
                if (NR==1)
                    #printf("%s\t%s\t%s\t%s\t%s\n", "id", "tpm", "mean_length", "mean_eff_length", "est_counts")
                    printf("%s\t%s\t%s\t%s\t%s\n", "name", "tpm", "mean_length", "mean_eff_length", "est_counts")
                else {
                    printf("%s\t%s\t%s\t%s\t%s\n", toupper($1), $5, $2, $3, $4)
                }
            }' "$filename"
            ;;
        2) # done
            awk '{
                if (NR==1)
                    printf("%s\t%s\t%s\n", "name", "fpkm", "fpkm_normalized")
                else {
                    printf("%s\t%s\t%s\n", toupper($1), $2, $3)
                }
            }' "$filename"
            ;;
        3) # done
            awk '{
                if (NR==1)
                    printf("%s\t%s\t%s\n", "name", "rpkm", "raw_count")
                else {
                    printf("%s\t%s\t%s\n", toupper($1), $3, $2)
                }
            }' "$filename"
            ;;
        4) # done
            awk '{
                if (NR==1)
                    #printf("%s\t%s\t%s\t%s\t%s\n", "id", "rpkm", "name", "raw_counts", "median_length_normalized")
                    printf("%s\t%s\t%s\t%s\t%s\n", "name", "rpkm", "name2", "raw_counts", "median_length_normalized")
                else {
                    split($1, arr, "|")
                    printf("%s\t%s\t%s\t%s\t%s\n", toupper(arr[2]), $4, toupper(arr[1]), $2, $3)
                }
            }' "$filename"
            ;;
        5) # done
            awk '{
                if (NR==1)
                    printf("%s\t%s\t%s\t%s\n", "name", "rpkm", "raw_counts", "median_length_normalized")
                else {
                    printf("%s\t%s\t%s\t%s\n", toupper($1), $4, $2, $3)
                }
            }' "$filename"
            ;;
        6) # done
            awk '{
                if (NR==1)
                    printf("%s\t%s\t%s\t%s\t%s\n", "name", "tpm", "length", "effective_length", "num_reads")
                else {
                    split($1, arr, ".")
                    printf("%s\t%s\t%s\t%s\t%s\n", toupper(arr[1]), $4, $2, $3, $5)
                }
            }' "$filename"
            ;;
        7) # done
            awk '{
                if (NR==1)
                    printf("%s\t%s\n", "name", "UNKNOWN")
                else {
                    printf("%s\t%s\n", toupper($1), $2)
                }
            }' "$filename"
            ;;
        8) # done
            awk '{
                if (NR==1)
                    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "name", "fpkm", "name2", "name3", "class_code", "nearest_ref_id", "tss_id", "locus", "length", "coverage", "fpkm_conf_lo", "fpkm_conf_hi", "fpkm_status")
                else {
                    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", toupper($1), $10, toupper($4), toupper($5), $2, $3, $6, $7, $8, $9, $11, $12, $13)
                }
            }' "$filename"
            ;;
        *)
            echo "ERROR: Unknown format number ($format_num)"
            exit
            ;;
    esac
}


# Main function to call directly in order to extract the datafiles into TSV files to be subsequently read into Pandas dataframes in Python
# It also writes the project and file metadata
extract_data() {

    # Sample call:
    #   extract_data "https://target-data.nci.nih.gov/" "${project_dir}data/" |& tee "${working_dir}extract_data_out_and_err.txt"

    # Parameters
    base_url=$1
    datadir=$2

    # Constant
    file_index="${datadir}all_files_in_tree.txt"

    # Create a directory to hold the TSV files
    tsv_dir="${datadir}tsv_files/"
    if [ ! -d "$tsv_dir" ]; then
        mkdir "$tsv_dir"
    else
        echo "ERROR: TSV file directory '$tsv_dir' already exists"
        exit 8
    fi

    # Set the JSON filename
    metadata_json="${datadir}metadata.json"

    # Append "tree" to the data directory
    datadir="${datadir}tree/"

    # Get the list of files that should have been downloaded using the full list of indexed files
    datafiles=$(get_datafile_list "$file_index")
    #ndatafiles=$(echo "$datafiles" | awk '{print NF}')
    ndatafiles=$(echo $datafiles | awk '{print NF}') # not quoting this on purpose!

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
            tsv_file="${tsv_dir}tsv_file_$(printf "%07i" $idatafile).tsv"

            # Determine the format number of the current file
            format_num=$(determine_file_format "$filename" "${datadir}../unique_headers.txt" "${datadir}../unique_basenames.txt")
            
            # Given a raw data file, process it using awk to make it a more uniform format prior to reading into Pandas
            #echo $tsv_file
            get_uniform_version_of_file "$format_num" "$filename" > "$tsv_file"

            # Save the data corresponding to each file so we can later create a JSON file of all the file data
            if [ "x$ndatafiles" != "x1" ]; then
                if [ "x$idatafile" == "x0" ]; then
                    filenames="[\"$filename\""
                    weblinks="[\"$weblink\""
                    idatafiles="[$idatafile"
                    tsv_files="[\"$tsv_file\""
                    format_nums="[$format_num"
                elif [ "x$idatafile" == "x$((ndatafiles-1))" ]; then
                    filenames="$filenames, \"$filename\"]"
                    weblinks="$weblinks, \"$weblink\"]"
                    idatafiles="$idatafiles, $idatafile]"
                    tsv_files="$tsv_files, \"$tsv_file\"]"
                    format_nums="$format_nums, $format_num]"
                else
                    filenames="$filenames, \"$filename\""
                    weblinks="$weblinks, \"$weblink\""
                    idatafiles="$idatafiles, $idatafile"
                    tsv_files="$tsv_files, \"$tsv_file\""
                    format_nums="$format_nums, $format_num"
                fi
            else
                filenames="[\"$filename\"]"
                weblinks="[\"$weblink\"]"
                idatafiles="[$idatafile]"
                tsv_files="[\"$tsv_file\"]"
                format_nums="[$format_num]"
            fi

        fi

        # Increate the datafile index
        idatafile=$((idatafile+1))
        
    done

    # Save the metadata and file data to a JSON file
    echo "
    {
        \"metadata\": {
            \"base_url\": \"$base_url\",
            \"datadir\": \"$datadir\",
            \"file_index\": \"$file_index\",
            \"working_dir\": \"$(pwd)\",
            \"ndatafiles\": $ndatafiles
        },
        \"filedata\": {
            \"filenames\": $filenames,
            \"weblinks\": $weblinks,
            \"idatafiles\": $idatafiles,
            \"tsv_files\": $tsv_files,
            \"format_nums\": $format_nums
        }
    }
    " > "$metadata_json"

}







# # Create a file containing three columns for each unique best gene name in the entire dataset: (1) the number of matches of symbols from the lookup table, (2) the number of matches of the IDs from the lookup table, (3) the unique best gene name itself
# show_whether_genes_are_known() {
#     # Sample call:
#     #   show_whether_genes_are_known "${project_dir}data/" "$working_dir"

#     # Parameters
#     datadir=$1
#     working_dir=$2

#     # Ensure we're not overwriting an already created output file
#     if [ -f "${datadir}whether_genes_are_known.txt" ]; then
#         echo "ERROR: Output file ${datadir}whether_genes_are_known.txt already exists"
#         exit
#     fi

#     # Temporary files to go in the working directory
#     tmp_symbol_file="${working_dir}tmp_symbol_results.txt"
#     tmp_id_file="${working_dir}tmp_id_results.txt"

#     # For every unique best gene name in the entire dataset...
#     {
#         while read -r name; do

#             # Export any matches in the symbols file
#             grep "$name" "${datadir}lookup_symbol_uppercase.txt" > "$tmp_symbol_file"

#             # Export any matches in the IDs file
#             grep "$name" "${datadir}lookup_id_uppercase.txt" > "$tmp_id_file"

#             # Output the number of matches in each file type, as well as the best gene name itself
#             echo -e "$(wc -l "$tmp_symbol_file" | awk '{print $1}')\t$(wc -l "$tmp_id_file" | awk '{print $1}')\t$name"

#         done < "${datadir}unique_best_gene_names_uppercase.txt"
#     } > "${datadir}whether_genes_are_known.txt"
# }


# # For each of the four types of gene names identified in show_whether_genes_are_known() above (i.e., how they do in the HGNC lookup table), do as best as possible using the HGNC lookup table to get the Ensembl IDs, creating a full-size lookup table with blank third columns if really no Ensembl ID is currently known
# create_partial_global_lookup_table() {
#     # Sample call:;
#     #   create_partial_global_lookup_table "${project_dir}data/"

#     # Parameter
#     datadir=$1

#     # Don't overwrite an already existing output file
#     if [ -f "${datadir}partial_global_lookup_table.txt" ]; then
#         echo "ERROR: Output file ${datadir}partial_global_lookup_table.txt already exists"
#         exit
#     fi

#     {

#         # Output the trivial lookup for these two sets, whose gene names are already the Ensembl IDs
#         awk '($1==0 && $2==1){split($3,arr,"|"); printf("lookup:\t%s\t%s\n", arr[2], arr[2])}' "${datadir}whether_genes_are_known.txt" # these are gene names that match to a single Ensembl ID, so we know the Ensembl ID for them (N=37822)
#         awk '($1==0 && $2==2){split($3,arr,"|"); printf("lookup:\t%s\t%s\n", arr[2], arr[2])}' "${datadir}whether_genes_are_known.txt" # these are gene names that match to 2 Ensembl IDs, but they are the same Ensembl ID obviously, so we know the Ensembl ID; there is no ambiguity for us (though there is for the HGNC database) (N=3)

#         # For every gene name that matches to a single symbol, see if an Ensembl ID exists in the HGNC lookup table and print out the lookup; otherwise make a note that no Ensembl ID exists
#         while read -r name; do
#             tail -n +2 "${datadir}gene_lookup_table.txt" | awk -v name="$name" '
#             BEGIN { count = 0 }
#             {
#                 upper_2 = toupper($2)
#                 if (upper_2==name && $3!="") {
#                     count++
#                     #printf("MATCH FOUND (%i):\t%s\t%s\n", count, toupper($3), upper_2)
#                     printf("lookup:\t%s\t%s\n", name, toupper($3))
#                 }
#             }
#             END {
#                 if (count==0) printf("unmatched:\t%s\n", name)
#             }
#             '
#         done < <(awk '($1==1 && $2==0){split($3,arr,"|"); print arr[2]}' "${datadir}whether_genes_are_known.txt") # these are gene names that match to a single HGNC symbol (N=27624)

#         # Output every gene name in this set as unmatched since it wasn't found at all in the HGNC lookup table, unless it's an Ensembl ID (remember the HGNC lookup table doesn't include all Ensembl IDs)
#         awk '($1==0 && $2==0){
#             split($3,arr,"|")
#             name = arr[2]
#             if (name ~ /^ENSG[0-9]{11}$/)
#                 printf("lookup:\t%s\t%s\n", name, name)
#             else
#                 printf("unmatched:\t%s\n", name)
#         }' "${datadir}whether_genes_are_known.txt" # these are gene names that don't match anything in the lookup table (N=63161)

#     } > "${datadir}partial_global_lookup_table.txt"
# }


# lookup() {
#     name=$1
#     base_url="https://rest.ensembl.org/"
#     endpoint_lookup_symbol="lookup/symbol/homo_sapiens/"
#     endpoint_lookup_id="lookup/id/"
#     endpoint_xrefs_symbol="xrefs/symbol/homo_sapiens/"
#     endpoint_xrefs_id="xrefs/id/"
#     endpoint_xrefs_name="xrefs/name/homo_sapiens/"
#     endpoint_archive_id="archive/id/"
#     curl "${base_url}${endpoint_lookup_symbol}${name}" -H 'Content-type:application/json' -H 'Accept:application/json'
#     #curl "${base_url}${endpoint_lookup_id}${name}" -H 'Content-type:application/json' -H 'Accept:application/json'
#     #curl "${base_url}${endpoint_xrefs_symbol}${name}" -H 'Content-type:application/json' -H 'Accept:application/json' # this is the good one, but no POST!
#     #curl "${base_url}${endpoint_xrefs_id}${name}" -H 'Content-type:application/json' -H 'Accept:application/json'
#     #curl "${base_url}${endpoint_xrefs_name}${name}" -H 'Content-type:application/json' -H 'Accept:application/json'
#     #curl "${base_url}${endpoint_archive_id}${name}" -H 'Content-type:application/json' -H 'Accept:application/json'
# }


# testing() {
#     names_in="ENSG00000259762 ENSG00000259937 NAMPTL MINOS1P3 MINOS1P4 AC000029.1 AC007253.1 A3GALT2P A1BG"
#     names_out="ENSG00000259765 ENSG00000259936 NAG18 MINOS1P1 AC000032.2 AC000003.1 AC007251.2"

#     for name in $names_in; do
#         echo -e "\n$name"
#         lookup "$name"
#         echo ""
#     done
# }


# base_url="https://rest.ensembl.org/"
# endpoint_lookup_symbol="lookup/symbol/homo_sapiens/"
# endpoint_xrefs_symbol="xrefs/symbol/homo_sapiens/"
# #curl "${base_url}${endpoint_lookup_symbol}" -H 'Content-type:application/json' -H 'Accept:application/json' -X POST -d '{ "symbols" : ["ENSG00000259762", "ENSG00000259937", "NAMPTL", "MINOS1P3", "MINOS1P4", "AC000029.1", "AC007253.1", "A3GALT2P", "A1BG"] }'
# curl "${base_url}${endpoint_xrefs_symbol}" -H 'Content-type:application/json' -H 'Accept:application/json' -X POST -d '{ "symbols" : ["ENSG00000259762", "ENSG00000259937", "NAMPTL", "MINOS1P3", "MINOS1P4", "AC000029.1", "AC007253.1", "A3GALT2P", "A1BG"] }'


# testing() {

#     #filename="/data/BIDS-HPC/private/projects/dmi/data/tree/Public/OS/Toronto/mRNA-seq/L3/expression/NCI-Meltzer/TARGET-40-0A4HLQ-01A-01R.gene.quantification.txt" # 1
#     #filename="/data/BIDS-HPC/private/projects/dmi/data/tree/Public/ALL/mRNA-seq/Phase2/L3/expression/StJude/TARGET-10-PARASZ-09A-01R.expression.txt" # 2
#     #filename="/data/BIDS-HPC/private/projects/dmi/data/tree/Public/AML/mRNA-seq/L3/expression/NCI-Meerzaman/AML_23196.gene.quantification.txt" # 3
#     #filename="/data/BIDS-HPC/private/projects/dmi/data/tree/Public/ALL/mRNA-seq/Phase1/L3/expression/BCCA/HS0825.gene.quantification.txt" # 4
#     #filename="/data/BIDS-HPC/private/projects/dmi/data/tree/Public/ALL/mRNA-seq/Phase1/L3/expression/StJude/SJBALL010.gene.quantification.txt" # 5
#     #filename="/data/BIDS-HPC/private/projects/dmi/data/tree/Public/OS/mRNA-seq/L3/expression/NCI-Meltzer/TARGET-40-0A4HLD-01A-01R.gene.quantification.txt" # 6
#     #filename="/data/BIDS-HPC/private/projects/dmi/data/tree/Public/NBL/gene_expression_array/L3/gene/Core/chla.org_NBL.HumanExon.Level-3.BER.core_gene.TARGET-30-PAAPFA-01A-01R.txt" # 7
#     filename="/data/BIDS-HPC/private/projects/dmi/data/tree/Public/CCSK/mRNA-seq/L3/expression/NCI-Khan/CCSK002.gene.fpkm.txt" # 8
#     format_num=8

#     awk '{
#         if (NR==1)
#             printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "name", "fpkm", "name2", "name3", "class_code", "nearest_ref_id", "tss_id", "locus", "length", "coverage", "fpkm_conf_lo", "fpkm_conf_hi", "fpkm_status")
#         else {
#             printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", toupper($1), $10, toupper($4), toupper($5), $2, $3, $6, $7, $8, $9, $11, $12, $13)
#         }
#     }' "$filename"

#     # awk -v format_num="$format_num" -v filename="$filename" '{
#     #     if (NR!=1) {
#     #         printf("%s\t%s\t%s\n", toupper($1), format_num, filename)
#     #     }
#     # }' "$filename"

# }


# Process a GTF file using awk for read-in using Python into a Pandas dataframe
process_gtf_file() {
    full_path=$1
    awk '
        BEGIN {
            printf("%s\t%s\t%s\t%s\n", "id", "version", "name", "biotype")
        }
        $3=="gene" {
            split($10, id, "\"")
            split($12, version, "\"")
            split($14, name, "\"")
            split($18, biotype, "\"")
            printf("%s\t%s\t%s\t%s\n", id[2], version[2], name[2], biotype[2])
        }
    ' "$full_path"
}


# Create processed TSV files for use by Pandas to generate a lookup table
make_tsv_files_for_lookup_table() {
    # Sample call:
    #   make_tsv_files_for_lookup_table /data/BIDS-HPC/private/projects/dmi/

    # Input (needs trailing forward slash at end of project directory path)
    project_dir=$1

    # Create the directory to contain the processed files for reading into Pandas, if it doesn't already exist
    if [ ! -d "${project_dir}data/processed_files_for_lookup" ]; then
        mkdir "${project_dir}data/processed_files_for_lookup"
    fi

    #### Process the GTF files if they don't already exist

    # This is the most up-to-date version of the human genome from Ensembl
    # ftp://ftp.ensembl.org/pub//release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
    if [ ! -f "${project_dir}data/processed_files_for_lookup/main_38.tsv" ]; then
        process_gtf_file "${project_dir}data/ensembl_ftp_site/grch38/Homo_sapiens.GRCh38.100.gtf" > "${project_dir}data/processed_files_for_lookup/main_38.tsv"
    fi

    # This is the most recent update to the prior version of the human genome from Ensembl
    # ftp://ftp.ensembl.org/pub//grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz (the most recent version actually links to release 87)
    if [ ! -f "${project_dir}data/processed_files_for_lookup/main_37.tsv" ]; then
        process_gtf_file "${project_dir}data/ensembl_ftp_site/grch37/Homo_sapiens.GRCh37.87.gtf" > "${project_dir}data/processed_files_for_lookup/main_37.tsv"
    fi

    # This is the genome release (GRCh38.p2) corresponding to Gencode release #22 (https://www.gencodegenes.org/human/release_22.html), which is what the GDC Data Portal uses as its data set, e.g. https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
    # ftp://ftp.ensembl.org/pub//release-80/gtf/homo_sapiens/Homo_sapiens.GRCh38.80.gtf.gz
    if [ ! -f "${project_dir}data/processed_files_for_lookup/main_38-80.tsv" ]; then
        process_gtf_file "${project_dir}data/ensembl_ftp_site/grch38/release-80/Homo_sapiens.GRCh38.80.gtf" > "${project_dir}data/processed_files_for_lookup/main_38-80.tsv"
    fi

    # This is the ena synonyms file for the most recent version
    # ftp://ftp.ensembl.org/pub/release-100/tsv/homo_sapiens/Homo_sapiens.GRCh38.100.ena.tsv.gz
    if [ ! -f "${project_dir}data/processed_files_for_lookup/xref_38-ena.tsv" ]; then
        {
            echo -e "id\tname"
            awk -v FS="\t" '{printf("%s\t%s\n", $3, $7)}' "${project_dir}data/ensembl_ftp_site/grch38/Homo_sapiens.GRCh38.100.ena.tsv" | awk 'NF==2{print}' | tail -n +2 | sort -u
        } > "${project_dir}data/processed_files_for_lookup/xref_38-ena.tsv"
    fi

    # This is the entrez synonyms file for the most recent version
    # ftp://ftp.ensembl.org/pub/release-100/tsv/homo_sapiens/Homo_sapiens.GRCh38.100.entrez.tsv.gz
    if [ ! -f "${project_dir}data/processed_files_for_lookup/xref_38-entrez.tsv" ]; then
        {
            echo -e "id\tname"
            awk -v FS="\t" '{printf("%s\t%s\n", $1, $4)}' "${project_dir}data/ensembl_ftp_site/grch38/Homo_sapiens.GRCh38.100.entrez.tsv" | tail -n +2 | sort -u
        } > "${project_dir}data/processed_files_for_lookup/xref_38-entrez.tsv"
    fi

    # This is the refseq synonyms file for the most recent version
    # ftp://ftp.ensembl.org/pub/release-100/tsv/homo_sapiens/Homo_sapiens.GRCh38.100.refseq.tsv.gz
    if [ ! -f "${project_dir}data/processed_files_for_lookup/xref_38-refseq.tsv" ]; then
        {
            echo -e "id\tname"
            awk -v FS="\t" '{printf("%s\t%s\n", $1, $4)}' "${project_dir}data/ensembl_ftp_site/grch38/Homo_sapiens.GRCh38.100.refseq.tsv" | tail -n +2 | sort -u
        } > "${project_dir}data/processed_files_for_lookup/xref_38-refseq.tsv"
    fi

    # This is the uniprot synonyms file for the most recent version
    # ftp://ftp.ensembl.org/pub/release-100/tsv/homo_sapiens/Homo_sapiens.GRCh38.100.uniprot.tsv.gz
    if [ ! -f "${project_dir}data/processed_files_for_lookup/xref_38-uniprot.tsv" ]; then
        {
            echo -e "id\tname"
            awk -v FS="\t" '{printf("%s\t%s\n", $1, $4)}' "${project_dir}data/ensembl_ftp_site/grch38/Homo_sapiens.GRCh38.100.uniprot.tsv" | tail -n +2 | sort -u
        } > "${project_dir}data/processed_files_for_lookup/xref_38-uniprot.tsv"
    fi

    # This is the uniprot synonyms file for the most recent version of GRCh37
    # ftp://ftp.ensembl.org/pub/grch37/release-85/tsv/homo_sapiens/Homo_sapiens.GRCh37.85.uniprot.tsv.gz
    if [ ! -f "${project_dir}data/processed_files_for_lookup/xref_37-uniprot.tsv" ]; then
        {
            echo -e "id\tname"
            awk -v FS="\t" '{printf("%s\t%s\n", $1, $4)}' "${project_dir}data/ensembl_ftp_site/grch37/Homo_sapiens.GRCh37.85.uniprot.tsv" | tail -n +2 | sort -u
        } > "${project_dir}data/processed_files_for_lookup/xref_37-uniprot.tsv"
    fi

    # This is the ena synonyms file for the most recent version of GRCh37
    # ftp://ftp.ensembl.org/pub/grch37/release-85/tsv/homo_sapiens/Homo_sapiens.GRCh37.85.ena.tsv.gz
    if [ ! -f "${project_dir}data/processed_files_for_lookup/xref_37-ena.tsv" ]; then
        {
            echo -e "id\tname"
            awk -v FS="\t" '{printf("%s\t%s\n", $3, $7)}' "${project_dir}data/ensembl_ftp_site/grch37/Homo_sapiens.GRCh37.85.ena.tsv" | awk 'NF==2{print}' | tail -n +2 | sort -u
        } > "${project_dir}data/processed_files_for_lookup/xref_37-ena.tsv"
    fi

}


process_all_gtf_files() {
    datadir="/home/weismanal/notebook/2020-06-01/dmi"
    processed_dir="/home/weismanal/notebook/2020-06-01/dmi/processed"
    for gtf_file in "$datadir"/*.gtf; do
        version=$(basename "$gtf_file" | awk '{split($1,arr,"."); split(arr[2],arr2,"GRCh"); print(arr2[2])}')
        release=$(basename "$gtf_file" | awk '{split($1,arr,"."); print(arr[3])}')
        processed_file="$processed_dir/main_$version-$release.tsv"
        process_gtf_file "$gtf_file" > "$processed_file"
    done
}


write_python_lines() {
    processed_dir="/home/weismanal/notebook/2020-06-01/dmi/processed"
    for tsv_file in "$processed_dir"/*.tsv; do
        echo "lookup = tc.add_tsv_file_to_lookup(lookup, '${tsv_file}')"
    done
}
