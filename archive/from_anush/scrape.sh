#!/bin/bash
#program scrapes website and prints all the names of directories and text files 

get_dir() {
    echo $1
    subdirs=$(curl https://target-data.nci.nih.gov/$1 2> /dev/null | awk -v RS=" " '{print}' | grep "^href=" | awk '{split($0,arr1,"href=\""); split(arr1[2],arr2,"\""); print arr2[1]}' | grep -v "^/\|^?\|^http" | sort -u | grep "/$")
    for subdir in $subdirs; do
        echo $1$subdir | grep "/VisCap_Female_Somatic/$\|/VisCap_Male_Somatic/$" > /dev/null
        if  [ ! $? -eq 0 ]; then
            get_dir $1$subdir
	    text_files=$(curl https://target-data.nci.nih.gov/$1$subdir 2> /dev/null | awk -v RS=" " '{print}' | grep "^href=" | awk '{split($0,arr1,"href=\""); split(arr1[2],arr2,"\""); print arr2[1]}' | grep -v "^/\|^?\|^http" | sort -u | grep "\.txt$")
	    for file in $text_files; do
		echo "$1$subdir$file"
	    done
	else
           echo $1$subdir
        fi
    done
}

get_dir Public/
