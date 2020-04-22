#!/bin/bash

get_dir() {
    local dir=$1
    echo $dir
    subdirs=$(curl https://target-data.nci.nih.gov/$dir 2> /dev/null | awk -v RS=" " '{print}' | grep "^href=" | awk '{split($0,arr1,"href=\""); split(arr1[2],arr2,"\""); print arr2[1]}' | egrep -v "^/|^\?|^http" | sort -u | grep "/$")
    for subdir in $subdirs; do
        #echo $dir$subdir | egrep "/VisCap_Female_Somatic/$|/VisCap_Male_Somatic/$" > /dev/null
        #if [ ! $? -eq 0 ]; then
        if [ $? -eq 0 ]; then
            get_dir $dir$subdir
        else
            echo $dir$subdir
        fi
    done
}

get_dir Public/