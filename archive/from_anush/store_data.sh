#!/bin/bash

#program that iterates over the patients adding to array and then calling obtain_array.py

old=" "
p=0

#look into all .txt files 
for file in $(find /home/kaovasiaar/check_outs/L3_data-links/Public -iname "*.txt" | sort); do
    #grab just the name of patient 
    new=$( basename $file | awk -v FS="." '{print $1}')
    #index if new patient (there are same patient names with different data - genes, exon, etc.) 
    if [ $new != $old ]; then 
	p=$(( p+1 ))

    fi
    #export the file and patient counter
    export file 
    export p
    old=$new

done
#run the python code
python obtain_array.py



