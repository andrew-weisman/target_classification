#!/bin/bash
#program that goes through all the files with the L3 directory and creates directories
#with the linked text files in the /data/kaovasiaar/L3_data/ directory

#create temp_file that containes directory and file names
sh /home/kaovasiaar/check_outs/dmi/scrape.sh > temp_file

#doesn't create the multiple VisCap directories due to it being so many, can edit
#scrape.sh script later to include these
#creates directories from website with L3 data  
grep "/L3/" temp_file | grep -v "\.txt$"|while read n;do mkdir -p /data/kaovasiaar/L3_data/"$n";done

#grabs text files from website and places them in previously created directories
grep "/L3/.*\.txt$" temp_file | while read n; do
    cd /data/kaovasiaar/L3_data/$(dirname $n)
    wget "https://target-data.nci.nih.gov/$n"
    
done


