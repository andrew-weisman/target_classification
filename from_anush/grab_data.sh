#!/bin/bash
#program grabs specific files listed in target-gene document from data directory and
#puts in personal directory under L3_data_links

#delete first two lines of file since useless for extracting
tail -n +2 /data/kaovasiaar/target-gene-expr-level3.txt > /data/kaovasiaar/target-gene-expr-level3-fixed.txt
rm -rf /data/kaovasiaar/target-gene-expr-level3.txt

for file in $(awk '{print $2}' /data/kaovasiaar/target-gene-expr-level3-fixed.txt); do
    full_path="/data/kaovasiaar/L3_data/$file"
    extra_dir=$(ls $full_path &> /dev/null || ls $(dirname $full_path))
    if [ "a$extra_dir" != "a" ]; then
	mydir=$(dirname $full_path)
	myfile=$(basename $full_path)
	full_path="$mydir/$extra_dir/$myfile"
    fi
    #ignore this part of path and use only path after Public/ including file name
    suffix=$(echo $full_path | awk -v FS="/data/kaovasiaar/L3_data/" '{print $2}')
    #new directory name containing directories after Public including new one missed
    new_dir=$(dirname $suffix)
    #puts it in own path
    mkdir -p /home/kaovasiaar/check_outs/L3_data-links/$new_dir
    #puts it in the path made
    new_path=/home/kaovasiaar/check_outs/L3_data-links/$suffix
    if [ ! -f $new_path ]; then
	ln -s $full_path $new_path
    fi
     
done
