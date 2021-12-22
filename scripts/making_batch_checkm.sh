#!/bin/bash

dir_size=200
dir_name="0"
n=$((`find . -maxdepth 1 -type f | wc -l`/$dir_size+1))
for i in `seq 1 $n`;
do
    mkdir -p "$dir_name$i";
    find . -maxdepth 1 -type f | head -n $dir_size | xargs -i mv "{}" "$dir_name$i"
done
