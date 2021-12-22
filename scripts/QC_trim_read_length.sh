#!/bin/bash

zcat $1| awk -v f=$3 '{if(NR%4==2) print length($1),f}' | sort -n | uniq -c > $2
#zcat $3| awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $4
