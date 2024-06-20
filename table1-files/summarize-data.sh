#!/bin/bash
set -e
enable -n echo

touch intro.data
rm intro.data

for a in fluA fluB fluC fluD; do 
    for b in 2019 2020 2021 2022 2023; do
        echo -n "$a $b gb " >> intro.data
        wc -l $a.$b.gblist | awk '{ print $1 }' >> intro.data
    done
    for b in all.20240220; do
        echo -n "$a $b gb " >> intro.data
        wc -l $a.$b.gblist | awk '{ print $1 }' >> intro.data
        echo -n "$a $b all " >> intro.data
        wc -l $a.$b.list | awk '{ print $1 }' >> intro.data
    done
done
