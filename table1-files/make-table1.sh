#!/bin/bash
set -e

touch table1.data
rm table1.data

for a in fluA fluB fluC fluD; do 
    for b in 2018 2019 2020 2021 2022 2023; do
        echo -n "$a $b gb " >> table1.data
        wc -l $a.$b.gblist | awk '{ print $1 }' >> table1.data
    done
    for b in all.20240220; do
        echo -n "$a $b gb " >> table1.data
        wc -l $a.$b.gblist | awk '{ print $1 }' >> table1.data
        echo -n "$a $b all " >> table1.data
        wc -l $a.$b.list | awk '{ print $1 }' >> table1.data
    done
done
perl table-nseqs-data2tex.pl table1.data > table1.tex
