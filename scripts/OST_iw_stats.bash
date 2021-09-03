#!/bin/bash
awk -F"\t" '{if ($3==0){print $1}}' integration_windows.nearest_lab_verified_lookup.txt |awk -F"|" '{print $1}' > integration_windows.lab-verified.txt
cat *aggregate_OST_calls.tsv|grep -v nreads|cut -f2|sort|uniq -c|awk '{print $2}'|awk -F"|" '{print $1}' > integration_windows.OST.txt
cat integration_windows.lab-verified.txt integration_windows.OST.txt |sort|uniq -c|awk  '{if ($1==2) {print $2}}' > integration_windows.OST.lab-verified.txt
wc -l integration_windows.annotated.bed
wc -l integration_windows.lab-verified.txt
wc -l integration_windows.OST.txt
wc -l integration_windows.OST.lab-verified.txt
