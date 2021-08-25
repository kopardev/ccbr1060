#!/bin/bash
# get detailed counts per sample
for i in {1..8};do echo $i;python create_detailed_window_count_matrix_per_sample.py test${i}/test${i}_p1.Chimeric.out.junction > test${i}.detailed_counts.txt;done
# get aggregate detailed counts table
echo -ne "sample\tintegration_window\t" > detailed_counts.txt
for i in {1..8};do while read a;do echo -ne "test${i}\t$a\n";done < test${i}.detailed_counts.txt;done|head -n1 |cut -f2- >> detailed_counts.txt
for i in {1..8};do while read a;do echo -ne "test${i}\t$a\n";done < test${i}.detailed_counts.txt;done|grep -v shRNA|sort -k2,2 -k1,1 >> detailed_counts.txt
