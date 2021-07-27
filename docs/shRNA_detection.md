## Finding consensus shRNA integration sites across all samples

`/data/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/bams_withoutChimericJunctions` has the results for running STAR with chimeric junctions reported for each sample as a tab-delimited file.

```bash
(base) kopardevn@biowulf 11:00 PM /data/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/bams_withoutChimericJunctions
 HEADNODE! % find . -name "*_p1.Chimeric.out.junction"
./test5/test5_p1.Chimeric.out.junction
./test4/test4_p1.Chimeric.out.junction
./test2/test2_p1.Chimeric.out.junction
./test3/test3_p1.Chimeric.out.junction
./test1/test1_p1.Chimeric.out.junction
./test8/test8_p1.Chimeric.out.junction
./test6/test6_p1.Chimeric.out.junction
./test7/test7_p1.Chimeric.out.junction
```

Now I filter each of these files to only look at chimeric junctions involving shRNA, i.e., either the donor or the accepter of the junction must be shRNA the report the results in `filtered.chimeric.junctions`. These junctions have also been filtered to show only those which have greater than 5 reads of support.

```bash
 % sort -k1,1nr filtered.chimeric.junctions|head
2089	chr5	73642754	test7
2087	chr5	73642754	test5
2018	chr5	73642754	test6
1908	chr5	73642754	test8
585	chr5	73643317	test5
506	chr5	73643317	test8
484	chr5	73643317	test7
475	chr5	73643317	test6
449	chr7	44881390	test5
398	chr7	44881390	test6
 % sort -k1,1nr filtered.chimeric.junctions|tail
6	chr9	9442130	test6
6	chr9	9442137	test7
6	chr9	9442158	test8
6	chr9	9442170	test7
6	chr9	9442234	test6
6	chr9	9442235	test6
6	chr9	9442241	test7
6	chr9	9442241	test8
6	chr9	9442242	test5
6	chrM	2046	test7
```

The columns here are:

1. Number of read supporting the junction per sample
2. Host chromosome
3. Host coordinate
4. Sample name

Next, we aggregate junctions accross multiple samples in `filtered.chimeric.junctions.aggregate_across_samples`. If a host chromosome + host coordinate combination is seen in multiple samples, then the counts are summed up and a single row is reported.

```bash
 % sort -k3,3nr filtered.chimeric.junctions.aggregate_across_samples|head
chr5	73642754	8191
chr5	73643317	2152
chr7	44881390	1870
chr6	33272627	1335
chr1	189437490	752
chr7	44881427	670
chr12	124914855	631
chr13	95133445	613
chr2	11332341	551
chr5	113546953	492
 % sort -k3,3nr filtered.chimeric.junctions.aggregate_across_samples|tail
chr6	52995836	6
chr6	52995864	6
chr6	52995871	6
chr6	52995891	6
chr9	73160812	6
chr9	9442121	6
chr9	9442130	6
chr9	9442158	6
chr9	9442170	6
chrM	2046	6
```

The columns here are:

1. Host chromosome
2. Host coordinate
3. Total number of reads supporting this junctions, all samples combined

`create_integration_windows.py` python script then takes the above file and determines plausible shRNA integration windows.

* Integration host sites within 100bp of each other are grouped together as a _site-list_
* The coordinate of the _site_ with the highest read support within a _site-list_ is used to represent the coordinate of the entire _site-list_
* 100bp buffer is applied to either side of this representative coordinate to create the ` integration_windows.bed` 

Overlapping windows are simply merged together with bedtools.

This leads to a total of 54 putative shRNA integration sites:

```bash
 % wc -l integration_windows.bed
54 integration_windows.bed
```

These are then annotated with `genes.bed6` for hg38 (alignment is done against hg38 by PJ)

```bash
 % head integration_windows.bed.annotated.txt
chr1	632158	632357	chr1	631074	632616	ENSG00000237973.1|MTCO1P12	.	+
chr1	28507280	28507479	chr1	28507366	28507571	ENSG00000274266.1|SNORA73A	.	+
chr1	71072449	71072648	chr1	71063291	71081289	ENSG00000132485.14|ZRANB2	.	-
chr1	71076470	71076669	chr1	71063291	71081289	ENSG00000132485.14|ZRANB2	.	-
chr1	71076694	71076893	chr1	71063291	71081289	ENSG00000132485.14|ZRANB2	.	-
chr1	189437275	189437589	.	-1	-1	.	-1	.
chr1	212615110	212615309	chr1	212565334	212620777	ENSG00000162772.17|ATF3	.	+
chr1	212620006	212620205	chr1	212565334	212620777	ENSG00000162772.17|ATF3	.	+
chr11	65504957	65505156	chr11	65497688	65506516	ENSG00000251562.8|MALAT1	.	+
chr12	57517324	57517523	chr12	57516588	57521737	ENSG00000175197.13|DDIT3	.	-
```

3 of these are in intergenic regions... others are withing gene boundaries.

```bash
 % awk -F"\t" '{if ($7=="."){print}}' integration_windows.bed.annotated.txt
chr1	189437275	189437589	.	-1	-1	.	-1	.
chr21	8258274	8258649	.	-1	-1	.	-1	.
chr5	113512841	113513040	.	-1	-1	.	-1	.
```

The genic 51 integration sites fall in 35 genes

```bash
 % cut -f7 integration_windows.bed.annotated.txt | sort | uniq |grep ENS
ENSG00000047188.16|YTHDC2
ENSG00000067141.17|NEO1
ENSG00000106511.6|MEOX2
ENSG00000124209.4|RAB22A
ENSG00000125257.16|ABCC4
ENSG00000132485.14|ZRANB2
ENSG00000134318.14|ROCK2
ENSG00000135046.14|ANXA1
ENSG00000135220.11|UGT2A3
ENSG00000137154.13|RPS6
ENSG00000141279.17|NPEPPS
ENSG00000146676.10|PURB
ENSG00000150991.15|UBC
ENSG00000158019.21|BABAM2
ENSG00000162772.17|ATF3
ENSG00000175197.13|DDIT3
ENSG00000200488.1|RN7SKP203
ENSG00000202198.1|7SK
ENSG00000210082.2|MT-RNR2
ENSG00000214944.10|ARHGEF28
ENSG00000223501.9|VPS52
ENSG00000230897.1|RPS18P12
ENSG00000231500.7|RPS18
ENSG00000237973.1|MTCO1P12
ENSG00000240877.3|RN7SL521P
ENSG00000251562.8|MALAT1
ENSG00000259001.3|RPPH1
ENSG00000261499.2|CH17-260O16.1
ENSG00000263740.2|RN7SL4P
ENSG00000265735.2|RN7SL5P
ENSG00000274266.1|SNORA73A
ENSG00000276168.1|RN7SL1
ENSG00000278771.1|RN7SL3
ENSG00000280441.3|CH507-528H12.1
ENSG00000282885.2|RP11-596C23.6
```

There are 15 lab-verified integration sites.

```bash
 % more lab_verified_sites.bed
chr1	71542253	71542253	.	.	+
chr1	189406937	189406937	.	.	-
chr12	132507964	132507964	.	.	-
chr13	95784002	95784002	.	.	+
chr14	84280651	84280651	.	.	+
chr15	73354871	73354871	.	.	-
chr17	36364840	36364840	.	.	+
chr2	11472501	11472501	.	.	+
chr2	28211890	28211890	.	.	-
chr4	69801325	69801325	.	.	-
chr5	72939142	72939142	.	.	-
chr5	112903625	112903625	.	.	-
chr6	33223194	33223194	.	.	+
chr7	15721280	15721280	.	.	-
chr7	44920734	44920734	.	.	+
 % more lab_verified_sites.annotated.bed
chr1	71076570	71076570	chr1:71076570-71076570|+|ENSG00000132485.14|ZRANB2	.	+
chr1	189437807	189437807	chr1:189437807-189437807|-|.	.	-
chr12	132023419	132023419	chr12:132023419-132023419|-|ENSG00000183495.14|EP400	.	-
chr13	95131748	95131748	chr13:95131748-95131748|+|ENSG00000125257.16|ABCC4	.	+
chr14	83814307	83814307	chr14:83814307-83814307|+|.	.	+
chr15	73062530	73062530	chr15:73062530-73062530|-|ENSG00000067141.17|NEO1	.	-
chr17	38208780	38208780	chr17:38208780-38208780|+|ENSG00000274487.2|NPEPPSP1	.	+
chr2	11332375	11332375	chr2:11332375-11332375|+|ENSG00000134318.14|ROCK2	.	+
chr2	27989023	27989023	chr2:27989023-27989023|-|ENSG00000158019.21|BABAM2	.	-
chr4	68935607	68935607	chr4:68935607-68935607|-|ENSG00000135220.11|UGT2A3	.	-
chr5	73643317	73643317	chr5:73643317-73643317|-|ENSG00000214944.10|ARHGEF28	.	-
chr5	113567928	113567928	chr5:113567928-113567928|-|ENSG00000047188.16|YTHDC2	.	-
chr6	33255417	33255417	chr6:33255417-33255417|+|ENSG00000223501.9|VPS52	.	+
chr7	15681655	15681655	chr7:15681655-15681655|-|ENSG00000106511.6|MEOX2	.	-
chr7	44881135	44881135	chr7:44881135-44881135|+|ENSG00000146676.10|PURB	.	+
```

12 of the 54 bioinformatically found integration sites are in intergenic regions

```bash
 % grep -c "\.|-1" integration_windows.nearest_lab_verified_site_lookup.txt
12
```

9 of the 15 lab-verified annotations sites are within the integration boundaries found bioinformatically

```bash
 % grep -c "|0" integration_windows.nearest_lab_verified_site_lookup.txt
9
 % grep "|0" integration_windows.nearest_lab_verified_site_lookup.txt
chr1:71076470-71076669	chr1:71076570-71076570|+|ENSG00000132485.14|ZRANB2|0
chr13:95131648-95131847	chr13:95131748-95131748|+|ENSG00000125257.16|ABCC4|0
chr15:73062431-73062630	chr15:73062530-73062530|-|ENSG00000067141.17|NEO1|0
chr2:11332241-11332440	chr2:11332375-11332375|+|ENSG00000134318.14|ROCK2|0
chr5:73643062-73643416	chr5:73643317-73643317|-|ENSG00000214944.10|ARHGEF28|0
chr5:113567827-113568026	chr5:113567928-113567928|-|ENSG00000047188.16|YTHDC2|0
chr6:33255317-33255516	chr6:33255417-33255417|+|ENSG00000223501.9|VPS52|0
chr7:15681556-15681755	chr7:15681655-15681655|-|ENSG00000106511.6|MEOX2|0
chr7:44881036-44881235	chr7:44881135-44881135|+|ENSG00000146676.10|PURB|0
```

