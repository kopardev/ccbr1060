def get_junction_files(wildcards):
    samples=GROUP2SAMPLES[wildcards.group]
    junction_files=list()
    for s in samples:
        f=join(RESULTSDIR,s,"STAR1p","withoutChimericJunctions",s+"_p1.Chimeric.out.junction")
        junction_files.append(f)
    return junction_files

def get_verified_sites_bed(wildcards):
    return CONFIG["lab_verified_sites"][wildcards.group]

localrules: get_filtered_chimeric_junctions,make_integration_windows

rule get_filtered_chimeric_junctions:
    input:
        get_junction_files
    output:
        junctions=join(RESULTSDIR,"shRNA_integration","{group}","filtered_chimeric_junctions.tsv"),
        aggregate_junctions=join(RESULTSDIR,"shRNA_integration","{group}","filtered_chimeric_junctions.aggregate_across_samples.tsv"),
    params:
        group="{group}",
        nreads_threshold=config['shRNA_detection']['nreads_threshold'],
    shell:"""
for f in {input};do
  sample=$(basename $f|awk -F"_p1" '{{print $1}}')
  grep chr_shRNA $f | awk -F"\\t" '{{if ($1!=$4){{read[$10]++;if (read[$10]==1){{print}}}}}}' > ${{f}}.shRNA
  awk -F"\\t" -v OFS="\\t" '{{if ($1=="chr_shRNA"){{print $4,$5}}else{{print $1,$2}}}}' ${{f}}.shRNA | \
    sort | uniq -c | awk -v OFS="\\t" -v s=$sample -v t=$threshold_nreads '{{if ($1>t){{print $1,$2,$3,s}}}}'
done | sort -k1,1 -k2,2n -k3,3 > {output.junctions}

cat {output.junctions} |awk '{{array[$2"_"$3]+=$1}}END{{for (i in array) print i"\\t"array[i]}}'|sed "s@_@\\t@g" > {output.aggregate_junctions}
"""

rule make_integration_windows:
    input:
        aggregate_junctions=rules.get_filtered_chimeric_junctions.output.aggregate_junctions,
        lab_verified_bed=get_verified_sites_bed
    output:
        iwbed=join(RESULTSDIR,"shRNA_integration","{group}","integration_windows.bed"),
        iwbedannotated=join(RESULTSDIR,"shRNA_integration","{group}","integration_windows.bed.annotated.txt"),
        iwbedannotatedbed=join(RESULTSDIR,"shRNA_integration","{group}","integration_windows.annotated.bed"),
        labedannotated=join(RESULTSDIR,"shRNA_integration","{group}","lab_verified_sites.hg38.bed.annotated.txt"),
        labedannotatedbed=join(RESULTSDIR,"shRNA_integration","{group}","lab_verified_sites.annotated.bed"),
    params:
        group="{group}",
        outdir=join(RESULTSDIR,"shRNA_integration","{group}"),
        script=join(SCRIPTSDIR,"create_integration_windows.py"),
        flanking_width=config['shRNA_detection']['integration_window_flanking_width'],
        genes_bed6=config['shRNA_detection']['genes_bed6'],
    envmodules: TOOLS["bedtools"]["version"], TOOLS["ucsc"]["version"]
    shell:"""
TMPDIR="/dev/shm/{params.group}"
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi
python {params.script} {input.aggregate_junctions} | grep SiteList | \
 awk -v OFS="\\t" -v t="{params.flanking_width}" '{{print $4,$6-t,$6+t-1}}' | sed "s@,@@g" | sort -k1,1 -k2,2n > $TMPDIR/integration_window.bed
bedtools merge -i $TMPDIR/integration_window.bed > $TMPDIR/integration_windows.bed.tmp && \
 mv $TMPDIR/integration_windows.bed.tmp {output.iwbed}
bedtools intersect -loj -wa -wb -a {output.iwbed} -b {params.genes_bed6} | \
 awk -F"\\t" '{{seen[$1":"$2"-"$3]++;if (seen[$1":"$2"-"$3]==1){{print}}}}' > {output.iwbedannotated}
awk -F"\\t" -v OFS="\\t" '{{print $1":"$2"-"$3,$7"|"$NF}}' {output.iwbedannotated} > {params.outdir}/integration_windows.annotation_lookup.txt
awk '{{print $1}}' $TMPDIR/integration_windows.annotation_lookup.txt | sed "s@:@\\t@g" | sed "s@-@\\t@g"  > $TMPDIR/integration_windows.annotation_lookup.txt.a
awk '{{print $2}}' $TMPDIR/integration_windows.annotation_lookup.txt >  $TMPDIR/integration_windows.annotation_lookup.txt.b
paste  $TMPDIR/integration_windows.annotation_lookup.txt.a  $TMPDIR/integration_windows.annotation_lookup.txt.b | \
 awk -v OFS="\\t" '{{print $1,$2,$3,$1":"$2"-"$3"|"$4,".","."}}' > {output.iwbedannotatedbed}
bedtools intersect -loj -wa -wb -a {input.lab_verified_bed} -b {params.genes_bed6} | \
 awk -F"\\t" '{{seen[$1":"$2"-"$3]++;if (seen[$1":"$2"-"$3]==1){{print}}}}' > {output.labedannotated}
awk -F"\\t" -v OFS="\\t" '{{print $1,$2,$3,$1":"$2"-"$3"|"$6"|"$10,".",$6}}' {output.labedannotated} > {output.labedannotatedbed}
bedtools closest -d -a {input.lab_verified_bed} -b {output.labedannotatedbed} | \
 awk -F"\\t" -v OFS="\\t" '{{print $1":"$2"-"$3,$7"|"$10}}' | \
 awk '{{seen[$1]++;if (seen[$1]==1) {{print}}}}' > {params.outdir}/integration_windows.nearest_lab_verified_site_lookup.txt
rm -rf $TMPDIR
"""

