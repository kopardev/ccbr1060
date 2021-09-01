def get_junction_files(wildcards):
    samples=GROUP2SAMPLES[wildcards.group]
    junction_files=list()
    for s in samples:
        f=join(RESULTSDIR,s,"STAR","withoutChimericJunctions",s+"_p1.Chimeric.out.junction")
        junction_files.append(f)
    return junction_files

def get_verified_sites_bed(wildcards):
    return CONFIG["lab_verified_sites"][wildcards.group]

localrules: get_filtered_chimeric_junctions,make_integration_windows,annotate_integration_windows,nearest_lab_verified_site_lookup,create_count_matrix

rule get_filtered_chimeric_junctions:
    input:
        get_junction_files
    output:
        junctions=join(RESULTSDIR,"shRNA_integration","{group}","filtered_chimeric_junctions.tsv"),
        aggregate_junctions=join(RESULTSDIR,"shRNA_integration","{group}","filtered_chimeric_junctions.aggregate_across_samples.tsv"),
        host_donor_acceptor_sites_bed=join(RESULTSDIR,"shRNA_integration","{group}","{group}.host_donor_acceptor_sites.bed")
    params:
        group="{group}",
        nreads_threshold=config['shRNA_detection']['nreads_threshold'],
    shell:"""
for f in {input};do
  sample=$(basename $f|awk -F"_p1" '{{print $1}}')
  grep chr_shRNA $f | awk -F"\\t" '{{if ($1!=$4){{read[$10]++;if (read[$10]==1){{print}}}}}}' > ${{f}}.shRNA
  awk -F"\\t" -v OFS="\\t" '{{if ($1=="chr_shRNA"){{print $4,$5}}else{{print $1,$2}}}}' ${{f}}.shRNA | \
    sort | uniq -c | awk -v OFS="\\t" -v s=$sample -v t={params.nreads_threshold} '{{if ($1>t){{print $1,$2,$3,s}}}}'
done | sort -k1,1n -k2,2 -k3,3n > {output.junctions}

cat {output.junctions} |awk '{{array[$2"_"$3]+=$1}}END{{for (i in array) print i"\\t"array[i]}}'|sed "s@_@\\t@g" > {output.aggregate_junctions}

cat {output.aggregate_junctions} | awk -v OFS="\\t" '{{for (i=0;i<$3;i++){{print $1,$2,$2+1}}}}' > {output.host_donor_acceptor_sites_bed}
"""

rule make_integration_windows:
    input:
        aggregate_junctions=rules.get_filtered_chimeric_junctions.output.aggregate_junctions,
        lab_verified_bed=get_verified_sites_bed,
        host_donor_acceptor_sites_bed=rules.get_filtered_chimeric_junctions.output.host_donor_acceptor_sites_bed
    output:
        iwbed=join(RESULTSDIR,"shRNA_integration","{group}","integration_windows.bed"),
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
 awk -v OFS="\\t" -v t="{params.flanking_width}" '{{print $4,$6-t,$6+t-1}}' | sed "s@,@@g" | sort -k1,1 -k2,2n > $TMPDIR/integration_windows.bed
bedtools merge -i $TMPDIR/integration_windows.bed > $TMPDIR/integration_windows.bed.tmp
mv $TMPDIR/integration_windows.bed.tmp $TMPDIR/integration_windows.bed
bedtools intersect -wa -a $TMPDIR/integration_windows.bed -b {input.host_donor_acceptor_sites_bed} | \
 sort|uniq -c| \
 awk -v OFS="\\t" '{{print $2,$3,$4,".",$1,"."}}' > {output.iwbed}
rm -rf $TMPDIR
"""

rule annotate_integration_windows:
    input:
        iwbed=rules.make_integration_windows.output.iwbed,
        lvbed=get_verified_sites_bed
    output:
        iwbedannotationlookup=join(RESULTSDIR,"shRNA_integration","{group}","integration_windows.bed.annotation_lookup.txt"),
        iwbedannotatedbed=join(RESULTSDIR,"shRNA_integration","{group}","integration_windows.annotated.bed"),
        lvbedannotationlookup=join(RESULTSDIR,"shRNA_integration","{group}","lab_verified_sites.bed.annotation_lookup.txt"),
        lvbedannotatedbed=join(RESULTSDIR,"shRNA_integration","{group}","lab_verified_sites.annotated.bed"),
    params:
        group="{group}",
        script1=join(SCRIPTSDIR,"bed_to_annotation_lookup.sh"),
        script2=join(SCRIPTSDIR,"annotation_lookup_to_annotated_bed.py"),
        genes_bed6=config['shRNA_detection']['genes_bed6'],
    shell:"""
TMPDIR="/dev/shm/{params.group}"
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

bash {params.script1} {input.iwbed} {params.genes_bed6} {output.iwbedannotationlookup}
bash {params.script1} {input.lvbed} {params.genes_bed6} {output.lvbedannotationlookup}

python {params.script2} {output.iwbedannotationlookup} > {output.iwbedannotatedbed}
python {params.script2} {output.lvbedannotationlookup} > {output.lvbedannotatedbed}

"""

rule nearest_lab_verified_site_lookup:
    input:
        iw=rules.annotate_integration_windows.output.iwbedannotatedbed,
        lv=rules.annotate_integration_windows.output.lvbedannotatedbed,
    output:
        iwnearestlvlookup=join(RESULTSDIR,"shRNA_integration","{group}","integration_windows.nearest_lab_verified_lookup.txt"),
    params:
        group="{group}",
    envmodules: TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"]
    shell:"""
TMPDIR="/dev/shm/{params.group}"
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi
bedSort {input.iw} {input.iw}
bedSort {input.lv} {input.lv}
bedtools closest -d -a {input.iw} -b {input.lv}|cut -f4,10,13 > {output.iwnearestlvlookup}
"""


rule create_count_matrix:
    input:
        iwbed=rules.make_integration_windows.output.iwbed,
        junctions=rules.get_filtered_chimeric_junctions.output.junctions,
        iwnearestlvlookup=rules.nearest_lab_verified_site_lookup.output.iwnearestlvlookup,
        nreads=rules.get_total_aligned_reads.output.nreads
    output:
        count_matrix=join(RESULTSDIR,"shRNA_integration","{group}","counts.tsv")
    params:
        group="{group}",
        outdir=join(RESULTSDIR,"shRNA_integration","{group}"),
        script1=join(SCRIPTSDIR,"make_count_matrix.py"),
        samplemanifest=config["samples"]
    envmodules: TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"]
    shell:"""
TMPDIR="/dev/shm/{params.group}"
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi
counts_files=""
for s in $(cut -f4 {input.junctions}|sort|uniq);do
    awk -F"\\t" -v OFS="\\t" -v s=$s '{{if ($NF==s) {{for (i=0;i<$1;i++){{print $2,$3,$3+1}}}} }}' {input.junctions} > {params.outdir}/${{s}}.chimeric_shRNA_donor_acceptor_sites.bed
    bedtools intersect -wa -a {input.iwbed} -b {params.outdir}/${{s}}.chimeric_shRNA_donor_acceptor_sites.bed | \
        awk -F"\\t" '{{print $1":"$2"-"$3"|"$5"|"$6}}' | sort | uniq -c | \
        awk -v OFS="\\t" '{{print $2,$1}}' > {params.outdir}/${{s}}.counts.tsv
    counts_files="$counts_files {params.outdir}/${{s}}.counts.tsv"
done

python {params.script1} {input.iwnearestlvlookup} {input.nreads} {params.samplemanifest} $counts_files > {output.count_matrix}

"""