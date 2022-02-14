def get_junction_files(wildcards):
    """
    @Description:
    Return a list of STAR chimeric.out.junctions files for all replicates of a sample.
    """
    samples=GROUP2SAMPLES[wildcards.group]
    junction_files=list()
    for s in samples:
        f=join(RESULTSDIR,s,"STAR","withoutChimericJunctions",s+"_p1.Chimeric.out.junction")
        junction_files.append(f)
    return junction_files

def get_verified_sites_bed(wildcards):
    return CONFIG["lab_verified_sites"][wildcards.group]

localrules: get_filtered_chimeric_junctions, make_integration_windows, annotate_integration_windows, nearest_lab_verified_site_lookup, create_count_matrix, annotate_junctions_init, make_excel_workbooks, make_counts_excel_workbook, plot_shRNA_coordinate_distribution, plot_per_iw_distribution, make_annotated_junctions_files_megalist

rule get_filtered_chimeric_junctions:
    """
    @Description
    Filter STAR chimeric junction output to keep:
    * one read per readid
    * chimeric junction involving shRNA as donor or acceptor
    * each junction should have atleast nreads_threshold number of read coverage
    @Input:
        List of STAR chimeric junction files for all replicates of a sample.
        Each chimeric junction file is tab-delimited with the following columns:
        1.	Chromosome of the donor.
        2.	First base of the intron of the donor (1-based).
        3.	Strand of the donor.
        4.	Chromosome of the acceptor.
        5.	First base of the intron of the acceptor (1-based).
        6.	Strand of the acceptor.
        7.	N/A—not used, but is present to be compatible with other tools. It will always be 1.
        8.	N/A—not used, but is present to be compatible with other tools. It will always be *.
        9.	N/A—not used, but is present to be compatible with other tools. It will always be *.
        10.	Read name.
        11.	First base of the first segment, on the + strand.
        12.	CIGAR of the first segment.
        13.	First base of the second segment.
        14.	CIGAR of the second segment.
    @Output:
        @param: junctions <tsv_file>
        TSV output file reporting shRNA specific chimeric junctions with the following tabs:
            1. read depth 
            2. non-shRNA chromosome in the chimeric junction
            3. coordinate on chromosome in col 2.
            4. replicate name
        @param: aggregate_junctions <tsv_file>
        TSV file derived from the above junctions output file with the following tabs:
            1. non-shRNA chromosome in the chimeric junction
            2. coordinate on chromosome in col 1.
            3. read coverage aggregated across all replicates in which this 
            non-shRNA+coordinate combination was found
        @param: host_donor_acceptor_sites_bed <bed_file>
        BED file derived from the above aggregate_junctions file. This give some idea where the
        shRNA could be integrating with the host genome.
    """
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
    """
    @Description:
    Use chimeric junctions aggregated accross multiple replicates and generate a bed file of
    plausible integration windows. Merge overlapping windows. Score windows by number of integration
    sites in the window.
    @Input:
        @param: aggregate_junctions <tsv_file>
        TSV file derived from the above with the following tabs:
            1. non-shRNA chromosome in the chimeric junction
            2. coordinate on chromosome in col 1.
            3. read coverage aggregated across all replicates in which this 
            non-shRNA+coordinate combination was found
        @param: lab_verified_bed <bed_file>
        User-defined bed file of lab-verified shRNA integration sites
        @param: host_donor_acceptor_sites_bed <bed_file>
        BED file derived from the above aggregate_junctions file. This gives some idea where the
        shRNA could be integrating with the host genome.
    @Output:
        @param: iwbed <bed_file>
        BED file of shRNA integration sites in host with 
        score = number of integration sites in each integration window.
    """
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
    envmodules: TOOLS["bedtools"]["version"], TOOLS["ucsc"]["version"]
    shell:"""
TMPDIR="/dev/shm/{params.group}"
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi
# create_integration_windows.py collapse integration sites which are within 100bp of each
# other into a single coordinate with aggregated site sites
# grep-ing for "Sitelist" in the output will give lines with following space-separated 
# relevant columns
# $4 = chr
# $6 = bestcoordinate or representative coordinate for all integration sites in the same 100 bp window
# $8 = number of sites in this sitelist .. number of integration sites in this 100 bp window
# $10 = cumulative sum of all reads of all sites in this 100 bp window
python {params.script} {input.aggregate_junctions} | grep SiteList | \
 awk -v OFS="\\t" -v t="{params.flanking_width}" '{{print $4,$6-t,$6+t-1}}' | sed "s@,@@g" | sort -k1,1 -k2,2n > $TMPDIR/integration_windows.bed
# merge intersecting integration windows
bedtools merge -i $TMPDIR/integration_windows.bed > $TMPDIR/integration_windows.bed.tmp
mv $TMPDIR/integration_windows.bed.tmp $TMPDIR/integration_windows.bed
# add score to the integration windows
# score = number of host-shrna integrations detected in the window
bedtools intersect -wa -a $TMPDIR/integration_windows.bed -b {input.host_donor_acceptor_sites_bed} | \
 sort|uniq -c| \
 awk -v OFS="\\t" '{{print $2,$3,$4,".",$1,"."}}' > {output.iwbed}
rm -rf $TMPDIR
"""

rule annotate_integration_windows:
    """
    @Description:
        Annotate the integration windows, both the bioinformatically detected ones and the
        lab-verified ones. "genes.bed" files created before running the pipeline for hg38
        is used for extracting the annotations.
    @Inputs:
        @param: iwbed <bed_file>
            Integration windows bed6 file. score = number of integration sites in the window.
        @param: lvbed <bed_file>
            Lab-verified integration sites in bed format.
    @Outputs:
        @param: iwbedannotationlookup <tsv_file>
            TSV file with columns:
            1: region to be annotated: <chr>:<start>-<end>|<score>|<strand> format
            2: annotation from genes.bed: <chr>:<start>-<end>|<ensembl_id>_<genename>|<strand> format
        @param: iwbedannotatedbed <bed_file>
            BED6 file with name column of the format:
            <chr>:<start>-<end>|<score>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <region to-be annotated>#<annotation1>#<annotation2> 
            i.e. "#"separated values with multiple annotations 2nd item onwards
        @param: lvbedannotationlookup <tsv_file>
            TSV file with columns:
            1: region to be annotated: <chr>:<start>-<end>|<score>|<strand> format
            2: annotation from genes.bed: <chr>:<start>-<end>|<ensembl_id>_<genename>|<strand> format
        @param: lvbedannotatedbed <bed_file>
            BED6 file with name column of the format:
            <chr>:<start>-<end>|<score>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <region to-be annotated>#<annotation1>#<annotation2> 
            i.e. "#"separated values with multiple annotations 2nd item onwards
    """
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
    """
    @Description:
    Find closest lab-verfied integration site to each site and report the distance to it in bp
    @Inputs:
        @param: iw <bed6_file>
        @param: lv <bed6_file>
        Both inputs are annotated BED6 TSV with "name" or 4th column as:
        <chr>:<start>-<end>|<score>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <region to-be annotated>#<annotation1>#<annotation2> 
        i.e. "#"separated values with multiple annotations 2nd item onwards   
    @Outputs:
        @param: iwnearestlvlookup <tsv_file>
        columns:
        1. iw in format: <chr>:<start>-<end>|<score>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <region to-be annotated>#<annotation1>#<annotation2> 
        2. lv in format: <chr>:<start>-<end>|<score>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <region to-be annotated>#<annotation1>#<annotation2> 
        3. distance in bp.
        Note of col3 (distance in bp) = "-1", then no lab-verified integration site found on the chromosome of interest
    """
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

rule annotate_junctions_init:
    input:
        get_junction_files
    output:
        junction_files_list=join(RESULTSDIR,"shRNA_integration","{group}","junction_files.lst")
    params:
        group="{group}",
    shell:"""
for i in {input};do
    echo $i
done > {output.junction_files_list} 
"""


rule annotate_junctions_and_aggregate_OST_calls:
    input:
        junction_files_list=rules.annotate_junctions_init.output.junction_files_list,
        iwbedannotatedbed=rules.annotate_integration_windows.output.iwbedannotatedbed
    output:
        annotated_junctions_files_list=join(RESULTSDIR,"shRNA_integration","{group}","annotated_junctions_files.lst"),
        detailed_counts_files_list=join(RESULTSDIR,"shRNA_integration","{group}","detailed_counts_files.lst"),
        aggregate_OST_calls_files=join(RESULTSDIR,"shRNA_integration","{group}","aggregate_OST_calls_files.lst")
    params:
        group="{group}",
        outdir=join(RESULTSDIR,"shRNA_integration","{group}"),
        shRNA_integration_windows_bed=config['shRNA_detection']['shRNA_integration_windows_bed'],
        script1=join(SCRIPTSDIR,"annotate_junctions.py"),
        script2=join(SCRIPTSDIR,"append_junction_type_info.py")
    envmodules: TOOLS["bedtools"]["version"]
    shell:"""
set -euxo pipefail
cd {params.outdir}
while read junction_file;do
bn=$(basename $junction_file)
sn=$(echo $bn|awk -F"_p1.Chimeric" '{{print $1}}')
out1="{params.outdir}/${{sn}}.annotated.junctions"
out2="{params.outdir}/${{sn}}.detailed_counts_matrix.tsv"
out3="{params.outdir}/${{sn}}.aggregate_OST_calls.tmp.tsv"
out4="{params.outdir}/${{sn}}.aggregate_OST_calls.tsv"
echo "$out1" >> {output.annotated_junctions_files_list}
echo "$out2" >> {output.detailed_counts_files_list}
echo "$out4" >> {output.aggregate_OST_calls_files}
python {params.script1} $junction_file {input.iwbedannotatedbed} {params.shRNA_integration_windows_bed} \
    $out1 $out2
echo -ne "nreads\\t" > $out3
head -n1 $out1 |cut -f22- >> $out3
grep yes $out1 |cut -f22-|sort|uniq -c |awk -v OFS="\\t" '{{print $1,$2,$3,$4,$5}}'|sort -k1,1nr >> $out3
python {params.script2} $out1 $out3 \
    $out4
rm -f $out3
done < {input.junction_files_list}
while read a;do ls -larth $a;done < {output.annotated_junctions_files_list}
while read a;do ls -larth $a;done < {output.detailed_counts_files_list}
while read a;do ls -larth $a;done < {output.aggregate_OST_calls_files}
"""

rule make_excel_workbooks:
    input:
        detailed_counts_files_list=rules.annotate_junctions_and_aggregate_OST_calls.output.detailed_counts_files_list,
        aggregate_OST_calls_files=rules.annotate_junctions_and_aggregate_OST_calls.output.aggregate_OST_calls_files,
    output:
        detailed_counts_excel=join(RESULTSDIR,"shRNA_integration","{group}","{group}.detailed_counts.xlsx"),
        aggregate_OST_calls_excel=join(RESULTSDIR,"shRNA_integration","{group}","{group}.aggregate_OST_calls.xlsx"),
    params:
        group="{group}",
        script=join(SCRIPTSDIR,"create_excel_workbook.py")
    shell:"""
set -euxo pipefail
python {params.script} {output.detailed_counts_excel} $(cat {input.detailed_counts_files_list} | tr '\\n' ' ')
python {params.script} {output.aggregate_OST_calls_excel} $(cat {input.aggregate_OST_calls_files} | tr '\\n' ' ')
"""

rule make_counts_excel_workbook:
    input:
        expand(join(RESULTSDIR,"shRNA_integration","{group}","counts.tsv"),group=GROUPS)
    output:
        join(RESULTSDIR,"shRNA_integration","integration_sites.counts.xlsx")
    params:
        script=join(SCRIPTSDIR,"create_excel_workbook.py")
    shell:"""
set -euxo pipefail
python {params.script} {output} {input}
"""

rule plot_shRNA_coordinate_distribution:
    input:
        annotated_junctions_files_list=rules.annotate_junctions_and_aggregate_OST_calls.output.annotated_junctions_files_list
    output:
        coordinate_tsv=join(RESULTSDIR,"shRNA_integration","{group}","shRNA_coordinates.tsv"),
        coordinate_pdf=join(RESULTSDIR,"shRNA_integration","{group}","{group}.shRNA_coordinate_distribution.pdf"),
    params:
        group="{group}",
        outdir=join(RESULTSDIR,"shRNA_integration","{group}"),
        rscript=join(SCRIPTSDIR,"plot_shRNA_coordinate_histogram.R")
    envmodules: TOOLS["R"]["version"]
    shell:"""
set -euxo pipefail
echo -ne "sample_name\\tshRNA_coordinate\\tshRNA_strand\\tOST\\n" > {output.coordinate_tsv}
while read a;do
	bn=$(basename $a)
	sn=$(echo $bn|awk -F"." '{{print $1}}')
	awk -F"\\t" -v s="$sn" -v OFS="\\t" '{{if ($1=="chr_shRNA") {{print s,$2,$3,$NF}}else{{print s,$5,$6,$NF}}}}' $a | \
        grep "yes\\|no\\|undetermined" | \
        sed "s@yes@OST@g" | sed "s@no@other@g" | sed "s@undetermined@other@g" | \
	    awk -F"\\t" -v OFS="\\t" '{{if ($3=="+") {{$3="-";print $1,$2,$3,$4}} else {{$3="+";print $1,$2,$3,$4}}}}'
done < {input.annotated_junctions_files_list} >> {output.coordinate_tsv}
Rscript {params.rscript} {output.coordinate_tsv} {output.coordinate_pdf}
"""

rule plot_per_iw_distribution:
    input:
        detailed_counts_files_list=rules.annotate_junctions_and_aggregate_OST_calls.output.detailed_counts_files_list,
    output:
        pdf=join(RESULTSDIR,"shRNA_integration","{group}","{group}.histograms_per_integration_window.pdf"),
    params:
        group="{group}",
        outdir=join(RESULTSDIR,"shRNA_integration","{group}"),
        rscript=join(SCRIPTSDIR,"plot_histograms_per_integration_window.R") 
    envmodules: TOOLS["R"]["version"]
    shell:"""
set -euxo pipefail
Rscript {params.rscript} {output.pdf} $(cat {input.detailed_counts_files_list} | tr '\\n' ' ')
"""

rule make_annotated_junctions_files_megalist:
    input:
        expand(join(RESULTSDIR,"shRNA_integration","{group}","annotated_junctions_files.lst"),group=GROUPS)
    output:
        out=join(RESULTSDIR,"annotated_junctions_files.megalst")
    shell:"""
set -euxo pipefail
cat {input} > {output}
"""

rule make_shRNA_chimeric_reads_bigwigs:
    input:
        inbam=rules.star.output.bam,
        annotated_junctions_files_list=rules.make_annotated_junctions_files_megalist.output.out
    output:
        outbam=join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","shRNA_chimeric_only","{sample}.shRNA_chimeric.bam"),
        fwdbw=join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","shRNA_chimeric_only","{sample}.shRNA_chimeric.fwd.bw"),
        revbw=join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","shRNA_chimeric_only","{sample}.shRNA_chimeric.rev.bw"),
    params:
        sample="{sample}",
        outdir=join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","shRNA_chimeric_only"),
        pyscript=join(SCRIPTSDIR,"extract_shRNA_chimeric_reads.py"),
        shscript=join(SCRIPTSDIR,"bam_to_strand_specific_bigwigs.bash"),
    envmodules: TOOLS["samtools"]["version"]
    shell:"""
set -euxo pipefail
cd {params.outdir}
sn="{params.sample}"
annotated_junctions_file=$(grep -m1 "${{sn}}.annotated.junctions" {input.annotated_junctions_files_list})
python {params.pyscript} $annotated_junctions_file {input.inbam} {output.outbam}
samtools index {output.outbam}
bash {params.shscript} --sample {params.sample}.shRNA_chimeric --bam {output.outbam} 
"""