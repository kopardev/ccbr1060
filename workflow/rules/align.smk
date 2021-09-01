localrules: get_total_aligned_reads

rule star:
    input:
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2
    output:
        junction=join(RESULTSDIR,"{sample}","STAR","withoutChimericJunctions","{sample}_p1.Chimeric.out.junction"),
        pergenecounts=join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","{sample}_p1.ReadsPerGene.out.tab"),
        bam=join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","{sample}_p1.Aligned.sortedByCoord.out.bam"),
        bai=join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","{sample}_p1.Aligned.sortedByCoord.out.bam.bai"),
    params:
        sample="{sample}",
        peorse=get_peorse,
        workdir=WORKDIR,
        outdir=join(RESULTSDIR,"{sample}","STAR"),
        starindexdir=STARINDEXDIR,
        alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
        gtf=GTF
    envmodules: TOOLS["star"]["version"], TOOLS["samtools"]["version"]
    threads: getthreads("STAR")
    shell:"""
tmpdir="/dev/shm/{params.sample}"
if [ -d /lscratch/${{SLURM_JOBID}} ];then tmpdir="/lscratch/${{SLURM_JOBID}}/{params.sample}";fi
# STAR tmpdir should not exist
if [ -d $tmpdir ];then rm -rf $tmpdir;fi
if [ ! -d {params.outdir} ];then mkdir -p {params.outdir};fi
if [ "{params.peorse}" == "PE" ];then
# paired-end
    overhang=$(zcat {input} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
    echo "sjdbOverhang for STAR: ${{overhang}}"

    # running STAR with "chimOutType Junctions"
    if [ ! -d {params.outdir}/withoutChimericJunctions ]; then mkdir -p {params.outdir}/withoutChimericJunctions;fi
    cd {params.outdir}/withoutChimericJunctions
    STAR --genomeDir {params.starindexdir} \
    --outSAMstrandField None  \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.3  \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --readFilesIn {input.R1} {input.R2} \
    --readFilesCommand zcat \
    --runThreadN {threads} \
    --outFileNamePrefix {params.sample}_p1. \
    --chimSegmentMin 20 \
    --chimMultimapNmax 10 \
    --chimOutType Junctions \
    --alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \
    --outSAMtype BAM SortedByCoordinate \
    --alignEndsProtrude 10 ConcordantPair \
    --outFilterIntronMotifs None \
    --sjdbGTFfile {params.gtf} \
    --outTmpDir $tmpdir \
    --sjdbOverhang $overhang \
    --twopassMode Basic

    # running STAR with "chimOutType WithinBAM"
    if [ ! -f {params.outdir}/withChimericJunctions ];then mkdir -p {params.outdir}/withChimericJunctions;fi
    cd {params.outdir}/withChimericJunctions
    STAR --genomeDir {params.starindexdir} \
    --outSAMstrandField None  \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.3  \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --readFilesIn {input.R1} {input.R2}\
    --readFilesCommand zcat \
    --runThreadN {threads} \
    --outFileNamePrefix {params.sample}_p1. \
    --chimSegmentMin 20 \
    --chimMultimapNmax 10 \
    --chimOutType WithinBAM \
    --alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \
    --outSAMtype BAM SortedByCoordinate \
    --alignEndsProtrude 10 ConcordantPair \
    --outFilterIntronMotifs None \
    --sjdbGTFfile {params.gtf} \
    --outTmpDir $tmpdir \
    --sjdbOverhang $overhang \
    --twopassMode Basic \
    --quantMode GeneCounts

    # index bam file
    samtools index {output.bam}
fi
rm -rf $tmpdir
"""

rule bam2bw:
    input:
        bam=rules.star.output.bam,
    output:
        fbw=join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","{sample}.fwd.bw"),
        rbw=join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","{sample}.rev.bw"),
    params:
        script=join(SCRIPTSDIR,"bam_to_strand_specific_bigwigs.bash"),
        sample="{sample}"
    shell:"""
set -e -x -o pipefail
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then TMPDIR="/lscratch/${{SLURM_JOB_ID}}";else TMPDIR="/dev/shm";fi
bn=$(dirname {output.fbw})
cd $bn
bash {params.script} --sample {params.sample} --bam {input.bam}
"""
    
rule get_total_aligned_reads:
    input:
        expand(join(RESULTSDIR,"{sample}","STAR","withChimericJunctions","{sample}_p1.ReadsPerGene.out.tab"),sample=SAMPLES)
    output:
        nreads=join(RESULTSDIR,"nreads.txt")
    params:
        strandedness=config["library_strand_info"]
    shell:"""
for f in {input};do
    sample=$(basename $f|awk -F"_p1" '{{print $1}}')
    if [ "{params.strandedness}" == "1" ];then
        nreads=$(grep -v "^N_" $f|awk '{{sum=$2+sum}}END{{print sum}}')
    elif [ "{params.strandedness}" == "2" ];then
        nreads=$(grep -v "^N_" $f|awk '{{sum=$3+sum}}END{{print sum}}')
    elif [ "{params.strandedness}" == "3" ];then
        nreads=$(grep -v "^N_" $f|awk '{{sum=$4+sum}}END{{print sum}}')
    fi
    echo -ne "$sample\t$nreads\n"
done > {output.nreads}
"""