

rule star1p:
    input:
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2
    output:
        junction=join(RESULTSDIR,"{sample}","STAR1p","withoutChimericJunctions","{sample}_p1.Chimeric.out.junction"),
        bam=join(RESULTSDIR,"{sample}","STAR1p","withChimericJunctions","{sample}_p1.Aligned.sortedByCoord.out.bam"),
        bai=join(RESULTSDIR,"{sample}","STAR1p","withChimericJunctions","{sample}_p1.Aligned.sortedByCoord.out.bam.bai"),
    params:
        sample="{sample}",
        peorse=get_peorse,
        workdir=WORKDIR,
        outdir=join(RESULTSDIR,"{sample}","STAR1p"),
        starindexdir=STARINDEXDIR,
        alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
        gtf=GTF
    envmodules: TOOLS["star"]["version"], TOOLS["samtools"]["version"]
    threads: getthreads("star1p")
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
    --sjdbOverhang $overhang

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
    --sjdbOverhang $overhang

    # index bam file
    samtools index {output.bam}
fi
rm -rf $tmpdir
"""

rule bam2bw:
    input:
        bam=rules.star1p.output.bam,
    output:
        bw=join(RESULTSDIR,"{sample}","STAR1p","withChimericJunctions","{sample}.bw"),
        rbw=join(RESULTSDIR,"{sample}","STAR1p","withChimericJunctions","{sample}.fwd.bw"),
        fbw=join(RESULTSDIR,"{sample}","STAR1p","withChimericJunctions","{sample}.rev.bw"),
    params:
        script=join(SCRIPTSDIR,"bam_to_strand_specific_bigwigs.bash"),
        sample="{sample}"
    shell:"""
set -e -x -o pipefail
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then TMPDIR="/lscratch/${{SLURM_JOB_ID}}";else TMPDIR="/dev/shm";fi
bn=$(dirname {output.bw})
cd $bn
bash {params.script} --sample {params.sample} --bam {input.bam}
"""
    