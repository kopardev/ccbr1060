rule create_rsem_index:
    input:
        reffa=REFFA,
        gtf=GTF
    output:
        ump=join(RSEMINDEXDIR,"ref.transcripts.ump")
    params:
        rsemindexdir=RSEMINDEXDIR,
    envmodules: TOOLS['rsem']['version']
    threads: getthreads("create_rsem_index")
    shell:"""
set -euf -o pipefail
cd {params.rsemindexdir}
rsem-prepare-reference -p {threads} --gtf {input.gtf} {input.reffa} ref
rsem-generate-ngvector ref.transcripts.fa ref.transcripts
"""

rule create_bed12:
    input:
        gtf=GTF
    output:
        bed12=join(RSEMINDEXDIR,"ref.bed12")
    params:
        rsemindexdir=RSEMINDEXDIR
    envmodules:
        TOOLS['ucsc']['version'], 
    shell:"""
set -euf -o pipefail
cd {params.rsemindexdir}
genesgtf={input.gtf}
bn=$(basename $genesgtf)
gtfToGenePred -genePredExt -geneNameAsName2 $genesgtf ${{bn}}.tmp
awk -v OFS="\\t" '{{print $2,$4,$5,$1,"0",$3,$6,$7,"0",$8,$9,$10}}' ${{bn}}.tmp > {output.bed12}
rm -f ${{bn}}.tmp
"""

rule star_for_rsem:
    input:
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2
    output:
        bam=join(RESULTSDIR,"{sample}","RSEM","{sample}.Aligned.out.bam"),
        tbam=join(RESULTSDIR,"{sample}","RSEM","{sample}.Aligned.toTranscriptome.out.bam"),
    params:
        sample="{sample}",
        peorse=get_peorse,
        workdir=WORKDIR,
        outdir=join(RESULTSDIR,"{sample}","RSEM"),
        starindexdir=STARINDEXDIR,
        alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
        gtf=GTF
    envmodules: TOOLS["star"]["version"], TOOLS["samtools"]["version"]
    threads: getthreads("star")
    shell:"""
set -exuf -o pipefail
tmpdir="/dev/shm/{params.sample}"
if [ -d /lscratch/${{SLURM_JOBID}} ];then tmpdir="/lscratch/${{SLURM_JOBID}}/{params.sample}";fi
# STAR tmpdir should not exist
if [ -d $tmpdir ];then rm -rf $tmpdir;fi
if [ ! -d {params.outdir} ];then mkdir -p {params.outdir};fi
if [ "{params.peorse}" == "PE" ];then
# paired-end
    overhang=$(zcat {input} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
    echo "sjdbOverhang for STAR: ${{overhang}}"

    cd {params.outdir}
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
    --outFileNamePrefix {params.sample}. \
    --chimSegmentMin 20 \
    --chimMultimapNmax 10 \
    --chimOutType WithinBAM \
    --alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \
    --outSAMtype BAM Unsorted \
    --alignEndsProtrude 10 ConcordantPair \
    --outFilterIntronMotifs None \
    --sjdbGTFfile {params.gtf} \
    --outTmpDir $tmpdir \
    --sjdbOverhang $overhang \
    --twopassMode Basic \
    --quantMode GeneCounts TranscriptomeSAM

    # index bam file
    samtools index {output.bam}
fi
rm -rf $tmpdir
"""

rule get_rsem_counts:
    input:
        bed12=rules.create_bed12.output.bed12,
        bam=rules.star_for_rsem.output.bam,
        tbam=rules.star_for_rsem.output.tbam        
    output:
        strandinfo=join(RESULTSDIR,"{sample}","RSEM","{sample}.strandinfo"),
        gcounts=join(RESULTSDIR,"{sample}","RSEM","{sample}.RSEM.genes.results"),
        tcounts=join(RESULTSDIR,"{sample}","RSEM","{sample}.RSEM.isoforms.results"),
    envmodules: TOOLS['rseqc']['version'], TOOLS['rsem']['version'], TOOLS["samtools"]["version"]
    threads: getthreads("get_rsem_counts")
    params:
        sample="{sample}",
        rsemdir=join(WORKDIR,"rsem")
    shell:"""
set -exuf -o pipefail
cd {params.rsemdir}
# inter strandedness
samtools index -@{threads} {input.bam}
infer_experiment.py -r {input.bed12} -i {input.bam} -s 1000000 > {output.strandinfo}
# Get strandedness to calculate Forward Probability
fp=$(tail -n1 {output.strandinfo} | awk '{{if($NF > 0.75) print "0.0"; else if ($NF<0.25) print "1.0"; else print "0.5";}}')
echo "Forward Probability Passed to RSEM: $fp"
rsemindex=$(echo {input.bed12}|sed "s@.bed12@@g")
rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  \
        --bam --paired-end -p {threads}  {input.tbam} $rsemindex {params.sample} --time \
        --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files --forward-prob=${{fp}} --estimate-rspd
mv {params.sample}.genes.results {output.gcounts}
mv {params.sample}.isoforms.results {output.tcounts}
"""