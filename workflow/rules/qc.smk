rule fastqc:
# """
# Run FASTQC on:
# * Raw fastqs
# * Trimmed fastqs
# """
    input:
        expand(join(WORKDIR,"fastqs","{sample}.R1.fastq.gz"),sample=SAMPLES),
        expand(join(WORKDIR,"fastqs","{sample}.R2.fastq.gz"),sample=SAMPLES),
        expand(rules.cutadapt.output.of1,sample=SAMPLES),
        expand(rules.cutadapt.output.of2,sample=SAMPLES),
    output:
        expand(join(QCDIR,"fastqc","{sample}.R1_fastqc.zip"), sample=SAMPLES),
        expand(join(QCDIR,"fastqc","{sample}.R2_fastqc.zip"), sample=SAMPLES),
        expand(join(QCDIR,"fastqc","{sample}.R1.trim_fastqc.zip"), sample=SAMPLES),
        expand(join(QCDIR,"fastqc","{sample}.R2.trim_fastqc.zip"), sample=SAMPLES),
    params:
        outdir=join(QCDIR,"fastqc"),
    threads: getthreads("fastqc")
    envmodules: TOOLS["fastqc"]["version"]
    shell: """
fastqc {input} -t {threads} -o {params.outdir};
    """ 

#########################################################

rule fastq_screen:
    """
    Quality-control step to screen for different sources of contamination.
    FastQ Screen compares your sequencing data to a set of different reference
    genomes to determine if there is contamination. It allows a user to see if
    the composition of your library matches what you expect.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        FastQ Screen report and logfiles
    """
    input:
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2,
    output:
        out1=join(QCDIR,"FQscreen","{sample}.R1.trim_screen.txt"),
        out2=join(QCDIR,"FQscreen","{sample}.R1.trim_screen.png"),
        out3=join(QCDIR,"FQscreen","{sample}.R2.trim_screen.txt"),
        out4=join(QCDIR,"FQscreen","{sample}.R2.trim_screen.png"),
    params:
        outdir = join(QCDIR,"FQscreen"),
        fastq_screen_config=FASTQ_SCREEN_CONFIG,
    threads: getthreads("fastq_screen")  
    envmodules: TOOLS["bowtie2"]["version"],TOOLS["fastq_screen"]["version"]
    shell: """
set -e -x -o pipefail
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then TMPDIR="/lscratch/${{SLURM_JOB_ID}}";else TMPDIR="/dev/shm";fi
fastq_screen --conf {params.fastq_screen_config} --outdir {params.outdir} \
    --threads {threads} --subset 1000000 \
    --aligner bowtie2 --force {input.R1} {input.R2}
    """

#########################################################

rule multiqc:
    input:
        expand(join(QCDIR,"fastqc","{sample}.R1_fastqc.zip"), sample=SAMPLES),
        expand(join(QCDIR,"fastqc","{sample}.R2_fastqc.zip"), sample=SAMPLES),
        expand(join(QCDIR,"fastqc","{sample}.R1.trim_fastqc.zip"), sample=SAMPLES),
        expand(join(QCDIR,"fastqc","{sample}.R2.trim_fastqc.zip"), sample=SAMPLES),
    output:
        join(RESULTSDIR,"QC","multiqc_report.html"),
    params:
        qcdir=QCDIR,
        workdir=WORKDIR,
    envmodules: TOOLS["multiqc"]["version"]
    shell:"""
set -e -x -o pipefail
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then TMPDIR="/lscratch/${{SLURM_JOB_ID}}";else TMPDIR="/dev/shm";fi
cd {params.qcdir}
multiqc --interactive --force {params.workdir}
"""

#########################################################