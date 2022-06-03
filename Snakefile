#import pandas as pd
#TMPDIR = "/scratch/rapeek"


rule all:
    input: "outputs/sbf1_stickle_baits.txt"

# assembly here: https://www.ncbi.nlm.nih.gov/assembly/GCA_000180675.1/
# genbank: GCA_000180675.1
rule download:
    input: "URL/GCA_000180675.1.fna"
    output: "file/GCA_000180675.1.fna"
    shell:'''
    curl -L {input} -o {output}
    '''

rule unzip:
    input: "GCA_000180675.1.fq.gz"
    output: "assembly/GCA_000180675.1.fastq"
    threads: 1
    resources:
        mem_mb=2000,
        time=2880
    shell:'''
    gunzip -c {input} > {output}
    '''

rule grepsbf1:
    input: "downloadedseq.fq"
    output: "fqfile_of_bait_locs.fq"
    shell:'''
    #do some thing here
    '''

rule bait_trim:
    input: "fqfile_of_bait_locs.fq"
    output: "filtered_flanking_regions_stickle"
    shell:'''
    #do some thing here
    '''

rule bcftools_call:
	input:
	    fa="data/genome.fa"
	    bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
	    bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
	output:
	    "calls/all.vcf"
        shell:
	    "bcftools mpileup -f {input.fa} {input.bam} | "
	    "bcftools call -mv - > {output}"
