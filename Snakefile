
rule all:
    input: "input_assembly/GCA_000180675.11_ASM18067v1.fq"
# assembly here: https://www.ncbi.nlm.nih.gov/assembly/GCA_000180675.1/
rule download_genome:
    conda: "env-wget.yml"
    output: "input_assembly/GCA_000180675.11_ASM18067v1_genomic.fna.gz"
    shell:'''
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/180/675/GCA_000180675.1_ASM18067v1/GCA_000180675.1_ASM18067v1_genomic.fna.gz -O {output}
    '''

# unzip
rule unzip:
    input: "input_assembly/GCA_000180675.11_ASM18067v1_genomic.fna.gz"
    output: "input_assembly/GCA_000180675.11_ASM18067v1.fq"
    threads: 1
    resources:
         mem_mb=2000,
         time=2880
    shell:'''
    gunzip -c {input} > {output}
    '''

#rule grepsbf1:
#    input: "downloadedseq.fq"
#    output: "fqfile_of_bait_locs.fq"
#    shell:'''
#    #do some thing here
#    '''

#rule bait_trim:
#    input: "fqfile_of_bait_locs.fq"
#    output: "filtered_flanking_regions_stickle"
#    shell:'''
#    #do some thing here
#    '''

#rule samtools_index:
#    input:
#        "sorted_reads/{sample}.bam"
#    output:
#        "sorted_reads/{sample}.bam.bai"
#    shell:
#        "samtools index {input}"

#rule bcftools_call:
#    input:
#        fa="data/genome.fa",
#	      bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
#	      bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
#    output:
#         "calls/all.vcf"
#    shell:
#        "bcftools mpileup -f {input.fa} {input.bam} | "
#        "bcftools call -mv - > {output}"
