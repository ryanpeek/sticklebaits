
rule all:
    input:
        "input_assembly/GCA_000180675.11_ASM18067v1_genomic.fna.gz",
        "output/baits_sbf1_120bp_flank.txt"

# assembly here: https://www.ncbi.nlm.nih.gov/assembly/GCA_000180675.1/
rule download_genome:
    conda: "env-wget.yml"
    output: "input_assembly/GCA_000180675.11_ASM18067v1_genomic.fna.gz"
    shell:'''
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/180/675/GCA_000180675.1_ASM18067v1/GCA_000180675.1_ASM18067v1_genomic.fna.gz -O {output}
    '''

# unzip
#rule unzip:
#    input: "input_assembly/GCA_000180675.11_ASM18067v1_genomic.fna.gz"
#    output: "input_assembly/GCA_000180675.11_ASM18067v1.fq"
#    threads: 1
#    resources:
#         mem_mb=2000,
#         time=2880
#    shell:"""
#    gunzip -c {input} > {output}
#    """

rule grep_sbf1:
    input: "input_assembly/GCA_000180675.11_ASM18067v1_genomic.fna.gz"
    output: "output/bait_sbf1_seqs_to_clip.fa.gz"
    conda: "env-seqkit.yml"
    params: "CCTGCAGG"
    shell:"""
    seqkit grep -s -i -p {params} {input} | seqkit seq -n -s -u -w 0 -o {output}
    """

rule bait_trim:
    input: "output/bait_sbf1_seqs_to_clip.fa.gz"
    output: "output/baits_sbf1_120bp_flank.txt"
    params: "CCTGCAGG"
    shell:"""
    zgrep -E -o ".{{0,0}}{params}.{{0,120}}" {input} > {output}
    """

# blast to see if there's a microbial match
# could use sourmash gather but not read by read
# probably blast to make sure there's not off target matches
# could grab some sequence data from SRA and check for snp variability

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
