# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: count.smk
# time: 2021/06/04


def find_all_library(wildcards):
    return [
        config['workspace'] + f'/samples/{wildcards.sample}/mapping/{library}_Aligned.out.bam'
        for library in config['samples'][wildcards.sample]['fastq']
    ]


rule count:
    output:
        config['workspace'] + '/samples/{sample}/count/{sample}_raw_count.txt'
    input:
        find_all_library
    log:
        config['workspace'] + '/log/aggregate/{sample}/{sample}_htseq_count.log'
    params:
        gtf=config['genome']['gtf']
    shell:
        'htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m intersection-nonempty'
        ' {input} {params.gtf} > {output} 2>{log}'