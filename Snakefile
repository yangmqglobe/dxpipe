# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: Snakefile
# time: 2021/06/03
include: 'rules/preprocess.smk'
include: 'rules/mapping.smk'
include: 'rules/aggregate.smk'
include: 'rules/comparison.smk'


import os

BASE_DIR = os.path.dirname(workflow.snakefile)


rule all:
    input:
        expand(
            config['workspace'] + '/aggregate/all_sample_{set}_{out}.txt',
            set=['all', 'ccds'], out=['raw_counts', 'fpkm', 'fpkm_qnorm']
        ),
        expand(
            config['workspace'] + '/aggregate/all_sample_{set}_pcaplot.{fmt}',
            fmt=config['plot_formats'], set=['all', 'ccds']
        ),
        expand(
            config['workspace'] + '/comparisons/{comparison}/{comparison}_{set}_result.txt',
            comparison=config['comparisons'], set=['all', 'ccds']
        ),
        expand(
            config['workspace'] + '/comparisons/{comparison}/{comparison}_gsea_{geneset}',
            comparison=config['comparisons'], geneset=config['genome']['geneset']
        )
