# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: aggregate.smk
# time: 2021/06/04


def find_all_library(wildcards):
    return [
        config['workspace'] + f'/samples/{wildcards.sample}/mapping/{library}_Aligned.out.bam'
        for library in config['samples'][wildcards.sample]['fastq']
    ]


rule count:
    output:
        config['workspace'] + '/aggregate/{sample}_raw_count.txt'
    input:
        find_all_library
    log:
        config['workspace'] + '/log/aggregate/{sample}/{sample}_htseq_count.log'
    params:
        gtf=config['genome']['gtf']
    shell:
        'htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m intersection-nonempty'
        ' {input} {params.gtf} > {output} 2>{log}'


rule merge_count:
    output:
        raw=config['workspace'] + '/aggregate/all_sample_all_raw_counts.txt',
        fpkm=config['workspace'] + '/aggregate/all_sample_all_fpkm.txt',
        qnorm=config['workspace'] + '/aggregate/all_sample_all_fpkm_qnorm.txt',
        raw_ccds=config['workspace'] + '/aggregate/all_sample_ccds_raw_counts.txt',
        fpkm_ccds=config['workspace'] + '/aggregate/all_sample_ccds_fpkm.txt',
        qnorm_ccds=config['workspace'] + '/aggregate/all_sample_ccds_fpkm_qnorm.txt'
    input:
        expand(
            config['workspace'] + '/aggregate/{sample}_raw_count.txt',
            sample=config['samples']
        )
    params:
        names=expand('{sample}', sample=config['samples'])
    run:
        import pandas as pd
        import numpy as np
        import gtfparse
        import qnorm

        counts = [
            pd.read_table(file, names=['gene_id', name], index_col=0, squeeze=True)
            for name, file in zip(params.names, input)
        ]
        counts = pd.concat(counts, axis=1)
        counts = counts.loc[~counts.index.str.startswith('__')]
        counts.to_csv(output.raw, sep='\t', index_label='')

        gtf = gtfparse.read_gtf(config['genome']['gtf'])
        genes = gtf.loc[
            (gtf['feature'] == 'gene') & (gtf['gene_type'] == 'protein_coding'), 'gene_id'
        ].values

        total = counts.loc[genes].sum(0)

        def exon_length(exon):
            exon = exon.sort_values(['start', 'end'])
            exon['group'] = (exon['start'] > exon['end'].shift()).cumsum()
            width = exon.groupby('group').agg({'start': 'min', 'end': 'max'})
            return (width['end'] - width['start']).sum()

        exons = gtf[gtf['feature'] == 'exon']
        length = exons.groupby('gene_id').apply(exon_length)

        fpkm = np.log2((counts * 10 ** 9 / total).div(length, axis=0) + 1)
        fpkm.to_csv(output.fpkm, sep='\t', index_label='')

        norm = qnorm.quantile_normalize(fpkm)
        norm.to_csv(output.qnorm, sep='\t', index_label='')

        ccds = gtf[(gtf['feature'] == 'transcript') & (gtf['tag'].str.contains('CCDS'))].copy()
        ccds['length'] = ccds['end'] - ccds['start']
        ccds = ccds.sort_values('length').drop_duplicates('gene_name', keep='last')

        gene_id = ccds['gene_id'].values
        gene_name = ccds['gene_name'].values

        counts.loc[gene_id].set_index(gene_name).to_csv(output.raw_ccds, sep='\t', index_label='')
        fpkm.loc[gene_id].set_index(gene_name).to_csv(output.fpkm_ccds, sep='\t', index_label='')
        norm.loc[gene_id].set_index(gene_name).to_csv(output.qnorm_ccds, sep='\t', index_label='')


rule comparisons:
    output:
        config['workspace'] + '/aggregate/comparisons.txt',
        config['workspace'] + '/aggregate/metadata.txt'
    input:
        rules.merge_count.output
    run:
        import pandas as pd

        meta = {sample: meta['meta'] for sample, meta in config['samples'].items()}
        meta = pd.DataFrame.from_dict(meta, orient='index')
        meta.to_csv(output[1], sep='\t')

        meta_cols = list({
            comparison['condition'] for comparison in config['comparisons'].values()
        })
        if len(meta_cols) == 0:
            raise ValueError('comparisons not defind!')
        meta = meta[meta_cols]

        meta.to_csv(output[0], sep='\t')


rule pcaplot:
    output:
        config['workspace'] + '/aggregate/all_sample_{set}_pcaplot.{fmt}',
    input:
        rules.comparisons.output[0],
        config['workspace'] + '/aggregate/all_sample_{set}_fpkm_qnorm.txt'
    params:
        script=lambda wildcards: f'{BASE_DIR}/tools/pcaplot.R'
    shell:
        'Rscript {params.script} {input} {output}'
