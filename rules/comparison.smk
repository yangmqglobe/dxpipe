# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: comparison.smk
# time: 2021/04/22
import os


rule plan:
    output:
        config['workspace'] + '/comparisons/{comparison}/{comparison}_plan.txt'
    run:
        import pandas as pd
        # construct metadata
        meta = {sample: meta['meta'] for sample, meta in config['samples'].items()}
        meta = pd.DataFrame.from_dict(meta, orient='index')
        # filter out unuse samples
        use = pd.Series(True, index=meta.index)
        if 'exclude' in config['comparisons'][wildcards.comparison]:
            for key, value in config['comparisons'][wildcards.comparison]['exclude'].items():
                if isinstance(value, list):
                    use = use & (~meta[key].isin(value))
                else:
                    use = use & (meta[key] != value)
        # construct comparison
        condition = config['comparisons'][wildcards.comparison]['condition']
        numerator = config['comparisons'][wildcards.comparison]['numerator']
        denominator = config['comparisons'][wildcards.comparison]['denominator']
        comparison = pd.Series('unuse', index=meta.index, name='condition')
        comparison[(meta[condition] == numerator) & use] = 'numerator'
        comparison[(meta[condition] == denominator) & use] = 'denominator'
        comparison.to_csv(output[0], sep='\t')


rule compare:
    output:
        config['workspace'] + '/comparisons/{comparison}/{comparison}_{set}_result.txt'
    input:
        rules.plan.output,
        config['workspace'] + '/aggregate/all_sample_{set}_raw_counts.txt'
    params:
        script=lambda wildcards: f'{BASE_DIR}/tools/compare.R'
    shell:
        'Rscript {params.script} {input} {output}'



rule gsea_cls:
    output:
        config['workspace'] + '/comparisons/{comparison}/{comparison}_gsea.cls'
    input:
        rules.plan.output
    run:
        import pandas as pd
        plan = pd.read_table(input[0], index_col=0)
        plan = plan[plan['condition'] != 'unuse'].copy()
        plan = plan.sort_values('condition', ascending=False)
        num = plan.shape[0]

        condition = config['comparisons'][wildcards.comparison]['condition']
        numerator = config['comparisons'][wildcards.comparison]['numerator']
        denominator = config['comparisons'][wildcards.comparison]['denominator']

        condition = [
            config['samples'][sample]['meta'][condition]
             for sample, position in plan.itertuples()
        ]

        with open(output[0], 'w') as f:
            f.write(f'{num}\t2\t1\n')
            f.write(f'# {numerator}\t{denominator}\n')
            f.write('\t'.join(condition))


rule gsea_res:
    output:
        config['workspace'] + '/comparisons/{comparison}/{comparison}_gsea.txt'
    input:
        res=rules.merge_count.output.qnorm_ccds,
        plan=rules.plan.output
    run:
        import pandas as pd
        plan = pd.read_table(input[1], index_col=0)
        plan = plan[plan['condition'] != 'unuse'].copy()
        plan = plan.sort_values('condition', ascending=False)

        res = pd.read_table(input[0], index_col=0)
        res['Description'] = 'na'

        res[['Description'] + plan.index.tolist()].to_csv(
            output[0], sep='\t', index_label='NAME'
        )


rule gsea:
    output:
        directory(config['workspace'] + '/comparisons/{comparison}/{comparison}_gsea_{geneset}')
    input:
        res=rules.gsea_res.output,
        plan=rules.gsea_cls.output
    log:
        config['workspace'] + '/log/comparisons/{comparison}/{comparison}_gsea_{geneset}.log'
    params:
        geneset=lambda wildcards: config['genome']['geneset'][wildcards.geneset]
    shell:
        'gsea-cli.sh GSEA -res {input.res} -collapse false -permute phenotype'
        ' -rpt_label {wildcards.comparison}_{wildcards.geneset} -out {output}'
        ' -cls {input.plan} -gmx {params.geneset} >{log} 2>&1'
