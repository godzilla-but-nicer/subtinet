configfile: 'workflow/config.yaml'

rule delta6_network:
    input:
        gml='data/arrietta-orritz-networks/arrietta-orritz-{size}.graphml',
        d6genes='data/genelists/delta6_168_cds_matched.csv',
        spor='data/genelists/SWxWWxRS_sporulation_genes.csv'
    output:
        'data/delta6-networks/unannotated/delta6-{size}.graphml'
    script:
        'scripts/delta6_network.py'

rule rnaseq_network:
    input:
        gml='data/delta6-networks/unannotated/delta6-{size}.graphml',
        csvs=expand('data/significant-deseq2/{treatment}_IPTG_vs_{treatment}_NT.FDR_5.xlsx',
                    treatment=config['rnaseq'])
    output:
        'data/delta6_networks/annotated/delta6-{size}-rnaseq.graphml'
    script:
        'scripts/rnaseq_network.py'

rule plot_degree_distribution:
    input:
        gml='data/delta6_networks/annotated/delta6-{size}-rnaseq.graphml'
    output:
        'plots/degree_distribution/delta6-{size}-degree.png'
    script:
        'scripts/plot_degree_distribution.py'