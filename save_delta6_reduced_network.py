import networkx as nx
import pandas as pd

# load graph and genes in delta6 genome
G = nx.read_graphml('data/net--child.graphml')
delta6 = pd.read_csv('data/delta6_168_cds_matched.csv')
# syntax for the locus tags is different for A-O network
keep_genes = delta6['locus_tag.168'].str.replace('_', '')
# thousands of lookups will be faster on a set than a list(n vs n*m)
gene_set = set(keep_genes)

# make list of genes to drop
drop_genes = []
for node in G.nodes(data=True):
    if node[1]['name'] not in gene_set:
        drop_genes.append(node[0]) # equiv to node[1]['SUID']

# drop em!
G.remove_nodes_from(drop_genes)

nx.write_graphml(G, 'data/arrietta-orritz-delta6-genome-reduced.graphml')