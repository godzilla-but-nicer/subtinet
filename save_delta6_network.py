import networkx as nx
import pandas as pd
from tqdm import tqdm

# load graph and genes in delta6 genome
G = nx.read_graphml('data/arrieta-orritz_data/Dataset_EV5/net.graphml')
delta6 = pd.read_csv('data/delta6_168_cds_matched.csv')
sporulation = pd.read_csv('data/SWxWWxRS_sporulation_genes.csv')
# syntax for the locus tags is different for A-O network
keep_genes = delta6['locus_tag.168'].str.replace('_', '')
spor_genes = sporulation['gene'].str.replace('_', '').values
d6_tags = delta6['locus_tag.d6']
name_dict = {t168 : d6 for t168, d6 in zip(keep_genes, d6_tags)}
# thousands of lookups will be faster on a set than a list(n vs n*m)
gene_set = set(keep_genes)
spor_set = set(spor_genes)

# make list of genes to drop
drop_genes = []
hits = 0
for node in tqdm(G.nodes(data=True)):
    if node[1]['name'] not in gene_set:
        drop_genes.append(node[0]) # equiv to node[1]['SUID']
    else:
        node[1]['locus_tag.d6'] = name_dict[node[1]['name']]
        # check if node is sporulation gene or not
        if node[1]['name'] in spor_set:
            hits += 1
            node[1]['sporulation'] = 1
        else:
            node[1]['sporulation'] = 0

# drop em!
G.remove_nodes_from(drop_genes)

nx.write_graphml(G, 'data/arrietta-orritz-delta6-genome.graphml')
print(hits)