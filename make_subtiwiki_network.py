import pandas as pd
import networkx as nx
import numpy as np
from sklearn.preprocessing import LabelEncoder
gene_list = pd.read_csv('data/SWxWWxRS_sporulation_genes.csv')
regulations = pd.read_csv('data/regulations.csv')

# First we'll make a simplified list of modes
# in this list anti termination is not really the same thing as activation
# I'm not sure I'm crazy about "indirect positive regulation" appearing as a link
activating_modes = ['activation', 'sigma factor', 'positive regulation', 
                  'transcriptional activation', 'transcription activation',
                  'antitermination', 'indirect positive regulation', 
                  'mRNA binding', 'processive antitermination']
repressing_modes = ['repression', 'transcription repression', 
                    'negative autoregulation', 'attenuation', 'auto-repression',
                    'anti-activation', 'termination', 'antisense RNA',
                    'negative regulation', 'autorepression', 
                    'mRNA destabilization', 'transcription termination']
unknown_modes = [mode for mode in regulations['mode'].unique() \
                    if mode not in activating_modes and mode not in repressing_modes]

# combine these to mapping dictionary and add column to df
mode_map = {mode : 'upreg' for mode in activating_modes }
mode_map.update({mode : 'downreg' for mode in repressing_modes })
mode_map.update({mode : 'unknown' for mode in unknown_modes})
regulations['simple_mode'] = regulations['mode'].map(mode_map)

# extract values for names of genes and regulators
regulator_tag = regulations['regulator locus tag'].dropna().values
target_tag = regulations['locus tag'].values
regulator_name = regulations['regulator'].dropna().values
target_name = regulations['gene'].values

# also extract the names and tags of listed sporulation genes
gene_name = gene_list['locus tag'].values
gene_tag = gene_list['gene'].values


# encode genes with integer labels to nx doesnt act up
le = LabelEncoder()
total_labels = np.hstack((regulator_tag, target_tag, gene_tag))
total_enc = le.fit_transform(total_labels)
target_enc = le.transform(target_tag)
regulator_enc = le.transform(regulator_tag)
gene_enc = le.transform(gene_tag)


# Get ede attributes and assemble list
attrs = [{'mode' : m, 'regulon': r, 'simple_mode' : s} for m, r, s in zip(regulations['mode'], regulations['regulon'], regulations['simple_mode'])]
edge_list = [(u, v, d) for u, v, d in zip(regulator_enc, target_enc, attrs)]

G = nx.DiGraph()
G.add_edges_from(edge_list)

# get names and locus tags in the node attribute dict
node_attrs = {}
spor_set = set(gene_enc)
for node in G.nodes():
    tag = le.inverse_transform(np.array([node]))
    #name = regulations[regulations['locus tag'] == str(tag)]['gene'].values
    if node in spor_set:
        spor = True
    else:
        spor = False
    node_attrs.update({node : {'locus_tag': tag[0], 'sporulation': spor}})
    #node_attrs.update({str(node) : {'locus_tag': tag, 'gene': name, 'sporulation': spor}})

nx.set_node_attributes(G, node_attrs)

nx.write_gexf(G, 'subtiwiki_network.gexf')

