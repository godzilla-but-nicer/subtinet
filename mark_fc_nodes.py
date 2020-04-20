from glob import glob
from tqdm import tqdm
import pandas as pd
import networkx as nx

# threshold below which we don't save changes. can be zero
FC_THRESHOLD = 0
# get paths to RNA seq data files
DATA_DIR = 'data/significant_DESeq2_xlsx_files/'
data_files = glob(DATA_DIR + '*')

# later we'll need names for the network files
data_files_no_path = [fin.split('/')[-1] for fin in data_files]
treatments = [fin.split('_')[0] for fin in data_files_no_path]

# load existing network to write changes to
delta6 = nx.read_graphml('data/arrietta-orritz-delta6-genome.graphml')

# for each treatment we will add the fold changes under each treatment
# as attributes. We will track only signifigant changes. We could also
# exclude treatments that are too small, I'm not sure how we want to do it
for fin, tr in zip(data_files, treatments):
    print('\nTreatment:', tr)
    exdf = pd.read_csv(fin)

    # name of fold-change data col depends on treatment as will
    # attreibute name for the network
    fc_col = tr + '_IPTG_vs_' + tr + '_NT_log2fc'
    attr_name = tr + '_log2_fold_change'

    # we want to limit the columns that we track so we can neatly iterate
    # obviously iterating over columns isn't ideal but you know. We also
    # allow filtering of rows by magnitude of change.
    # NOTE that we are doing matching on the strain 168 names and had to drop
    # some values
    keep_cols = ['locus_tag.168', fc_col]
    changes = (exdf[(exdf[fc_col] > FC_THRESHOLD) |
                    (exdf[fc_col] < -FC_THRESHOLD)][keep_cols].dropna())

    # iterate over nodes in network and rows in changes df. If we find a match
    # in the d6 tags we change the fold change to the value found in exp.
    for node in tqdm(delta6.nodes(data=True)):
        hit_flag = False
        for _, tag, fc in changes.itertuples():
            # this is only needed for the 168 tags
            tag = tag.replace('_', '')
            if tag == node[1]['name']:
                node[1][attr_name] = fc
                hit_flag = True
        if hit_flag == False:
            node[1][attr_name] = 0.0

nx.write_graphml(
    delta6, 'data/arrietta-orritz-delta6-rnaseq-annotated.graphml')
