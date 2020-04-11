import numpy as np
import pandas as pd
import os
import glob
from tqdm import tqdm

RAW_DATA_DIR = 'data/DESeq2_output'
CLEANED_DATA_DIR = 'data/named_DESeq2_xlsx_files'

# later we will make our dataframe smaller
DROP_COLS_ALLOWED = ['Goe3_IPTG_r1', 'Goe3_IPTG_r2', 'Goe3_IPTG_r3', 'Goe3_NT_r1',
                     'Goe3_NT_r2', 'Goe3_NT_r3', 'pDR110_IPTG_r1', 'pDR110_IPTG_r2',
                     'pDR110_IPTG_r3', 'pDR110_NT_r1', 'pDR110_NT_r2', 'pDR110_NT_r3', 'sigF_IPTG_r1',
                     'sigF_IPTG_r2', 'sigF_IPTG_r3', 'sigF_NT_r1', 'sigF_NT_r2',
                     'sigF_NT_r3', 'sigG_IPTG_r1', 'sigG_IPTG_r2', 'sigG_IPTG_r3',
                     'sigG_NT_r1', 'sigG_NT_r2', 'sigG_NT_r3', 'SP10_IPTG_r1',
                     'SP10_IPTG_r2', 'SP10_IPTG_r3', 'SP10_NT_r1', 'SP10_NT_r2',
                     'SP10_NT_r3', 'Goe3_IPTG.x', 'Goe3_NT.x', 'pDR110_IPTG.x',
                     'pDR110_NT.x', 'sigF_IPTG.x', 'sigF_NT.x', 'sigG_IPTG.x', 'sigG_NT.x',
                     'SP10_IPTG.x', 'SP10_NT.x']

# We need to map the names used in the RNA seq data to the names used on the network
cds = pd.read_csv('data/delta6_168_cds_matched.csv')
tags = cds[['locus_tag.d6', 'locus_tag.168']]
map = {tag_d6: tag_168 for _, tag_d6, tag_168 in tags.itertuples()}

# sometimes our directory will be gone
if not os.path.exists(CLEANED_DATA_DIR):
    os.mkdir(CLEANED_DATA_DIR)

for dir in tqdm(os.listdir(RAW_DATA_DIR)):
    treat_dir = os.path.join(RAW_DATA_DIR, dir)
    data_file = glob.glob(treat_dir + '/*.xlsx')[0]

    # save less data, change names to match conventions
    df = pd.read_excel(data_file)
    df.columns = df.columns.to_series().apply(
        lambda x: x.strip())  # column names have spaces at end randomly

    drop_cols = [col for col in DROP_COLS_ALLOWED if col in df.columns]
    df = df.drop(drop_cols, axis=1)
    df['locus_tag.168'] = df['Name'].map(map)
    df['locus_tag.d6'] = df['Name']

    output_name = data_file.split('/')[-1]
    df.to_csv(CLEANED_DATA_DIR + '/' + output_name)
