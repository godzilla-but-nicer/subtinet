import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# load graph
G = nx.read_graphml(snakemake.input.gml)

# extract degrees
degree = [k[1] for k in G.degree()]
nonzero = [k for k in degree if k > 0]
num_zeros = len([z for z in degree if z == 0])
sorted_degree = sorted(nonzero)

# report some stuff. Might need to note number of zeros somewhere
print('(Min, Median, Max) :', (min(degree), np.median(degree), max(degree)))
print('Number of nodes: {}'.format(len(degree)))
print('Number of nonzero degree nodes: {0} ({1:.2f})'.format(len(nonzero), len(nonzero) / len(degree)))

# calculate bins
n_bins = 20
bins = np.logspace(np.log10(sorted_degree[0]), np.log10(sorted_degree[-1]), n_bins, base=10)
counts, bin_edges = np.histogram(sorted_degree, bins)
norm_counts = np.array(counts) / np.diff(np.array(bin_edges))

# count CCDF
ccdf = np.cumsum(norm_counts[::-1])
rev_bins = bin_edges[::-1]

# draw plot
fig, ax = plt.subplots()
ax.scatter(rev_bins[1:], ccdf)
ax.set_xscale('log')
ax.set_yscale('log')
# ax.set_ylim((10**-3, 10**5))
# ax.set_xlim((10**-0.1, 10**34))
ax.set_xlabel(r'$\log{k}$')
ax.set_ylabel(r'$\log{count}$')
plt.savefig(snakemake.output[0])