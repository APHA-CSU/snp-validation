import pandas as pd

snps = pd.read_csv('blah.tab', delimiter='\t')

mutations = pd.read_csv('AF2122.refseq2simseq.map.txt', delimiter='\t')

mutated_pos = set(mutations['ref_start'].values)
snps_pos = set(snps['POS'].values)

missed = mutated_pos - snps_pos


print(mutations)