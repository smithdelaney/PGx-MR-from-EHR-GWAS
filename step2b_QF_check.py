import pandas as pd
import numpy as np

#Step 1: Import the data
df = pd.read_csv('/Users/Delaney_Smith_1/Desktop/warfarin_by_genes_FDR.csv', sep='\t') #output from step 1 script
probes = pd.read_csv('/Users/Delaney_Smith_1/Desktop/warfarin_hits_step2a.csv', sep='\t') #ouput from step 2a script
probes = probes.rename(columns={'probes': 'ProbeID'})

# Step 2: Merge
merged_df = pd.merge(df, probes, on='ProbeID', how='outer')

#Step 3:calculate F statistic for weak instrument bias
ev = (merged_df['b_ivw']**2)/(merged_df['se_ivw']**2)
n = merged_df['n_iv']
rv = (merged_df['se_ivw']**2)
merged_df['F'] = (ev/n)/(rv/(n-1))

#Step 4: display which probes fail the F statistic test or show significant evidence of heterogeneity
print((merged_df['q_p'] < 0.05).sum()) 
print((merged_df['F'] < 10).sum())

#Step 5: remove probes that do not pass quality control and save file
merged_df = merged_df[~(merged_df['q_p'] < 0.05)]
merged_df = merged_df[~(merged_df['F'] < 10)]
merged_df.to_csv('/Users/Delaney_Smith_1/Desktop/warfarin_step2b.tsv', sep='\t', index=False) #output file name and path