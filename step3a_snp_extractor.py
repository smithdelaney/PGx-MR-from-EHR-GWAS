import pandas as pd

#Read output from step 2b in as 'df' and .snps output from SMR-IVW in as 'snps'
df = pd.read_csv('/Users/Delaney_Smith_1/Desktop/warfarin_step2b.tsv', sep='\t') #output from step 2b
snps = pd.read_csv('/Users/Delaney_Smith_1/Desktop/EHR_PGx_MR/mqtl_ivw_warfarin_new_1000G.snps', sep='\t') #.snps file

# 'probes' is a list of probe IDs from step 2b
probes = df['ProbeID'].tolist()

#create new dataframe
all_snps = pd.DataFrame()

#extract all SNPs assosciated with probes from step 2b
for probe in probes:
    probe_df = snps[snps['ProbeID'] == probe]
    all_snps = pd.concat([all_snps, probe_df])

#print output and save
print(all_snps)
all_snps.to_csv('/Users/Delaney_Smith_1/Desktop/warfarin_step3a.csv', sep = '\t', index=False) #output file name