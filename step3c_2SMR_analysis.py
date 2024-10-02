import rpy2.robjects as robjects
import pandas as pd
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import DataFrame
import sys
import statsmodels.stats.multitest as mt

#load R object from step 3b
robjects.r['load']('/Users/Delaney_Smith_1/Desktop/warfarin_2SMR_results/warfarin_results.Rdata') #update name to your R data object
loaded_env = robjects.r['comb']
data = dict(zip(loaded_env.names, list(loaded_env)))

mr_sig = {}

#optional: check all MR method p-vals (Egger, IVW, simple median)
for key in data: 
    df_list = data[key]
    df_mr = df_list[0]
    df_mr = pandas2ri.rpy2py_dataframe(df_mr)
    df_clean = df_mr.iloc[:-2]
    all_below_threshold = (df_clean['pval'] < 0.05/len(data)).all()
    
    if all_below_threshold:
            mr_sig[key] = data[key]

print(len(mr_sig.keys()))
print(mr_sig.keys())

q_sig = {}

#optional: second check for heterogeneity
for key in data:
    #df_list = mr_sig[key]
    df_list = data[key]
    df_q = df_list[1]
    df_q = pandas2ri.rpy2py_dataframe(df_q)
    all_above_threshold = (df_q['Q_pval'] > 0.05/len(mr_sig)).all()
    if all_above_threshold:
            #q_sig[key] = mr_sig[key]
            q_sig[key] = data[key]

print(len(q_sig.keys()))

pl_sig = {}

#check for directional pleitropy
for key in q_sig: 
    df_list = q_sig[key]
    df_pl = df_list[3]
    df_pl = pandas2ri.rpy2py_dataframe(df_pl)
    all_above_threshold = (df_pl['pval.pl'] > 0.05/len(q_sig)).all()
    if all_above_threshold:
            pl_sig[key] = q_sig[key]
            #print (key)
            #print(pl_sig[key])

output = pd.DataFrame()
input = pd.read_csv('/Users/Delaney_Smith_1/Desktop/warfarin_step2b.tsv', sep = '\t') #file generated in step 2b

for key in pl_sig:
      row = input[input['ProbeID'] == key]
      output = pd.concat([output, row], ignore_index=True)
      
output.to_csv("/Users/Delaney_Smith_1/Desktop/warfarin_step3c.csv", sep = '\t') #filepath to output file

print(pl_sig.keys())