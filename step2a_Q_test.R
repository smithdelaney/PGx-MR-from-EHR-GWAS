#R script to conduct cochran's Q test on hypothesis corrected probes form step 1

#Read in output from step1 as 'df' and .snps output file from SMR-IVW as 'snps'
df <- read.delim('/Users/Delaney_Smith_1/Desktop/warfarin_by_genes_FDR.csv',sep='\t' ) #probes from step 1
snps <- read.delim('/Users/Delaney_Smith_1/Desktop/EHR_PGx_MR/mqtl_ivw_warfarin_new_1000G.snps', sep='\t') #.snsps file

#select probeID
probes <- df$ProbeID

Qs <- list()
Q_ps <- list()

#Perform Cochran's Q-test for heterogeneity
for (probe in probes) {
  probe_df <- snps[snps$ProbeID == probe, ]
  se <- sqrt((probe_df$se.out^2) / (probe_df$beta.exp^2))
  b <- probe_df$beta.out / probe_df$beta.exp
  print(probe)
  res <- meta::metagen(b, se)
  Q_p  <- stats::pchisq(res$Q, res$df.Q, lower.tail=FALSE)
  Qs <- c(Qs, res$Q)
  Q_ps <- c(Q_ps, Q_p)
}

#Save results
probes <- cbind(probes, q_val = Qs)
probes <- cbind(probes, q_p = Q_ps)
write.table(probes, file='/Users/Delaney_Smith_1/Desktop/warfarin_hits_step2a.csv', sep='\t', row.names=FALSE) #set name to output file
head(df)