#R script to perform 2SMR tests on probes from step 3a

#Run these commands the first time you use the script
#install.packages("meta")
#library(meta)

#run these commands each time, unless saved to environment
library(TwoSampleMR)
library(ggplot2)

#makes directory, update name before running
dir_path <- "/Users/Delaney_Smith_1/Desktop/warfarin_2SMR_results" 
dir.create(dir_path)

#the step 3a csv file needs a new header as shown below (manually change before running code)
#'gene	SNP	effect_allele	other_allele	beta	se	beta.out	se.out'
#read in the file and extract probes
df <- read.delim('/Users/Delaney_Smith_1/Desktop/warfarin_step3a.csv', sep='\t') #output of step 3a here
p <- df$gene
probes <- unique(p)

comb <- list()

#perform test for each probe
for (probe in probes) {
  
  temp <- list()

  probe_path <- paste0("/Users/Delaney_Smith_1/Desktop/warfarin_2SMR_results/" , probe) #set output directory
  dir.create(probe_path)
  
  probe_df <- df[df$gene == probe, ]
  exp_dat <- format_data(probe_df, type = "exposure") #read in exposure data probe by probe

  out_dat <- read_outcome_data( #read in outcome data for each probe
    snps = exp_dat$SNP,
    filename = '/Users/Delaney_Smith_1/Desktop/warfarin_step3a.csv', #hit snps file again
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta.out",
    se_col = "se.out",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  
  #harmonize the exposure and outcome data
  dat <- harmonise_data( 
    exposure_dat = exp_dat, 
    outcome_dat = out_dat,
    action = 1
  )

  #perform mr egger regression, mr simple median, mr weighted median, mr ivw
  res <- mr(dat, method_list = c("mr_egger_regression_bootstrap", "mr_simple_median", "mr_weighted_median", "mr_ivw", "mr_sign"))
  new_res <- res[, c("method", "b", "se", "pval")]
  names(new_res)[names(new_res) == "method"] <- "MR_method"
  temp <- append(temp, list(new_res))
  
  #heterogeneity testing (again) via cochran's Q
  het <- mr_heterogeneity(dat)
  new_het <- het[, c("method", "Q", "Q_df", "Q_pval")]
  names(new_het)[names(new_het) == "method"] <- "Het_method"
  temp <- append(temp, list(new_het))
  new_het <- het[, c("method", "Q", "Q_df", "Q_pval")]
  
  #egger regression test for directional pleitropy
  pl <- mr_pleiotropy_test(dat) 
  new_pl <- pl[, c("egger_intercept", "se", "pval")]
  names <- c("egger_intercept", "se.pl", "pval.pl")
  names(new_pl) <- names
  temp <- append(temp, list(new_pl))
  print(temp)

  #visualize results using a scatter plot
  scatter <- mr_scatter_plot(res, dat)
  scatter_name = paste0(probe_path, "/", probe, "_scatter.png")
  ggsave(scatter[[1]], file = scatter_name, width = 7, height = 7)

  #save results to overarching dataframe  
  comb[[probe]] <- temp
  
}

#save all results for all probes
save(comb, file = "/Users/Delaney_Smith_1/Desktop/warfarin_2SMR_results/warfarin_results.Rdata") #for use in python, save object in directory



