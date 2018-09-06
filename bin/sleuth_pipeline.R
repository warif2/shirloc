# Title     : sleuth_pipeline.R
# Objective : quantify transcript abundance differences between two polysome profiling rna-seq fractions
# Created by: leweezard
# Created on: 9/2/18

# Import sleuth library
library(sleuth)

# Read in argument listing all paths of comparisons to perform
analysis_path <- commandArgs(trailingOnly = TRUE)

# Open metadata file stored in paths
s2c <- read.csv(paste(analysis_path,'/sleuth_metadata.txt', sep = ""), sep = ',')
s2c$path <- as.character(s2c$path)

# Filter Function
mod_filter <- function(row, min_reads = 5, min_prop = 0.47) {
  mean(row >= min_reads) >= min_prop
}

# Initialize sleuth object
so <- sleuth_prep(s2c, read_bootstrap = TRUE, extra_bootstrap_summary = TRUE, filter_fun = mod_filter, transformation_function = function(x) log2(x + 0.5))

# Fit LRT 'full' model
so <- sleuth_fit(so, ~fraction, 'full')

# Fit LRT 'reduced' model
so <- sleuth_fit(so, ~1, 'reduced')

# Perform the LRT
so <- sleuth_lrt(so, 'reduced', 'full')

# Perfrom the Wald test
so <- sleuth_wt(so, paste0('fraction'))

# Generate results table
sleuth_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_wt <- sleuth_results(so, 'fraction', show_all = TRUE)
res <- merge(sleuth_lrt, sleuth_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], on = 'target_id', sort = FALSE)

# Sort result table
res <- res[order(res[,1]),]

# Extract normalized counts
norm_counts <- kallisto_table(so, normalized = TRUE, include_covariates = TRUE)

# Write results into output file
write.csv(res, file = paste(analysis_path,'/sleuth_output.csv', sep = ""), quote = FALSE, row.names = FALSE)
write.csv(norm_counts, file = paste0(analysis_path,'/normalized_counts.csv'), quote = FALSE, row.names = FALSE)