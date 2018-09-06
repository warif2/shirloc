# Title     : sleuth_summarize.R
# Objective : consolidate outputs from all fractions into single file
# Created by: leweezard
# Created on: 9/4/18

# Read in argument listing all paths of comparisons to perform.
args <- commandArgs(trailingOnly = TRUE)

# Store path arguments into variables main_path and sl_out_path.
main_path <- args[1]
sl_out_path <- unlist(strsplit(args[2], ':'))

# Metric values of interest to consolidate
metrics <- c('b','se_b','pval','qval')

# Combine values from all fractions into single tables for output.
for (met in metrics){

    # Get label for fractional comparison for first path in index
    frac_comp <- tail(unlist(strsplit(sl_out_path[1], '/')), n = 1)

    # Open fractional comparison and store in df
    sl_df <- read.csv(paste0(sl_out_path[1],'/sleuth_output.csv'))

    # Subset df for values of interest
    met_df <- subset(sl_df, select = c('target_id', met))
    names(met_df)[names(met_df) == met] <- frac_comp

    # Open sleuth outputs and store information into respective tables.
    for (i in 2:length(sl_out_path)) {

        # Get label for fractional comparison
        frac_comp <- tail(unlist(strsplit(sl_out_path[i], '/')), n = 1)

        # Open fractional comparison and store in df
        sl_df <- read.csv(paste0(sl_out_path[i],'/sleuth_output.csv'))

        # Subset df for values of interest
        sl_subset <- subset(sl_df, select = c('target_id', met))
        met_df <- merge(met_df, sl_subset, by = 'target_id')
        names(met_df)[names(met_df) == met] <- frac_comp
    }

    # Save file into summary directory
    out_name = paste0(main_path,'/summary/',met,'.csv')
    write.csv(met_df, file = out_name, quote = FALSE, row.names = FALSE)
}

for (met in metrics){

    # Open consolidated file store in df
    combine_df <- read.csv(paste0(main_path,'/summary/',met,'.csv'))

    # Subset df for values of interest
    met_df <- subset(sl_df, select = c('target_id', met))
    names(met_df)[names(met_df) == met] <- frac_comp

    # Open sleuth outputs and store information into respective tables.
    for (i in 2:length(sl_out_path)) {

        # Get label for fractional comparison
        frac_comp <- tail(unlist(strsplit(sl_out_path[i], '/')), n = 1)

        # Open fractional comparison and store in df
        sl_df <- read.csv(paste0(sl_out_path[i],'/sleuth_output.csv'))

        # Subset df for values of interest
        sl_subset <- subset(sl_df, select = c('target_id', met))
        met_df <- merge(met_df, sl_subset, by = 'target_id')
        names(met_df)[names(met_df) == met] <- frac_comp
    }

    # Save file into summary directory
    out_name = paste0(main_path,'/summary/',met,'.csv')
    write.csv(met_df, file = out_name, quote = FALSE, row.names = FALSE)
}