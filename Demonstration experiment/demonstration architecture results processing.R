# This script contains the code to bundle the results of the archtitecture demonstration experiment

# Global parameters
alpha <- 0.05 # Control the fdr at the alpha level.

# Function to calculate the fdr based on the predicted values and the real values
calc_fdr <- function(truth, pred)
{
  truth_table <- table(Predictions = pred, 
                       Ground_truth = truth)
  TP <- truth_table[2,2]
  TN <- truth_table[1,1]
  FN <- truth_table[1,2]
  FP <- truth_table[2,1]
  
  # The fdr is the expected value of false discovery proportion
  fdr_val <- FP / (FP + TP)
  
  # Return fdr
  return(fdr_val)
}

# Function to calculate the sensitivity based on the predicted values and the real values
calc_sensitivity <- function(truth, pred)
{
  truth_table <- table(Predictions = pred, 
                       Ground_truth = truth)
  TP <- truth_table[2,2]
  TN <- truth_table[1,1]
  FN <- truth_table[1,2]
  FP <- truth_table[2,1]
  
  # The sensitivity is the true positive rate
  sensitivity_val <- TP / (TP + FN)
  
  # Return sensitivity
  return(sensitivity_val)
}

# Get all files containing the results
results_files <- list.files("benchmark_results/", pattern = "raw_results_demonstration")

# Loop over all files and merge them into one big data.frame
results.df <- data.frame("method" = character(),
                          "LFC" = numeric(),
                          "Fraction DE" = numeric(),
                          "n" = numeric(),
                          "sensitivity" = numeric(),
                          "FDR" = numeric(),
                          "runtime s" = numeric())
# Loop over each file
for(i in 1:length(results_files))
{
  # Read in single result file
  result_file <- results_files[i]
  result <- readRDS(paste0("benchmark_results/", result_file))
  
  # Loop over each method
  for(j in 1:(length(result)-1))
  {
    # Get the results for a single method
    results_method_j <- result[[j]]
    
    # Process each separate result
    for(k in 1:(length(results_method_j)-1))
    {
      results_method_j_sim_k <- results_method_j[[k]]
      
      # Adjust the pvalues with the BH95 method and handle each method separately
      if(results_method_j$settings$name == "DESeq2_test")
      {
        # Some pvalues are NA because a sample for that row contains an outlier 
        # (inner workings of DESeq2)
        adjusted_pvalues <- p.adjust(p = results_method_j_sim_k$pvalue,
                                     method = "fdr")
        method_name <- "DESeq2"
        trueDE_features <- results_method_j_sim_k$ground_truth[,1]
      }
      else if(results_method_j$settings$name == "edgeR_test")
      {
        adjusted_pvalues <- p.adjust(p = results_method_j_sim_k$table$PValue,
                                     method = "fdr")
        method_name <- "EdgeR"
        # EdgeR filters out some genes with insufficient expression strength but 
        # the ground truth still has all genes in it. Retain only the kept genes.
        trueDE_features <- results_method_j_sim_k$genes$DE.ind
      }
      else if(results_method_j$settings$name == "wilcoxon_rank_sum_test")
      {
        adjusted_pvalues <- p.adjust(p = results_method_j_sim_k$table$pvalue,
                                     method = "fdr")
        method_name <- "Wilcoxon rank sum test"
        trueDE_features <- results_method_j_sim_k$table$DE.ind
      }
      else
      {
        stop("Did you add a method to your result output but not add a way to retrieve its results?")
      }
      
      # The lines associated with an adjusted pvalue lower than alpha are 
      # classified as differentially abundant (filter out the NAs)
      DE_prediction <- (adjusted_pvalues < alpha)[!is.na(adjusted_pvalues)]
      DE_truth <- as.logical(trueDE_features[!is.na(adjusted_pvalues)])
      
      # Calculate senitivity and fdr
      # Sensitivity
      sensitivity <- calc_sensitivity(truth = DE_truth, pred = DE_prediction)
      # False discovery rate
      fdr <- calc_fdr(truth = DE_truth, pred = DE_prediction)
      
      # Get a single value for the runtime
      runtime <- as.numeric(unique(results_method_j_sim_k$runtime))
      
      # Construct the row to be added to the data.frame
      result_row <- data.frame("method" = method_name,
                               "LFC" = result$simulation_settings$lfc.thrld,
                               "Fraction DE" = result$simulation_settings$pDE,
                               "n" = result$simulation_settings$n,
                               "sensitivity" = sensitivity,
                               "FDR" = fdr,
                               "runtime s" = runtime)
      
      # Add row to data.frame
      results.df <- rbind(results.df, result_row)
    }
  }
}

# Save the data frame holding all the results
saveRDS(object = results.df,
        file ="benchmark_results/architecture_demonstration_results.RDS")