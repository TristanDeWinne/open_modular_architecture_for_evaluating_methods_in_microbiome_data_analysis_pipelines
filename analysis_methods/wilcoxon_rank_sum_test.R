# This script contains the implementation of the willcoxon rank sum test to do
# differential abundance testing

wilcoxon_rank_sum_test <- function(method_settings, sim_data)
{
  ##############################################################################
  # method_settings (list):
  #   A list containing the settings to be used by the wilcoxon rank sum test.
  #   The following are some arguments that can be passed:
  #     - group (string) name of grouping variable in the sample_data of the 
  #       phyloseq object. This group should either contain only 2 groups or 
  #       more, but then "compara" argument also needs to be specified.
  #     - compare (vector of strings) a vector containing 2 strings to specify 
  #       the subgroups to compare.
  #
  #     For more option to specify see:
  #     https://rdrr.io/github/TBrach/MicrobiomeX/man/test_differential_abundance_Wilcoxon.html
  ##############################################################################
  # Output (list): 
  #   A list containing the results, should at least have pvalues.
  ##############################################################################
  
  #  Normalize data via rarefying
  sim_data_normalized <- rarefy_even_depth(sim_data)
  
  # Retrieve optional parameters from method settigs
  if("compare" %in% names(method_settings)) 
    {compare <- method_settings$compare} 
  else 
    {compare <- NULL}
  
  # Perform  wilcoxon rank sum test
  results <- test_differential_abundance_Wilcoxon(
    physeq = sim_data_normalized,
    group = method_settings$group,
    compare = compare)
  
  # Return the results
  return(results)
}
