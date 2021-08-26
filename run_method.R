# This script contains the main function for calling the methods to benchmark
# And some additional functions used in the main program pretaining to setting 
# Up the benchmark

get_method_settings <- function(benchmark_settings)
{
  ##############################################################################
  # benchmark_settings (list):
  #   A list containing the settings of the benchmark i.e. the methods its 
  #   settings, the simulation settings and the metrics to use.
  ##############################################################################
  # Output (list): 
  #   A list containing only the settings concerning the methods.
  ##############################################################################
  
  # Each method should contain the string "method_setting" in the benchmark 
  # settings list
  method_settings <-  benchmark_settings[grepl(pattern = "method_setting", 
                                               x = names(benchmark_settings))]
  
  # If there are no method settings, throw a warning
  if(length(method_settings) == 0)
  {
    warning("The benchmark settings do not contain any method settings.")
  }
  
  # Return the method settings
  return(method_settings)
}

ground_truth <- function(data.phyloseq, truth_column_name = "DE.ind")
{
  ##############################################################################
  # data.phyloseq (phyloseq):
  #   A phyloseq data object containing a taxa table with the ground truth.
  # truth_column (string):
  #   A string with the column name containing the ground truth in the taxa table.
  ##############################################################################
  # Output (data.frame): 
  #   A data frame containing the ground truth wrt DE for the genes.
  ##############################################################################
  
  DE.id <- rank_names(data.phyloseq) == truth_column_name
  truth_column <- as.data.frame(tax_table(data.phyloseq)[,DE.id])
  
  return(truth_column)
}

# This is the main function for performing the analysis with the specified methods
# in the benchmark
perform_methods <- function(methods_settings=NULL, sim_datasets=NULL)
{
  ##############################################################################
  # methods_settings (list):
  #   A list with the methods and their settings.
  # sim_datasets (list):
  #   A list containing the simulated datasets as phyloseq objects.
  ##############################################################################
  # Output (list): 
  #   A list containing the results obtained with each method. 
  #   This should at least have a p-value for each feature 
  ##############################################################################
  
  # Check for valid input
  if(is.null(methods_settings) | is.null(sim_datasets))
  {
    stop("Method settings and/or simulated datasets are missing.")
  }
  
  # Initialize results collector object
  results.list <- list()
  
  # Perform analysis of simulated datasets for each specified method
  for(i in 1:length(methods_settings))
  {
    # Get the method settings
    analysis_method <- methods_settings[[i]]
    
    # Construct the name to save the results in, the iterator "i" is added as a 
    # suffix to make sure the name is unique
    method_results_name <- paste0(analysis_method$name, "_", i)
    
    # Source the analysis method. The name of the function should be contained
    # in "method_name.R"
    source(paste0("analysis_methods/", analysis_method$name, ".R"))
    
    # Perform the analysis method on the datasets iteratively
    for(j in 1:length(sim_datasets))
    {
      # Get the individual dataset to analyze
      sim_dataset <- sim_datasets[[j]]
      
      # The function to perform the method should be "analysis_method$name"
      # Construct the expression "analysis_method$name(methods_settings, sim_dataset)"
      analysis_method.expr <- parse(text = 
                                      paste0(analysis_method$name, 
                                             "(analysis_method, sim_dataset)"))
      
      # Run the method on the simulated dataset 
      # Analysis method method should output a list with results
      t_before <- Sys.time() # Used in runtime measurement
      results.list[[method_results_name]][[j]] <- eval(analysis_method.expr)
      t_after <- Sys.time() # Used in runtime measurement
      
      # Add runtime measurement to results
      results.list[[method_results_name]][[j]][["runtime"]] <- difftime(t_after,
                                                                        t_before, 
                                                                        units = "secs")
      # Add ground truth to the results, this is contained in the tax_table of 
      # the simulated data
      results.list[[method_results_name]][[j]][["ground_truth"]] <- ground_truth(sim_dataset)
    }
    
    # Add the method settings to the results list. This way the method and its 
    # settings are nicely linked with the results it generated
    results.list[[method_results_name]][["settings"]] <- analysis_method
  }
  
  # Return method results
  return(results.list)
}