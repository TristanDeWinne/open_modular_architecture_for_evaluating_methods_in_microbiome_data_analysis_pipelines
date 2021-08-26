# This script contains the general function to perform the simulation
# It will source the appropriate .R file to perform the simulation and call its
# functions to simulate N datasets

simulate_data <- function(simulation_settings=NULL)
{
  ##############################################################################
  # simulation_settings (list):
  #   A list containing the simulation settings. Mandatory main component 
  #   "sim_setting" with sub components:
  #     -"name" the name of the simulation method to call. This should be the 
  #       name of the script and the function within the script should also be 
  #       the same. e.g. the method SPsimSeq is contained in spsimseq_test.R
  #       the function that will be called has as name "spsimseq_test"
  #     - "N" the number of datasets to generate.
  #     - "n" the number(s) of samples to generate for each dataset.
  #
  #     - additonal settings possible
  ##############################################################################
  # Output (list): 
  #   A list containing the simulated datasets as phyloseq objects.
  ##############################################################################
  
  # Check input list for mandatory components
  if(!all(c("name", "N", "n") %in% names(simulation_settings)))
  {stop("Simulation setting does not contain a name/N/n component.")}
  
  # Initialize dataset collector list
  simulated_datasets <-  vector(mode = "list", length = simulation_settings$N)
  
  # Simulate N datasets
  for(i in 1:simulation_settings$N)
  {
    simulation_method <- simulation_settings$name
    
    # Source the simulation method. This should be contained in 
    # "simulation_name.R"
    source(paste0("simulation_methods/", simulation_method, ".R"))
    
    # The function to perform the simulation should be "simulation_name"
    # Construct the expression "simulation_name(simulation_settings)"
    simulation_method.expr <- parse(text = paste0(simulation_method, 
                                                  "(simulation_settings)"))
    
    # Simulation method should output a phyloseq object
    simulated_datasets[[i]] <- eval(simulation_method.expr)
  }
  
  # Return the simulated datasets as a list of phyloseq datasets 
  return(simulated_datasets)
}














