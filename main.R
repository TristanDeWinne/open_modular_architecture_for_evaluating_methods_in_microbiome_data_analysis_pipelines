# This is the main script of the benchmark architecture
# Test version 22/07/2021

# Load packages
library(phyloseq)
library(stringr)

# Set seed for reproducibility
set.seed(1234)

# Source all architecture scripts
source("read_settings.R")
source("simulate_data.R")
source("log_handling.R")
source("run_method.R")
source("helper_functions.R")

# Global parameters
save_data <- TRUE

# Get all benchmarks planned from folder benchmark/settings
benchmarks_planned <- list.files("benchmark_settings/")

# Perform each planned benchmark (this is the main loop)
for(i in 1:length(benchmarks_planned))
{
  # Benchmark i file location
  benchmark_file <- benchmarks_planned[i]
  benchmark_log_file <- str_replace(benchmark_file, "\\.txt", "_log\\.txt")
  write_to_log("Start benchmark.", benchmark_log_file)
  
  # Read settings contained in the benchmark file
  write_to_log("Read in settings (start).", benchmark_log_file)
  benchmark_setting <- read_settings(paste0("benchmark_settings/",
                                            benchmark_file))
  write_to_log("Read in settings (end).", benchmark_log_file)
  
  ##                            SIMULATION STEP                               ##
  
  # Simulate data according to the benchmark setting
  write_to_log("Simulation of data (start).", benchmark_log_file)
  simulated_data <- simulate_data(benchmark_setting$sim_setting)
  write_to_log("Simulation of data (end).", benchmark_log_file)
  
  # If save_data is set to TRUE then save the simulated data 
  # (overwriting existing files if equal names)
  if(save_data)
  {
    write_to_log("Save simulated data (start).", benchmark_log_file)
    saveRDS(simulated_data,
            file = paste0("temp/sim_data_", 
                          str_replace(benchmark_file, "\\.txt", ".RDS")))
    write_to_log("Save simulated data (end).", benchmark_log_file)
  }
  
  
  ##                            METHOD STEP                                   ##
  
  # Retrieve all method settings from the benchmark settings list
  methods_settings <- get_method_settings(benchmark_setting)
  
  # Run the methods on the simulated data and add simulation settings to the results
  write_to_log("Run methods (start).", benchmark_log_file)
  benchmark_results <- perform_methods(methods_settings = methods_settings,
                                       sim_datasets = simulated_data)
  benchmark_results[["simulation_settings"]] <- benchmark_setting[["sim_setting"]]
  write_to_log("Run methods (end).", benchmark_log_file)

  # Save the results from the benchmark
  write_to_log("Save benchmark results in folder benchmark_results (start).", 
               benchmark_log_file)
  saveRDS(object = benchmark_results, 
          file = paste0("benchmark_results/", 
                        "raw_results_", 
                        str_replace(benchmark_file, "\\.txt", ".RDS")))
  write_to_log("Save benchmark results in folder benchmark_results (end).", 
               benchmark_log_file)
  
  
  # Add new line to log file to mark end of the benchmark run
  write_to_log("End benchmark.\n", benchmark_log_file)
}


