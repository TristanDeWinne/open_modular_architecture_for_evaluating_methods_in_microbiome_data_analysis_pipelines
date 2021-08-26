# This script contains the code for simulating data with the spsimseq method

# Load SPsimSeq and install it if it isn't already
if(!require("SPsimSeq"))
{
  library(BiocManager)
  BiocManager::install("SPsimSeq")
}

spsimseq_test <-function(sim_settings)
{
  ##############################################################################
  # sim_settings (list):
  #   A list containing the simulation settings. With following components:
  #     - source_data (string) a parameter giving the name of the data to load in
  #       this should be located in the folder source_data. It should at least 
  #       contain a count matrix. It should be saved in a RDS format.
  #     - n (integer) number of samples to generate. This is the same as tot.samples.
  #     - min_nonzero_sample (integer) Each gene requires a minimum number of samples.
  #       genes not fulfilling this requirement are filtered out.
  #     - n_genes (integer)[optional] the number of genes to retain from the source data.
  #
  #   Additional arguments see the bioconductor page of the SPsimSeq package
  #   http://www.bioconductor.org/packages/release/bioc/manuals/SPsimSeq/man/SPsimSeq.pdf
  #
  ##############################################################################
  # Output (phyloseq): 
  #   A phyloseq object containing the simulated data with SPsimSeq.
  ##############################################################################
  
  # Add .RDS extension if not done in the benchmark file
  source_data_name <- ifelse(".RDS" %in% sim_settings$source_data, 
                        sim_settings$source_data,
                        paste0(sim_settings$source_data, ".RDS"))
  
  input_data <- readRDS(paste0("source_data/", source_data_name))
  
  # Check source data for valid format
  if(!class(input_data) == "phyloseq")
  {
    stop(paste0("The source data should be a phyloseq object.", 
                "This is easily constructed with functions provided by the phyloseq package.", 
                "For more information see https://joey711.github.io/phyloseq/import-data.html"))
  }
  
  # Required format: Rows are otus and columns are samples
  if(!taxa_are_rows(input_data)) {input_data <- t(input_data)}
  

  
  # filter genes with sufficient expression (important step) 
  input_data <- filter_taxa(input_data,
                            function(x){sum(x>0) >= sim_settings$min_nonzero_sample},
                            prune = TRUE)
  

  # Subset source data to have n_genes genes if the parameter is defined
  if("n_genes" %in% names(sim_settings))
  {
    input_data <- prune_taxa(x = input_data,                            
                             taxa=taxa_names(input_data)[sample(nsamples(input_data), sim_settings$n_genes)])
  }
  
  # Define the group variable here to avoid clustering in the SPsimSeq main function
  group <- sample_data(input_data)[[sim_settings$group]]
  
  # tot.samples and n mean the same thing within the context of SPsimSeq.
  # if tot.samples is specified then use that, otherwise use n
  tot.samples <- ifelse(any("tot.samples" %in% names(sim_settings)),
                        sim_settings$tot.samples,
                        sim_settings$n)
  
  # Add the new variables to the settings list and then pass it to the general 
  # function run_w_args, this helps general function avoids hardcoding 
  # ifelse structures for each argument to add them
  sim_settings[["s.data"]] <- input_data
  sim_settings[["group"]] <- group
  sim_settings[["group.config"]] <- rep(1/length(unique(group)), length(unique(group))) # a balanced design
  sim_settings[["result.format"]] <- "list"
  
  # simulate data
  sim.data <- run_w_args("SPsimSeq", sim_settings)
  
  # Convert to phyloseq format
  sim.data_otu <- otu_table(sim.data[[1]]$counts, taxa_are_rows = TRUE)
  sim.data_sample <- sample_data(sim.data[[1]]$colData)
  sim.data_taxa <- tax_table(as.matrix(sim.data[[1]]$rowData))
  
  simulated_data <- phyloseq(sim.data_otu,
                             sim.data_sample,
                             sim.data_taxa)
  
  return(simulated_data)
}
