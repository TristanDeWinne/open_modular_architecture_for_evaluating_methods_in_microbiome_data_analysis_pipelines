# This script contains a quick test implementation of DESeq2 to test the architecture
# while building it

# Load DESeq2 and install it if it isn't already
if(!require("DESeq2"))
{
  BiocManager::install("DESeq2")
  library("DESeq2")
}


# This function uses edgeR to analyze a dataset
DESeq2_test <- function(DESeq2_settings, sim_data)
{
  ##############################################################################
  # DESeq2_settings (list):
  #   A list of settings for the DESeq2 method. This can include the following:
  #     - design (string) the formula of the design matrix as a string.
  # 
  # sim_data (phyloseq):
  #   A phyloseq data object containing the simulated data.
  #
  #     For more information on DESeq2 see:
  #     https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
  ##############################################################################
  # Output (list): 
  #   A list containing the results, should at least have pvalues.
  ##############################################################################
  
  dds <- phyloseq_to_deseq2(physeq = sim_data,
                            design = as.formula(DESeq2_settings$design))
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  
  result <- results(dds)
  
  return(result)
}
