# This script contains a quick test implementation of edgeR to test the architecture
# while building it

# Load edgeR and install it if it isn't already
if(!require("edgeR"))
{
  BiocManager::install("edgeR")
  library("edgeR")
}


# This function uses edgeR to analyze a dataset
edgeR_test <- function(edgeR_settings, sim_data)
{
  ##############################################################################
  # edgeR_settings (list):
  #   A list of settings for the edgeR method. This can include the following:
  #     - group (string) the name of grouping variable.
  # 
  # sim_data (phyloseq):
  #   A phyloseq data object containing the simulated data.
  #
  #     For more information on edgeR see:
  #     https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  ##############################################################################
  # Output (list): 
  #   A list containing the results, should at least have pvalues.
  ##############################################################################
  
  DGE <- phyloseq_to_edgeR(sim_data, group=edgeR_settings$group)
  keep <- filterByExpr(DGE)
  DGE <- DGE[keep, , keep.lib.sizes=FALSE]
  DGE <- calcNormFactors(DGE)
  DGE <- estimateTagwiseDisp(estimateCommonDisp(DGE))
  
  result <- exactTest(DGE)
  
  return(result)
}
