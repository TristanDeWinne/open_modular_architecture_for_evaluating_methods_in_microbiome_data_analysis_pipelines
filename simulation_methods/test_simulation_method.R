# This is a test simulation script
# It will generate random normal draws

mat_list_to_phyloseq <- function(data.list_of_mat)
{
  # Constructing a phyloseq object based on a list of matrixes
  otu.matrix <- matrix()
  sample.df <- data.frame(exp_group = character())
  for(i in 1:length(data.list_of_mat))
  {
    if(i==1)
    {
      otu.matrix <- data.list_of_mat[[i]]
    }
    else
    {
      otu.matrix <- rbind(otu.matrix, data.list_of_mat[[i]])
    }
    sample.df <- rbind(sample.df,
                       data.frame(exp_group = 
                                    rep(paste0("exp_grp_", i),
                                        times = nrow(data.list_of_mat[[i]])))
    )
  }
  
  data.otu <- otu_table(otu.matrix, taxa_are_rows = FALSE)
  data.sample <- sample_data(sample.df)
  data.phyloseq  <-  phyloseq(data.otu, data.sample)
  
  return(data.phyloseq)
}

test_simulation_method <- function(settings.list)
{
  # Collector object of generated data
  test_simulated_data <- vector(mode = "list", length = 0)
  
  # Take random normal draws
  for(exp_group in 1:length(settings.list$n))
  {
    n_exp_group <- settings.list$n[exp_group]
    n_genes <- length(settings.list$mu[exp_group,])
    
    simulated_data_exp_group <- matrix(NA, nrow = n_exp_group, ncol = n_genes)
    colnames(simulated_data_exp_group) <- paste0("otu", 1:n_genes)
    
    for (gene in 1:length(settings.list$mu[exp_group,]))
    {
      mu_exp_grp_gene <- settings.list$mu[exp_group,gene]
      sd_exp_grp_gene <- settings.list$sd[exp_group,gene]
      
      simulated_data_exp_group[,gene] <- as.integer(rnorm(n_exp_group,
                                               mean = mu_exp_grp_gene,
                                               sd = sd_exp_grp_gene))
      
    }
    
    test_simulated_data[[paste0("exp_grp_", exp_group)]] <- simulated_data_exp_group
  }
  
  # Return simulated data as phyloseq object
  test_simulated_data.phylo <- mat_list_to_phyloseq(test_simulated_data)
  return(test_simulated_data.phylo)
}
