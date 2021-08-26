# This R script contains additional functions which do not have a clear place 
# to belong to. 

# the following piece of code is taken from
# https://joey711.github.io/phyloseq-extensions/edgeR.html

################################################################################
#' Convert phyloseq OTU count data into DGEList for edgeR package
#' 
#' Further details.
#' 
#' @param physeq (Required).  A \code{\link{phyloseq-class}} or
#'  an \code{\link{otu_table-class}} object. 
#'  The latter is only appropriate if \code{group} argument is also a 
#'  vector or factor with length equal to \code{nsamples(physeq)}.
#'  
#' @param group (Required). A character vector or factor giving the experimental
#'  group/condition for each sample/library. Alternatively, you may provide
#'  the name of a sample variable. This name should be among the output of
#'  \code{sample_variables(physeq)}, in which case
#'  \code{get_variable(physeq, group)} would return either a character vector or factor.
#'  This is passed on to \code{\link[edgeR]{DGEList}},
#'  and you may find further details or examples in its documentation.
#'  
#' @param method (Optional). The label of the edgeR-implemented normalization to use.
#'  See \code{\link[edgeR]{calcNormFactors}} for supported options and details. 
#'  The default option is \code{"RLE"}, which is a scaling factor method 
#'  proposed by Anders and Huber (2010).
#'  At time of writing, the \link[edgeR]{edgeR} package supported 
#'  the following options to the \code{method} argument:
#'  
#'  \code{c("TMM", "RLE", "upperquartile", "none")}.
#'
#' @param ... Additional arguments passed on to \code{\link[edgeR]{DGEList}}
#' 
#' @examples
#' 
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}



################################################################################





# The code for the wilcoxon rank sum test is taken from 
# https://rdrr.io/github/TBrach/MicrobiomeX/

#######################################
### Wilcoxon Test
#######################################

#' Wilcoxon rank sum differential abundance test
#'
#' @description Perform Wilcoxon rank sum test on individual taxa in a phylseq object.
#'              If \code{block} is specified, then it is used as a blocking factor.
#'              If \code{block} is NULL, then standard Wilcoxon test is applied.
#'
#' @param physeq  Phyloseq object.
#' @param group   Name of column in sample metadata to perform group-wise comparisons on.
#' @param compare List of 2 groups to be compared, which are present in the \code{group} column.
#'                Can be ignored or set to \code{NULL} when there are only 2 levels in the factor.
#' @param block   Name of column in sample metadata to control for the group-wise comparisons.
#' @param excludeZeros logical. Should \code{abundance=0} for a taxon in a sample be ignored when performing Wilcoxon test?
#' @param p.adjust.method Method for adjusting P-values for multiple hypothesis testing.
#'
#' @return Returns a list containing the following:
#' \itemize{
#'   \item \code{table}: data frame containing the prevalence test results, with the following columns: \cr
#'             \code{Taxon, teststat, W, pvalue, Med_1, Med_2, padj, Significance, Direction, tax_table(phyloseq)}. \cr
#'             The teststatistic is based on the standardized teststatistic, equation provided by \code{multtest::mt.minP}.
#'   \item \code{test_name}: Name of test used - in this case \code{"Wilcoxon rank sum test"}.
#'   \item \code{block} : Name of variable used as blocking factor.
#' }
#'
#' @export
#'
#' @importFrom coin wilcox_test
#'
#' @examples
#' \dontrun{
#'
#' # If there are only 2 levels for \code{Diet}
#' res <- test_differential_abundance_Wilcoxon(physeq, group = "Diet", block = "Gender")
#'
#' # If there are multiple levels for \code{Diet}
#' res <- test_differential_abundance_Wilcoxon(physeq, group = "Diet", compare=c("Control", "High_Fat"), block = "Gender")
#' res <- test_differential_abundance_Wilcoxon(physeq, group = "Diet", compare=c("Control", "High_Fibre"), block = "Gender")
#' }
test_differential_abundance_Wilcoxon <-
  function(physeq,
           group,
           compare = NULL,
           block = NULL,
           excludeZeros = FALSE,
           p.adjust.method = "fdr") {
    
    if (taxa_are_rows(physeq)) {
      physeq <- t(physeq)
    }
    
    CT <- as(otu_table(physeq), "matrix")
    
    group_fac <- as.factor(sample_data(physeq)[[group]])
    
    if (!is.null(compare)) {
      group_var_levels <- compare
      subset <- group_fac %in% group_var_levels
    } else {
      group_var_levels <- levels(group_fac)
      subset = NULL
    }
    
    if (length(group_var_levels) > 2) {
      stop(paste0("test_differential_prevalence() can only generate results for comparing 2 groups - you asked for ", paste(group_var_levels, collapse=",")))
    }
    
    i <- group_var_levels[1]
    j <- group_var_levels[2]
    
    # Apply to every taxon
    #
    res_mat <- apply(CT, 2, function(taxon_counts) {
      x <- taxon_counts[group_fac == i]
      if (excludeZeros) {
        x <- x[x != 0]
      }
      Median_grp1 <- median(x, na.rm = T) # NA in case all 0
      
      y <- taxon_counts[group_fac == j]
      if (excludeZeros) {
        #if(all(y == 0)){y[1] <- ceiling(mean(taxon_counts))+1}
        y <- y[y != 0]
      }
      Median_grp2 <- median(y, na.rm = T)
      
      # Do not worry about the subsetting of x and y above under excludeZeros.
      # The data used by the actual test is untouched - it is the original data.
      
      data <- cbind("RA" = taxon_counts, sample_data(physeq))
      data[[group]] <- as.factor(data[[group]])
      
      formula <- paste0("RA ~ ", group)
      if (!is.null(block)) {
        formula <- paste0(formula, " | ", block)
        data[[block]] <- as.factor(data[[block]])
      }
      formula <- as.formula(formula)
      
      pval <- NA
      W <- NA
      standStat <- NA
      
      if (length(x) != 0 && length(y) != 0) {
        wt <-
          coin::wilcox_test(
            formula = formula,
            data = data,
            subset = subset,
            conf.int = FALSE,
            distribution = "asymptotic",
            alternative = "two.sided"
          )
        pval <- coin::pvalue(wt)
        W <- coin::statistic(wt)
        Ranks <- rank(c(x, y))
        n1 <- length(x)
        n2 <- length(y)
        standStat <-
          -1 * ((sum(Ranks[1:n1]) - n1 * (n1 + n2 + 1) / 2) / sqrt(n1 * n2 * (n1 +
                                                                                n2 + 1) / 12))
      }
      
      c(
        teststat = standStat,
        W = W,
        pvalue = pval,
        Median_grp1 = Median_grp1,
        Median_grp2 = Median_grp2
      )
    })
    
    res_mat <- as.data.frame(t(res_mat))
    padj <-
      p.adjust(res_mat$pval, method = p.adjust.method)
    significance <- sapply(padj, function(p) {
      sig = " NS"
      if (!is.nan(p) & !is.na(p)) {
        if (p <= 0.10) {
          sig = "."
        }
        if (p <= 0.05)  {
          sig = "*"
        }
        if (p <= 0.01)  {
          sig = "**"
        }
        if (p <= 0.001) {
          sig = "***"
        }
      }
      sig
    })
    
    direction <- rep(i, length(padj))
    direction[res_mat$teststat > 0] <- j
    
    DF <-
      data.frame(
        Taxon = as.character(rownames(res_mat)),
        res_mat,
        padj = padj,
        Significance = significance,
        Direction = direction,
        tax_table(physeq)
      )
    
    my_results <- vector(mode="list")
    my_results$table <- DF
    my_results$test_name <- "Wilcoxon rank sum test"
    my_results$block <- block
    
    return(my_results)
  }


# This is a general function to run a function with all possible arguments given 
# filled in. These arguments are contained in a list and should have the name of 
# the arguments of the function
run_w_args <- function(function_name, args.list)
{
  # get list of function arguments
  function_args <- formalArgs(function_name)
  
  # Run over each argument name and look if it is defined, if so then make an 
  # expression to add it to the eventual function call
  function_arg.str <- list()
  for(i in 1:length(function_args))
  {
    func_arg <- function_args[i]
    
    if(func_arg %in% names(args.list))
    {
      function_arg.str[[i]] <- paste0(func_arg, "=", paste0("args.list$", func_arg))
    }
  }
  
  # Construct function call
  function_arg.str <- unlist(function_arg.str)
  function_expression.str <- paste0(function_name,
                                    "(", 
                                    paste(function_arg.str, collapse = ","),
                                    ")")
  function_expression.expr <- parse(text = function_expression.str)
  
  # Return evaluation of expression of function call with the proper arguments filled in
  return(eval(function_expression.expr))
}



