# Compute associations between deconvolution results (multideconv) and clinical metadata
#'
#' This function tests for associations between deconvolution results (multideconv)
#' and available clinical traits. It uses Pearson correlation for continuous (numeric) traits,
#' and ANOVA for categorical traits. Results are visualized as a labeled heatmap and violin plots.
#' All plots are saved in the `Results/` directory.
#'
#' deconv.results A numeric matrix or data frame of deconvolution matrix. Rows represent samples, columns represent celltypes modules.
#' coldata A data frame containing clinical traits (both categorical and numerical) for the same samples. Row names should match those of `deconv.results`.
#' @param pval A numeric threshold (default = 0.05) to determine statistical significance.
#'        Only associations with p-values below this threshold are considered significant in the heatmap.
#' @param width A numeric value indicating the width (in inches) of the output heatmap plot (default = 20).
#' @param height A numeric value indicating the height (in inches) of the output heatmap plot (default = 8).
#'
#' @return This function saves the following to the `Results/` directory:
#' A labeled heatmap showing Pearson correlations and ANOVA test p-values.
#' Individual violin plots for significant categorical trait associations.
#'
#' The function does not return an object to the R environment.

metada.association = function(deconv.results, coldata, pval = 0.05, width = 20, height = 8){
  ###Association with categorical variables
  coldata_categorical = coldata %>%
    dplyr::select(dplyr::where(is.character)|dplyr::where(is.factor))
  
  if(ncol(coldata_categorical) != 0){
    data = cbind(deconv.results, coldata_categorical)
    pvals = data.frame()
    fvals = data.frame()
    
    for(i in 1:ncol(deconv.results)){
      contador = 1
      for (j in (ncol(deconv.results)+1):ncol(data)) {
        module <- names(data)[i]
        trait <- names(data)[j]
        
        df <- data.frame(
          value = data[, i],
          trait = as.factor(data[, j])
        )
        
        # Global ANOVA (rstatix, so it's compatible downstream)
        res.aov <- rstatix::anova_test(data = df, dv = value, between = trait)
        
        pvals[i, contador] = res.aov$p
        fvals[i, contador] = res.aov$F
        contador = contador + 1
        
        # Only continue if ANOVA significant
        if(res.aov$p < pval) {
          
          # Pairwise Tukey
          pwc <- df %>%
            rstatix::tukey_hsd(value ~ trait) %>%
            rstatix::add_xy_position(x = "trait")
          
          pdf(paste0("Results/ANOVA_", module, "-", trait, ".pdf"))
          print(
            ggpubr::ggboxplot(df, x = "trait", y = "value", fill = "trait", add = "jitter") +
              ggpubr::stat_pvalue_manual(pwc, hide.ns = TRUE) +
              labs(
                title = "One-way ANOVA with Tukey HSD",
                subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
                caption  = rstatix::get_pwc_label(pwc),
                x = paste0("Clinical trait: ", trait),
                y = paste0("Values for ", module)
              ) +
              theme(
                axis.text.x = element_text(size = 20, angle = 0),
                axis.title.x = element_text(size = 20),
                legend.position = "none",
                axis.title.y = element_text(size = 20, angle = 90),
                plot.title = element_text(size = 20),
                plot.subtitle = element_text(size = 15, face = "bold")
              )
          )
          dev.off()
        }
      }
    }
    
    rownames(pvals) = colnames(deconv.results)
    colnames(pvals) = colnames(coldata_categorical)
    rownames(fvals) = colnames(deconv.results)
    colnames(fvals) = colnames(coldata_categorical)
    
    pvals = as.matrix(pvals)
    textMatrix2 = paste("ANOVA\n(", signif(pvals, 2), ")", sep = "")
    dim(textMatrix2) = dim(pvals)
  }
  
  ###Association with quantitative variables
  coldata_quantitative = coldata %>%
    dplyr::select(dplyr::where(is.numeric))
  
  if(ncol(coldata_quantitative)!=0){
    moduleTraitCor = WGCNA::cor(deconv.results, coldata_quantitative, method = "p", use = "pairwise.complete.obs")
    moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nrow(deconv.results))
    
    #### Replace pvalues for significance labels
    # Define cutoffs and corresponding labels
    breaks <- c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf)  # Define intervals
    labels <- c("****", "***", "**", "*", "")  # Corresponding significance labels
    
    # Replace p-values with significance levels
    moduleTraitPvalue <- matrix(
      cut(moduleTraitPvalue, breaks = breaks, labels = labels, right = FALSE),
      nrow = nrow(moduleTraitPvalue),
      dimnames = dimnames(moduleTraitPvalue)
    )
    
    
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(", moduleTraitPvalue, ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    
    if(ncol(coldata_categorical)!=0){
      textMatrix = cbind(textMatrix, textMatrix2)
      moduleTraitPvalue = cbind(moduleTraitPvalue, pvals)
      simulated_corr = matrix(stats::runif(n=nrow(pvals)*ncol(pvals), min=-0.1, max=0.1), nrow = nrow(pvals), ncol = ncol(pvals))
      colnames(simulated_corr) = colnames(pvals)
      moduleTraitCor = cbind(moduleTraitCor, simulated_corr)
    }
    
    idx = which(moduleTraitPvalue==""|moduleTraitPvalue>pval)
    for (i in idx) {
      textMatrix[i] = NA
    }
    
    moduleTraitCor_sig <- moduleTraitCor
    moduleTraitCor_sig[moduleTraitPvalue > pval] <- NA
    
    pdf("Results/deconv_metadata", width = width, height = height)
    par(mar = c(25, 15, 3, 3))
    WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                          xLabels = colnames(moduleTraitCor),
                          yLabels = rownames(moduleTraitCor),
                          ySymbols = rownames(moduleTraitCor),
                          colorLabels = FALSE,
                          colors = WGCNA::blueWhiteRed(50),
                          textMatrix = textMatrix,
                          setStdMargins = FALSE,
                          cex.text = 0.7,
                          zlim = c(-1,1),
                          main = paste0("Clinical associations with cell types\nOnly showing significant associations (pvalue < ", pval, ")"))
    dev.off()
    return(moduleTraitCor_sig)
  }
  
}



## Function to calculate the correlation between a cluster column and clinical data (categorical or continuous)
# Apply Fisher test for categorical and Anova for continuous

perform_statistical_tests <- function(data, cluster_col = "cluster") {
  
  # Separate the cluster column
  clusters <- data[[cluster_col]]
  variables <- data[, !names(data) %in% cluster_col]
  
  # Initialize the results dataframe
  results <- data.frame(
    variable = character(),
    type = character(),
    test = character(),
    statistic = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over each variable
  for (var_name in names(variables)) {
    var <- variables[[var_name]]
    
    # Determine if the variable is continuous or categorical
    is_continuous <- is.numeric(var) && length(unique(var)) > 10
    
    if (is_continuous) {
      # ANOVA test for continuous variables
      tryCatch({
        aov_result <- aov(var ~ clusters)
        aov_summary <- summary(aov_result)
        
        results <- rbind(results, data.frame(
          variable = var_name,
          type = "continuous",
          test = "ANOVA",
          statistic = aov_summary[[1]]$`F value`[1],
          p_value = aov_summary[[1]]$`Pr(>F)`[1]
        ))
      }, error = function(e) {
        message(paste("ANOVA error for", var_name, ":", e$message))
      })
      
    } else {
      # Fisher's exact test for categorical variables
      tryCatch({
        # Create contingency table
        cont_table <- table(var, clusters)
        
        # Fisher’s exact test (with simulation for large tables)
        test_result <- fisher.test(cont_table, simulate.p.value = TRUE)
        
        results <- rbind(results, data.frame(
          variable = var_name,
          type = "categorical",
          test = "Fisher exact",
          statistic = NA,  # Fisher test has no F or Chi-square statistic
          p_value = test_result$p.value
        ))
      }, error = function(e) {
        message(paste("Fisher test error for", var_name, ":", e$message))
      })
    }
  }
  
  # Add significance columns
  results$significant <- results$p_value < 0.05
  results$significance_stars <- cut(results$p_value,
                                    breaks = c(0, 0.001, 0.01, 0.05, 1),
                                    labels = c("***", "**", "*", "ns"),
                                    include.lowest = TRUE)
  
  return(results)
}


