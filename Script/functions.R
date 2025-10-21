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


## Heatmap function 
plot_heatmap<- function(relationship_result, 
                                           file_name = "features_traits_heatmap.pdf",
                                           width = 10, 
                                           height = 12,
                                           cluster_rows = FALSE,
                                           cluster_cols = FALSE) {
  # Extraire les données
  cor_matrix <- relationship_result[[1]]  # Matrice: traits (lignes) x features (colonnes)
  sig_list   <- relationship_result[[2]]  # Liste des features significatives par trait
  
  # Transposer pour avoir features en lignes et traits en colonnes
  cor_matrix_t <- t(cor_matrix)
  
  # Créer une matrice de texte vide pour afficher uniquement les corrélations significatives
  text_matrix <- matrix("", 
                        nrow = nrow(cor_matrix_t), 
                        ncol = ncol(cor_matrix_t),
                        dimnames = dimnames(cor_matrix_t))
  
  # Remplir UNIQUEMENT les corrélations significatives avec 2 décimales
  for (trait in colnames(cor_matrix_t)) {
    sig_features <- sig_list[[trait]]
    if (is.null(sig_features) || length(sig_features) == 0) next
    
    sig_features <- intersect(sig_features, rownames(cor_matrix_t))
    for (feat in sig_features) {
      val <- cor_matrix_t[feat, trait]
      text_matrix[feat, trait] <- sprintf("%.2f", val)
    }
  }
  
  # Extraire les types cellulaires et méthodes
  feature_names <- rownames(cor_matrix_t)
  
  # Fonction pour extraire le type cellulaire
  extract_cell_type <- function(name) {
    # Liste des méthodes connues
    methods <- c("DeconRNASeq", "Epidish", "DWLS", "Quantiseq", "EPIC", "xCell", "MCP", "TIMER")
    
    # Retirer la méthode au début SEULEMENT si elle existe
    name_clean <- name
    for (method in methods) {
      if (startsWith(name, paste0(method, "_"))) {
        name_clean <- sub(paste0("^", method, "_"), "", name)
        break
      }
    }
    
    # Si "Subgroup" ou "Iteration" présent, garder uniquement la partie avant
    if (grepl("_Subgroup|_Iteration", name_clean)) {
      cell_type <- sub("_Subgroup.*", "", name_clean)
      cell_type <- sub("_Iteration.*", "", cell_type)
    } else {
      # Sinon, comportement pour les noms techniques avec préfixes
      # Extraire après le dernier underscore si présent
      if (grepl("_", name_clean)) {
        cell_type <- sub(".*_", "", name_clean)
      } else {
        cell_type <- name_clean
      }
    }
    
    # Nettoyer les suffixes .Subgroup/.Iteration qui restent (sécurité)
    cell_type <- gsub("\\.Subgroup.*", "", cell_type)
    cell_type <- gsub("\\.Iteration.*", "", cell_type)
    
    # Remplacer les points par des espaces pour la lisibilité
    cell_type <- gsub("\\.", " ", cell_type)
    
    return(cell_type)
  }
  
  
  # Extraire les types cellulaires
  cell_types <- sapply(feature_names, extract_cell_type, USE.NAMES = FALSE)
  
  # Créer le dataframe d'annotation (sans la méthode)
  annotation_row <- data.frame(
    Cell_Types = cell_types,
    row.names = feature_names
  )
  
  # Définir les couleurs pour les types cellulaires
  unique_cells <- unique(cell_types)
  n_cells <- length(unique_cells)
  
  # Générer suffisamment de couleurs pour tous les types cellulaires
  if (n_cells <= 12) {
    # Utiliser Set3 qui a 12 couleurs max
    if (n_cells < 3) {
      cell_colors <- RColorBrewer::brewer.pal(3, "Set3")[1:n_cells]
    } else {
      cell_colors <- RColorBrewer::brewer.pal(n_cells, "Set3")
    }
  } else {
    # Utiliser rainbow pour un grand nombre de types
    cell_colors <- rainbow(n_cells, s = 0.6, v = 0.9)
  }
  names(cell_colors) <- unique_cells
  
  # Liste des couleurs d'annotation
  annotation_colors <- list(
    Type_Cellulaire = cell_colors
  )
  
  # Générer la heatmap
  pdf(file_name, width = width, height = height)
  pheatmap::pheatmap(
    cor_matrix_t,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    display_numbers = text_matrix,
    number_color = "black",
    fontsize_number = 8,
    na_col = "white",
    main = "Clinical associations with myeloid lineage cell types\nOnly showing significant associations (pvalue < 0.05)",
    angle_col = 45,
    fontsize_row = 8,
    fontsize_col = 10,
    border_color = NA,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors
  )
  dev.off()
  
  
  return(invisible(cor_matrix_t))
}


## Survival curve per group cells
compute.survival.by.features <- function(features, survival.data, PFS, PFS_event, time_unit, p.value = 0.05, thres = 0.5, max_factors = Inf) {
  library(survival)
  n_features <- ncol(features)
  significant_combinations <- list() # To store significant feature combinations
  # Generate all possible combinations of the features
  contador <- 1
  
  for (n in 1:min(n_features, max_factors)) {
    combinations <- combn(1:n_features, n, simplify = FALSE)
    
    for (comb in combinations) {
      # Create a formula dynamically based on the combination
      formula <- as.formula(paste("Surv(time, status) ~", paste(colnames(features)[comb], collapse = " + ")))
      
      # Prepare the data for survival analysis
      data_for_model <- data.frame("time" = survival.data[,PFS],
                                   "status" = survival.data[,PFS_event])
      quantiles <- quantile(features[,comb], thres)
      
      # Binarize the Cox model output to draw two KM lines (linear predictors are used to stratify between high-risk and low-risk groups)
      data_for_model$coxHL <- factor(ifelse(features[,comb] >= quantiles, paste0('High_', colnames(features)[comb]), paste0("Low_", colnames(features)[comb])))
      
      feature_name <- colnames(features)[comb]
      
      # Perform Kaplan-Meier based on coxHL
      km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)
      
      pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval
      
      if (!is.na(pval) && pval < p.value) {
        significant_combinations[[contador]] <- formula
        names(significant_combinations)[contador] <- paste0("Formula_", contador)
        
        ### Set colors of each curve with the correct order
        strata_names <- gsub("coxHL=", "", names(km_fit$strata))
        legend_labels <- strata_names
        legend.labs <- legend_labels
        
        pdf(paste0("Results/SurvPlot_", names(significant_combinations)[contador]), width = 10, height = 5, onefile = FALSE)
        print(survminer::ggsurvplot(km_fit,
                                    data = data_for_model,
                                    size = 1,
                                    palette = c("#E7B800", "#2E9FDF"),
                                    conf.int.style = "step",
                                    pval = TRUE,
                                    risk.table = TRUE,
                                    risk.table.col = "strata",
                                    legend.labs = legend.labs,
                                    risk.table.height = 0.3,
                                    ggtheme = theme_grey(),
                                    title = paste0("Cox PH for ", paste("Surv(time, status) ~", paste(colnames(features)[comb], collapse = " + "))),
                                    #subtitle = direction,
                                    xlab = paste0("Time to death/recurrence/progression (", time_unit, ")")
        ))
        dev.off()
        contador <- contador + 1
      }
    }
  }
  
  if (length(significant_combinations) == 0) {
    print("No significant combinations found.")
  } else {
    return(significant_combinations)
  }
}

## Wilcox test 
cell.groups.wilcox.test <- function(cell.groups, coldata, trait, pval = 0.05) {
  sig <- c()
  coldata[, trait] <- as.factor(coldata[, trait])
  
  if (length(unique(coldata[, trait])) != 2) {
    stop("Wilcoxon test requires a binary trait (exactly two groups).")
  }
  
  for (j in 1:ncol(cell.groups[[1]])) {
    data <- data.frame(Value = cell.groups[[1]][, j], Trait = coldata[, trait])
    
    res.test <- stats::wilcox.test(Value ~ Trait, data = data, exact = FALSE)
    
    if (round(res.test$p.value, 5) <= pval) {
      cat("Significant p-value after Wilcoxon test for", colnames(cell.groups[[1]])[j], "\n")
      
      # Save boxplot
      pdf(paste0("Results/Wilcoxon_", trait, "_", colnames(cell.groups[[1]])[j], ".pdf"),
          width = 12, height = 9)
      print(
        ggplot(data, aes(x = Trait, y = Value, fill = Trait)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "gray30") +
          geom_jitter(width = 0.15, size = 2.5, alpha = 0.8, color = "black") +
          stat_summary(fun = median, geom = "point", shape = 23,
                       size = 3, fill = "white") +
          scale_fill_brewer(palette = "Set2") +
          labs(
            title = paste0("Wilcoxon Test - ", colnames(cell.groups[[1]])[j]),
            subtitle = paste0("P-value: ", round(res.test$p.value, 5)),
            x = paste0("Clinical trait: ", trait),
            y = "Cell group score"
          ) +
          theme_minimal(base_size = 16) +
          theme(
            axis.text.x = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            axis.title.y = element_text(size = 20, angle = 90),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 15, face = "bold")
          ) +
          scale_fill_discrete(name = trait)
      )
      dev.off()
      
      sig <- c(sig, j)
    }
  }
  # Collect significant cell groups
  if (length(sig) == 0) {
    message("No significant cell groups (p-value < ", pval, ") after Wilcoxon test.")
    return(invisible(NULL))
  } else {
    cell.groups.sig <- list()
    cell.groups.sig[[1]] <- cell.groups[[1]][, sig, drop = FALSE]
    cell.groups.sig[[2]] <- cell.groups[[2]][sig]
    cell.groups.sig[[3]] <- cell.groups[[3]][sig]
    return(cell.groups.sig)
  }
}