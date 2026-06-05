library(Seurat)
library(harmony)
library(mclust)
library(aricode)
library(ggplot2)
library(dplyr)
library(patchwork)
library(openxlsx)
library(broom)

setwd("/Users/martinacavallini/Library/CloudStorage/Box-Box/DevelopmentalProteomics/Submission/NatureNeuroscience/UploadtoGithub/Analysis code")

### Clustering performance for CBCs---------------------------------------------
## IPL-CSP gene sets---------------

CBC_obj <- readRDS("./Saved_objects/CBC_object_allgenes.rds") #Load datasets
CBC_subset_IPLCSP <- readRDS("./Saved_objects/CBC_object_IPLCSP.rds")

IPL_CSP <- read.xlsx("./Data/IPL_CSP_list.xlsx") #Load IPL-CSP genes
genes_vector_IPLCSP <- unlist(strsplit(as.character(IPL_CSP$Gene.Names), split = " "))
genes_vector_IPLCSP <- unique(trimws(genes_vector_IPLCSP))

TP_list <- read.xlsx("./Data/Uniprot_TPlist_20250718.xlsx") #Load surface protein list from Uniprot
genes_vector_TPlist <- unlist(strsplit(as.character(TP_list$Gene.Names), split = " "))
genes_vector_TPlist <- unique(trimws(genes_vector_TPlist))
available_genes <- genes_vector_TPlist[genes_vector_TPlist %in% rownames(CBC_obj)]

#Remove IPL-CSP genes from random sample gene sets
availablegenes_notcsp <- available_genes[!(available_genes %in% genes_vector_IPLCSP)]
#availablegenes_notcsp <- available_genes[!(available_genes %in% genes_vector_Mods)]


resolutions <- c(0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
n_random_samples <- 30

results_df <- data.frame() #Initialize results
seurat_objects <- list() #Initialize Seurat object list

#Function to calculate ARI and NMI for different resolutions
calculate_metrics_for_resolutions <- function(seu_obj, sample_name) {
  ari_values <- c()
  nmi_values <- c()
  
  for (res in resolutions) {
    seu_obj <- FindClusters(seu_obj, resolution = res, verbose = FALSE)
    cluster_labels <- seu_obj$seurat_clusters
    reference_labels <- seu_obj$cell_type
    ari_value <- ARI(cluster_labels, reference_labels)
    nmi_value <- NMI(cluster_labels, reference_labels)
    ari_values <- c(ari_values, ari_value)
    nmi_values <- c(nmi_values, nmi_value)
  }
  
  data.frame(
    resolution = resolutions,
    ari = ari_values,
    nmi = nmi_values,
    sample = sample_name,
    line_type = ifelse(sample_name == "csp_genes", "solid", "dotted")
  )
}

#ARI, NMI for IPL-CSP genes
csp_results <- calculate_metrics_for_resolutions(CBC_subset_IPLCSP, "csp_genes")
results_df <- rbind(results_df, csp_results)
seurat_objects[["csp_genes"]] <- CBC_subset_IPLCSP

#ARI, NMI for random gene sets
for (i in 1:n_random_samples) {
  cat("Processing random sample", i, "of", n_random_samples, "\n")
  
  set.seed(i)
  genes_vector <- sample(availablegenes_notcsp, 1169)
  
  subset_seu <- subset(CBC_obj, features = genes_vector)
  subset_seu <- NormalizeData(subset_seu, verbose = FALSE)
  subset_seu <- FindVariableFeatures(subset_seu, verbose = FALSE)
  subset_seu <- ScaleData(subset_seu, verbose = FALSE)
  subset_seu <- RunPCA(subset_seu, npcs = 50, verbose = FALSE)
  subset_seu <- RunHarmony(subset_seu, group.by.vars = "orig.ident", verbose = FALSE)
  subset_seu <- RunUMAP(subset_seu, reduction = "harmony", dims = 1:10, verbose = FALSE)
  subset_seu <- FindNeighbors(subset_seu, reduction = "harmony", dims = 1:10, verbose = FALSE)
  
  sample_name <- paste0("random_", i)
  sample_results <- calculate_metrics_for_resolutions(subset_seu, sample_name)
  results_df <- rbind(results_df, sample_results)
  seurat_objects[[sample_name]] <- subset_seu
  
  cat("Added results for", sample_name, "\n")
}


statistical_analysis <- function(results_df) {
  
  final_results <- data.frame()
  
  for (res in resolutions) {
    csp_ari <- results_df[results_df$sample == "csp_genes" & results_df$resolution == res, "ari"]
    csp_nmi <- results_df[results_df$sample == "csp_genes" & results_df$resolution == res, "nmi"]
    
    random_ari <- results_df[results_df$sample != "csp_genes" & results_df$resolution == res, "ari"]
    random_nmi <- results_df[results_df$sample != "csp_genes" & results_df$resolution == res, "nmi"]
    
    # One-sample t-tests
    ari_test <- t.test(random_ari, mu = csp_ari, alternative = "two.sided")
    nmi_test <- t.test(random_nmi, mu = csp_nmi, alternative = "two.sided")
    
    # Effect sizes
    ari_cohens_d <- (csp_ari - mean(random_ari)) / sd(random_ari)
    nmi_cohens_d <- (csp_nmi - mean(random_nmi)) / sd(random_nmi)
    
    # 95% Confidence intervals for the difference
    ari_diff <- csp_ari - mean(random_ari)
    nmi_diff <- csp_nmi - mean(random_nmi)
    
    ari_se <- sd(random_ari) / sqrt(length(random_ari))
    nmi_se <- sd(random_nmi) / sqrt(length(random_nmi))
    
    ari_ci_lower <- ari_diff - 1.96 * ari_se
    ari_ci_upper <- ari_diff + 1.96 * ari_se
    nmi_ci_lower <- nmi_diff - 1.96 * nmi_se
    nmi_ci_upper <- nmi_diff + 1.96 * nmi_se
    
    cat("Resolution", res, ":\n")
    cat("  ARI: CSP =", round(csp_ari, 4), 
        "vs Random =", round(mean(random_ari), 4), "±", round(sd(random_ari), 4), "\n")
    cat("       Difference =", round(ari_diff, 4), 
        "95% CI: [", round(ari_ci_lower, 4), ",", round(ari_ci_upper, 4), "]\n")
    cat("       Cohen's d =", round(ari_cohens_d, 3), "| p =", format.pval(ari_test$p.value), "\n")
    
    cat("  NMI: CSP =", round(csp_nmi, 4), 
        "vs Random =", round(mean(random_nmi), 4), "±", round(sd(random_nmi), 4), "\n")
    cat("       Difference =", round(nmi_diff, 4), 
        "95% CI: [", round(nmi_ci_lower, 4), ",", round(nmi_ci_upper, 4), "]\n")
    cat("       Cohen's d =", round(nmi_cohens_d, 3), "| p =", format.pval(nmi_test$p.value), "\n\n")
    
    final_results <- rbind(final_results, data.frame(
      resolution = res,
      csp_ari = csp_ari,
      random_ari_mean = mean(random_ari),
      ari_difference = ari_diff,
      ari_cohens_d = ari_cohens_d,
      ari_pvalue = ari_test$p.value,
      ari_ci_lower = ari_ci_lower,
      ari_ci_upper = ari_ci_upper,
      csp_nmi = csp_nmi,
      random_nmi_mean = mean(random_nmi),
      nmi_difference = nmi_diff,
      nmi_cohens_d = nmi_cohens_d,
      nmi_pvalue = nmi_test$p.value,
      nmi_ci_lower = nmi_ci_lower,
      nmi_ci_upper = nmi_ci_upper
    ))
  }
  
  # Multiple testing correction
  final_results$ari_pvalue_adj <- p.adjust(final_results$ari_pvalue, method = "BH")
  final_results$nmi_pvalue_adj <- p.adjust(final_results$nmi_pvalue, method = "BH")
  
  correction_table <- final_results[, c("resolution", "ari_pvalue", "ari_pvalue_adj", 
                                        "nmi_pvalue", "nmi_pvalue_adj")]
  print(correction_table)
  
  return(final_results)
}

final_results <- statistical_analysis(results_df)


#Plot results
random_summary <- results_df %>%
  filter(sample != "csp_genes") %>%
  group_by(resolution) %>%
  summarise(
    mean_ari = mean(ari),
    sd_ari = sd(ari),
    se_ari = sd(ari) / sqrt(n()),
    ci_lower_ari = mean_ari - 1.96 * se_ari,
    ci_upper_ari = mean_ari + 1.96 * se_ari,
    mean_nmi = mean(nmi),
    sd_nmi = sd(nmi),
    se_nmi = sd(nmi) / sqrt(n()),
    ci_lower_nmi = mean_nmi - 1.96 * se_nmi,
    ci_upper_nmi = mean_nmi + 1.96 * se_nmi,
    .groups = 'drop'
  )

csp_data <- results_df %>% filter(sample == "csp_genes")

add_significance_stars <- function(p_value) {
  if (p_value < 0.001) return("***")
  else if (p_value < 0.01) return("**")
  else if (p_value < 0.05) return("*")
  else if (p_value < 0.1) return("†")
  else return("ns")
}

random_summary$ari_significance <- sapply(final_results$ari_pvalue_adj, add_significance_stars)
random_summary$nmi_significance <- sapply(final_results$nmi_pvalue_adj, add_significance_stars)

p_ari <- ggplot() +
  # Random samples confidence interval (95% CI)
  geom_ribbon(data = random_summary, 
              aes(x = resolution, ymin = ci_lower_ari, ymax = ci_upper_ari),
              alpha = 0.3, fill = "deepskyblue") +
  # Random samples mean ± SD
  geom_ribbon(data = random_summary, 
              aes(x = resolution, ymin = mean_ari - sd_ari, ymax = mean_ari + sd_ari),
              alpha = 0.2, fill = "lightgray") +
  # Individual random samples (thin lines)
  geom_line(data = filter(results_df, sample != "csp_genes"), 
            aes(x = resolution, y = ari, group = sample), 
            color = "gray", alpha = 0.5, size = 0.5) +
  # Random samples mean line
  geom_line(data = random_summary, 
            aes(x = resolution, y = mean_ari), 
            color = "black", size = 1, linetype = "dashed") +
  # CSP genes line
  geom_line(data = csp_data, 
            aes(x = resolution, y = ari), 
            color = "red", size = 2) +
  geom_point(data = csp_data, 
             aes(x = resolution, y = ari), 
             color = "red", size = 4) +
  # Add significance annotations
  geom_text(data = random_summary, 
            aes(x = resolution, y = max(random_summary$ci_upper_ari) + 0.02, 
                label = ari_significance),
            size = 4, vjust = 0) +
  labs(x = "Resolution", y = "ARI Value", 
       title = "CSP Genes vs Random Genes: ARI Comparison",
       subtitle = "Blue = Random genes (mean ± 95% CI), Red = CSP genes\n*p<0.05, **p<0.01, ***p<0.001, †p<0.1, ns=not significant") +
  ylim(0, 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

p_nmi <- ggplot() +
  # Random samples confidence interval (95% CI)
  geom_ribbon(data = random_summary, 
              aes(x = resolution, ymin = ci_lower_nmi, ymax = ci_upper_nmi),
              alpha = 0.3, fill = "deepskyblue") +
  # Random samples mean ± SD
  geom_ribbon(data = random_summary, 
              aes(x = resolution, ymin = mean_nmi - sd_nmi, ymax = mean_nmi + sd_nmi),
              alpha = 0.2, fill = "lightgray") +
  # Individual random samples
  geom_line(data = filter(results_df, sample != "csp_genes"), 
            aes(x = resolution, y = nmi, group = sample), 
            color = "gray", alpha = 0.5, size = 0.5) +
  # Random samples mean line
  geom_line(data = random_summary, 
            aes(x = resolution, y = mean_nmi), 
            color = "black", size = 1, linetype = "dashed") +
  # CSP genes line
  geom_line(data = csp_data, 
            aes(x = resolution, y = nmi), 
            color = "red", size = 2) +
  geom_point(data = csp_data, 
             aes(x = resolution, y = nmi), 
             color = "red", size = 4) +
  # Add significance annotations
  geom_text(data = random_summary, 
            aes(x = resolution, y = max(random_summary$ci_upper_nmi) + 0.02, 
                label = nmi_significance),
            size = 4, vjust = 0) +
  labs(x = "Resolution", y = "NMI Value", 
       title = "CSP Genes vs Random Genes: NMI Comparison",
       subtitle = "Blue = Random genes (mean ± 95% CI), Red = CSP genes") +
  ylim(0, 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

combined_plot <- p_ari / p_nmi
print(combined_plot) 


#Optional: create UMAP plots to visualize clustering with random gene sets
seu_obj <- seurat_objects[["random_8"]] #change to any of random_1 through random_30
DimPlot(seu_obj, reduction = "umap", group.by = "cell_type", label = T, repel = T, label.size = 4) + NoLegend()

#Optional: save list of seurat objects
saveRDS(seurat_objects, "./Saved_objects/CBC_seurat_objects_IPLCSP_sets.rds", compress = "gzip")


## Module gene sets---------------

CBC_obj <- readRDS("./Saved_objects/CBC_object_allgenes.rds") #Load datasets
CBC_subset_Mods <- readRDS("./Saved_objects/CBC_object_Mods.rds")

CBC_mods <- read.xlsx("./Data/CBC_modules.xlsx") #Load CBC module genes
genes_vector_Mods <- CBC_mods$gene_name

TP_list <- read.xlsx("./Data/Uniprot_TPlist_20250718.xlsx") #Load surface protein list from Uniprot
genes_vector_TPlist <- unlist(strsplit(as.character(TP_list$Gene.Names), split = " "))
genes_vector_TPlist <- unique(trimws(genes_vector_TPlist))
available_genes <- genes_vector_TPlist[genes_vector_TPlist %in% rownames(CBC_obj)]

#Remove module genes from random sample gene sets
availablegenes_notcsp <- available_genes[!(available_genes %in% genes_vector_Mods)]


resolutions <- c(0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
n_random_samples <- 30

results_df <- data.frame() #Initialize results
seurat_objects <- list() #Initialize Seurat object list

#Function to calculate ARI and NMI for different resolutions
calculate_metrics_for_resolutions <- function(seu_obj, sample_name) {
  ari_values <- c()
  nmi_values <- c()
  
  for (res in resolutions) {
    seu_obj <- FindClusters(seu_obj, resolution = res, verbose = FALSE)
    cluster_labels <- seu_obj$seurat_clusters
    reference_labels <- seu_obj$cell_type
    ari_value <- ARI(cluster_labels, reference_labels)
    nmi_value <- NMI(cluster_labels, reference_labels)
    ari_values <- c(ari_values, ari_value)
    nmi_values <- c(nmi_values, nmi_value)
  }
  
  data.frame(
    resolution = resolutions,
    ari = ari_values,
    nmi = nmi_values,
    sample = sample_name,
    line_type = ifelse(sample_name == "csp_genes", "solid", "dotted")
  )
}

#ARI, NMI for IPL-CSP genes
csp_results <- calculate_metrics_for_resolutions(CBC_subset_IPLCSP, "csp_genes")
results_df <- rbind(results_df, csp_results)
seurat_objects[["csp_genes"]] <- CBC_subset_IPLCSP

#ARI, NMI for random gene sets
for (i in 1:n_random_samples) {
  cat("Processing random sample", i, "of", n_random_samples, "\n")
  
  set.seed(i)
  genes_vector <- sample(availablegenes_notcsp, 248)
  
  subset_seu <- subset(CBC_obj, features = genes_vector)
  subset_seu <- NormalizeData(subset_seu, verbose = FALSE)
  subset_seu <- FindVariableFeatures(subset_seu, verbose = FALSE)
  subset_seu <- ScaleData(subset_seu, verbose = FALSE)
  subset_seu <- RunPCA(subset_seu, npcs = 50, verbose = FALSE)
  subset_seu <- RunHarmony(subset_seu, group.by.vars = "orig.ident", verbose = FALSE)
  subset_seu <- RunUMAP(subset_seu, reduction = "harmony", dims = 1:10, verbose = FALSE)
  subset_seu <- FindNeighbors(subset_seu, reduction = "harmony", dims = 1:10, verbose = FALSE)
  
  sample_name <- paste0("random_", i)
  sample_results <- calculate_metrics_for_resolutions(subset_seu, sample_name)
  results_df <- rbind(results_df, sample_results)
  seurat_objects[[sample_name]] <- subset_seu
  
  cat("Added results for", sample_name, "\n")
}


statistical_analysis <- function(results_df) {
  
  final_results <- data.frame()
  
  for (res in resolutions) {
    csp_ari <- results_df[results_df$sample == "csp_genes" & results_df$resolution == res, "ari"]
    csp_nmi <- results_df[results_df$sample == "csp_genes" & results_df$resolution == res, "nmi"]
    
    random_ari <- results_df[results_df$sample != "csp_genes" & results_df$resolution == res, "ari"]
    random_nmi <- results_df[results_df$sample != "csp_genes" & results_df$resolution == res, "nmi"]
    
    # One-sample t-tests
    ari_test <- t.test(random_ari, mu = csp_ari, alternative = "two.sided")
    nmi_test <- t.test(random_nmi, mu = csp_nmi, alternative = "two.sided")
    
    # Effect sizes
    ari_cohens_d <- (csp_ari - mean(random_ari)) / sd(random_ari)
    nmi_cohens_d <- (csp_nmi - mean(random_nmi)) / sd(random_nmi)
    
    # 95% Confidence intervals for the difference
    ari_diff <- csp_ari - mean(random_ari)
    nmi_diff <- csp_nmi - mean(random_nmi)
    
    ari_se <- sd(random_ari) / sqrt(length(random_ari))
    nmi_se <- sd(random_nmi) / sqrt(length(random_nmi))
    
    ari_ci_lower <- ari_diff - 1.96 * ari_se
    ari_ci_upper <- ari_diff + 1.96 * ari_se
    nmi_ci_lower <- nmi_diff - 1.96 * nmi_se
    nmi_ci_upper <- nmi_diff + 1.96 * nmi_se
    
    cat("Resolution", res, ":\n")
    cat("  ARI: CSP =", round(csp_ari, 4), 
        "vs Random =", round(mean(random_ari), 4), "±", round(sd(random_ari), 4), "\n")
    cat("       Difference =", round(ari_diff, 4), 
        "95% CI: [", round(ari_ci_lower, 4), ",", round(ari_ci_upper, 4), "]\n")
    cat("       Cohen's d =", round(ari_cohens_d, 3), "| p =", format.pval(ari_test$p.value), "\n")
    
    cat("  NMI: CSP =", round(csp_nmi, 4), 
        "vs Random =", round(mean(random_nmi), 4), "±", round(sd(random_nmi), 4), "\n")
    cat("       Difference =", round(nmi_diff, 4), 
        "95% CI: [", round(nmi_ci_lower, 4), ",", round(nmi_ci_upper, 4), "]\n")
    cat("       Cohen's d =", round(nmi_cohens_d, 3), "| p =", format.pval(nmi_test$p.value), "\n\n")
    
    final_results <- rbind(final_results, data.frame(
      resolution = res,
      csp_ari = csp_ari,
      random_ari_mean = mean(random_ari),
      ari_difference = ari_diff,
      ari_cohens_d = ari_cohens_d,
      ari_pvalue = ari_test$p.value,
      ari_ci_lower = ari_ci_lower,
      ari_ci_upper = ari_ci_upper,
      csp_nmi = csp_nmi,
      random_nmi_mean = mean(random_nmi),
      nmi_difference = nmi_diff,
      nmi_cohens_d = nmi_cohens_d,
      nmi_pvalue = nmi_test$p.value,
      nmi_ci_lower = nmi_ci_lower,
      nmi_ci_upper = nmi_ci_upper
    ))
  }
  
  # Multiple testing correction
  final_results$ari_pvalue_adj <- p.adjust(final_results$ari_pvalue, method = "BH")
  final_results$nmi_pvalue_adj <- p.adjust(final_results$nmi_pvalue, method = "BH")
  
  correction_table <- final_results[, c("resolution", "ari_pvalue", "ari_pvalue_adj", 
                                        "nmi_pvalue", "nmi_pvalue_adj")]
  print(correction_table)
  
  return(final_results)
}

final_results <- statistical_analysis(results_df)


#Plot results
random_summary <- results_df %>%
  filter(sample != "csp_genes") %>%
  group_by(resolution) %>%
  summarise(
    mean_ari = mean(ari),
    sd_ari = sd(ari),
    se_ari = sd(ari) / sqrt(n()),
    ci_lower_ari = mean_ari - 1.96 * se_ari,
    ci_upper_ari = mean_ari + 1.96 * se_ari,
    mean_nmi = mean(nmi),
    sd_nmi = sd(nmi),
    se_nmi = sd(nmi) / sqrt(n()),
    ci_lower_nmi = mean_nmi - 1.96 * se_nmi,
    ci_upper_nmi = mean_nmi + 1.96 * se_nmi,
    .groups = 'drop'
  )

csp_data <- results_df %>% filter(sample == "csp_genes")

add_significance_stars <- function(p_value) {
  if (p_value < 0.001) return("***")
  else if (p_value < 0.01) return("**")
  else if (p_value < 0.05) return("*")
  else if (p_value < 0.1) return("†")
  else return("ns")
}

random_summary$ari_significance <- sapply(final_results$ari_pvalue_adj, add_significance_stars)
random_summary$nmi_significance <- sapply(final_results$nmi_pvalue_adj, add_significance_stars)

p_ari <- ggplot() +
  # Random samples confidence interval (95% CI)
  geom_ribbon(data = random_summary, 
              aes(x = resolution, ymin = ci_lower_ari, ymax = ci_upper_ari),
              alpha = 0.3, fill = "deepskyblue") +
  # Random samples mean ± SD
  geom_ribbon(data = random_summary, 
              aes(x = resolution, ymin = mean_ari - sd_ari, ymax = mean_ari + sd_ari),
              alpha = 0.2, fill = "lightgray") +
  # Individual random samples (thin lines)
  geom_line(data = filter(results_df, sample != "csp_genes"), 
            aes(x = resolution, y = ari, group = sample), 
            color = "gray", alpha = 0.5, size = 0.5) +
  # Random samples mean line
  geom_line(data = random_summary, 
            aes(x = resolution, y = mean_ari), 
            color = "black", size = 1, linetype = "dashed") +
  # CSP genes line
  geom_line(data = csp_data, 
            aes(x = resolution, y = ari), 
            color = "red", size = 2) +
  geom_point(data = csp_data, 
             aes(x = resolution, y = ari), 
             color = "red", size = 4) +
  # Add significance annotations
  geom_text(data = random_summary, 
            aes(x = resolution, y = max(random_summary$ci_upper_ari) + 0.02, 
                label = ari_significance),
            size = 4, vjust = 0) +
  labs(x = "Resolution", y = "ARI Value", 
       title = "Module Genes vs Random Genes: ARI Comparison",
       subtitle = "Blue = Random genes (mean ± 95% CI), Red = CSP genes\n*p<0.05, **p<0.01, ***p<0.001, †p<0.1, ns=not significant") +
  ylim(0, 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

p_nmi <- ggplot() +
  # Random samples confidence interval (95% CI)
  geom_ribbon(data = random_summary, 
              aes(x = resolution, ymin = ci_lower_nmi, ymax = ci_upper_nmi),
              alpha = 0.3, fill = "deepskyblue") +
  # Random samples mean ± SD
  geom_ribbon(data = random_summary, 
              aes(x = resolution, ymin = mean_nmi - sd_nmi, ymax = mean_nmi + sd_nmi),
              alpha = 0.2, fill = "lightgray") +
  # Individual random samples
  geom_line(data = filter(results_df, sample != "csp_genes"), 
            aes(x = resolution, y = nmi, group = sample), 
            color = "gray", alpha = 0.5, size = 0.5) +
  # Random samples mean line
  geom_line(data = random_summary, 
            aes(x = resolution, y = mean_nmi), 
            color = "black", size = 1, linetype = "dashed") +
  # CSP genes line
  geom_line(data = csp_data, 
            aes(x = resolution, y = nmi), 
            color = "red", size = 2) +
  geom_point(data = csp_data, 
             aes(x = resolution, y = nmi), 
             color = "red", size = 4) +
  # Add significance annotations
  geom_text(data = random_summary, 
            aes(x = resolution, y = max(random_summary$ci_upper_nmi) + 0.02, 
                label = nmi_significance),
            size = 4, vjust = 0) +
  labs(x = "Resolution", y = "NMI Value", 
       title = "Module Genes vs Random Genes: NMI Comparison",
       subtitle = "Blue = Random genes (mean ± 95% CI), Red = CSP genes") +
  ylim(0, 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

combined_plot <- p_ari / p_nmi
print(combined_plot) 


#Optional: create UMAP plots to visualize clustering with random gene sets
seu_obj <- seurat_objects[["random_8"]] #change to any of random_1 through random_30
DimPlot(seu_obj, reduction = "umap", group.by = "cell_type", label = T, repel = T, label.size = 4) + NoLegend()

#Optional: save list of seurat objects
saveRDS(seurat_objects, "./Saved_objects/CBC_seurat_objects_Module_sets.rds", compress = "gzip")
