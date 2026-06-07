# STEP 1
BiocManager::install("DESeq2", force = TRUE)
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("BiocManager")
install.packages("dplyr")
install.packages("pak")
pak::pak("antpiron/RedRibbon")
install.packages("remotes")
remotes::install_github("antpiron/RedRibbon", force = TRUE)
"scales" %in% loadedNamespaces()
install.packages("pheatmap")
install.packages("data.table")
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("AnnotationDbi")
BiocManager::install("fgsea")
install.packages("msigdbr")
BiocManager::install("sva")

library(DESeq2)
library(sva)
library(data.table)
library(ggplot2)
library(ggrepel)

hc_count <- fread("HC_counts_an.txt", header = TRUE)
md_count <- fread("MD_counts_an.txt", header = TRUE)
metadata_full <- read.csv("Metadata_rin_tin_age_sex.csv") 
metadata_hc_md <- subset(metadata_full, 
                            group %in% c("HC", "MD") & 
                              region_acr == "CgGr")
hc_columns <- grep("CgGr", colnames(hc_count), ignore.case = TRUE)
hc_count_cggr <- hc_count[, ..hc_columns]
hc_counts_cggr <- cbind(hc_count[, 1, with = FALSE], hc_count_cggr)
md_columns <- grep("CgGr", colnames(md_count), ignore.case = TRUE)
md_count_cggr <- md_count[, ..md_columns]
md_counts_cggr <- cbind(md_count[, 1, with = FALSE], md_count_cggr)

metadata_hc_md$sex <- factor(metadata_hc_md$sex)
metadata_hc_md$group <- factor(metadata_hc_md$group)
metadata_hc_md$age <- as.numeric(metadata_hc_md$age)
metadata_hc_md$rin <- as.numeric(metadata_hc_md$rin)
metadata_hc_md$tin <- as.numeric(metadata_hc_md$tin)
metadata_hc_md$sex <- relevel(metadata_hc_md$sex, ref = "m")
metadata_hc_md$group <- relevel(metadata_hc_md$group, ref = "HC")

counts_all <- cbind(hc_count_cggr, md_count_cggr)
counts_all_with_genes <- cbind(hc_count[, 1, with = FALSE], counts_all)
rownames(counts_all_with_genes) <- counts_all_with_genes$gene_id
counts_all_with_genes$gene_id <- NULL
counts_matrix <- as.matrix(counts_all_with_genes)
rownames(counts_matrix) <- rownames(counts_all_with_genes)

metadata_hc_md$age_scaled <- scale(metadata_hc_md$age, center = TRUE, scale = TRUE)
metadata_hc_md$rin_scaled <- scale(metadata_hc_md$rin, center = TRUE, scale = TRUE)
metadata_hc_md$tin_scaled <- scale(metadata_hc_md$tin, center = TRUE, scale = TRUE)

p1 <- ggplot(metadata_hc_md, aes(x = group, y = rin, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "RIN distribution by diagnosis",
       y = "RIN value") +
  theme_minimal() +
  theme(legend.position = "none")
print(p1)

p2 <- ggplot(metadata_hc_md, aes(x = group, y = tin, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "TIN distribution by diagnosis",
       y = "TIN value") +
  theme_minimal() +
  theme(legend.position = "none")
print(p2)

cor_rin_tin <- cor(metadata_hc_md$rin, metadata_hc_md$tin)
cat("Pearson correlation between RIN and TIN:", round(cor_rin_tin, 3), "\n")

p3 <- ggplot(metadata_hc_md, aes(x = rin, y = tin, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = paste0("RIN vs TIN correlation (r = ", round(cor_rin_tin, 3), ")"),
       x = "RIN", y = "TIN") +
  theme_minimal()
print(p3)

t_test_rin <- t.test(rin ~ group, data = metadata_hc_md)
print(t_test_rin)
t_test_tin <- t.test(tin ~ group, data = metadata_hc_md)
print(t_test_tin)

print(t_test_rin$p.value)
print(t_test_tin$p.value)

identify_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  return(x < lower_bound | x > upper_bound)
}

metadata_hc_md$rin_outlier <- identify_outliers(metadata_hc_md$rin)
metadata_hc_md$tin_outlier <- identify_outliers(metadata_hc_md$tin)

cat("RIN outliers:", sum(metadata_hc_md$rin_outlier), "\n")
cat("TIN outliers:", sum(metadata_hc_md$tin_outlier), "\n")

tin_outliers <- metadata_hc_md[identify_outliers(metadata_hc_md$tin), , drop = FALSE]
tin_outliers[, c("group", "sex", "rin", "tin")]

table(RIN_low = metadata_hc_md$rin < 6, Group = metadata_hc_md$group)

dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = metadata_hc_md,
                              design = ~ tin_scaled + age_scaled + sex + group + sex:group)

keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep,]

dat  <- counts(dds, normalized = FALSE)
mod  <- model.matrix(~ tin_scaled + age_scaled + sex * group, data = colData(dds))
mod0 <- model.matrix(~ 1, data = colData(dds)) 

svseq <- svaseq(dat, mod, mod0)

if (svseq$n.sv > 0) {
  cat("Number of hidden batch factors:", svseq$n.sv, "\n")
  
  sv_matrix <- svseq$sv
  colnames(sv_matrix) <- paste0("SV", 1:svseq$n.sv)
  colData(dds) <- cbind(colData(dds), sv_matrix)
  
  sv_terms <- paste(colnames(sv_matrix), collapse = " + ")
  dynamic_formula <- as.formula(paste("~", sv_terms, "+ tin_scaled + age_scaled + sex * group"))
  
  design(dds) <- dynamic_formula
  cat("New design formula:", paste(deparse(dynamic_formula), collapse = ""), "\n")
} else {
  cat("No hidden batch effects were detected in the matrix. We're sticking with the standard design.\n")
  design(dds) <- ~ tin_scaled + age_scaled + sex * group
}

dds <- DESeq(dds)

# STEP 2

rld <- rlog(dds, blind = TRUE)
pca_data <- assay(rld)
pca_matrix <- t(pca_data)

pca <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
summary(pca)
variance_explained <- summary(pca)$importance[2, ] * 100

plot(pca$x[,1], pca$x[,2], 
     xlab = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
     main = "PCA: All samples")

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3],
  Sample = colnames(pca_data),
  Group = colData(dds)$group,
  Sex = colData(dds)$sex,
  Age = colData(dds)$age,
  RIN = colData(dds)$rin,
  TIN = colData(dds)$tin
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Group), linetype = "dashed") +
  labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       title = "PCA: Colored by Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(pca_df, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Sex), linetype = "dashed") +
  scale_color_manual(values = c("f" = "#F8766D", "m" = "#00BFC4")) +
  labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       title = "PCA: Colored by Sex") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Sex)) +
  geom_point(size = 4) +
  labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       title = "PCA: Color = Group, Shape = Sex") +
  theme_minimal() +
  theme(legend.position = "bottom")

library(patchwork)

p1 <- ggplot(pca_df, aes(x = Age, y = PC1, color = Group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "PC1 vs Age") +
  theme_minimal()

p2 <- ggplot(pca_df, aes(x = RIN, y = PC1, color = Group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "PC1 vs RIN") +
  theme_minimal()

p3 <- ggplot(pca_df, aes(x = TIN, y = PC1, color = Group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "PC1 vs TIN") +
  theme_minimal()

combined_plot <- p1 / (p2 | p3)
combined_plot

pc_correlations <- data.frame(
  Variable = c("Age", "RIN", "TIN"),
  Correlation = c(
    cor(pca_df$PC1, pca_df$Age),
    cor(pca_df$PC1, pca_df$RIN),
    cor(pca_df$PC1, pca_df$TIN)
  )
)

ggplot(pc_correlations, aes(x = Variable, y = Correlation, fill = Variable)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  labs(title = "Correlation of covariates with PC1",
       y = "Pearson correlation") +
  theme_minimal()

pdf("EDA_Report.pdf", width = 10, height = 8)

print(ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 4) +
        stat_ellipse(aes(group = Group)) +
        labs(title = "PCA: Group"))
print(ggplot(pca_df, aes(x = PC1, y = PC2, color = Sex)) +
        geom_point(size = 4) +
        stat_ellipse(aes(group = Sex)) +
        labs(title = "PCA: Sex"))
print(ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Sex)) +
        geom_point(size = 4) +
        labs(title = "PCA: Combined"))

var_df <- data.frame(PC = 1:10, 
                     Variance = summary(pca)$importance[2, 1:10] * 100)
print(ggplot(var_df, aes(x = PC, y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_line() +
        geom_point() +
        labs(title = "Scree plot: Variance explained by each PC") + 
        theme_minimal())

print(ggplot(pca_df, aes(x = Age, y = PC1, color = Group)) +
        geom_point() + geom_smooth(method = "lm"))

print(ggplot(pca_df, aes(x = RIN, y = PC1, color = Group)) +
        geom_point() + geom_smooth(method = "lm"))

print(ggplot(pca_df, aes(x = TIN, y = PC1, color = Group)) +
        geom_point() + geom_smooth(method = "lm"))

dev.off()

interpret_pca <- function(pca_df, variance_explained) {
  
  group_sep <- abs(mean(pca_df$PC1[pca_df$Group == "MD"]) - 
                         mean(pca_df$PC1[pca_df$Group == "HC"]))
  
  sex_sep <- abs(mean(pca_df$PC1[pca_df$Sex == "m"]) - 
                   mean(pca_df$PC1[pca_df$Sex == "f"]))
  
  cat("Total variance explained by PC1+PC2:", 
      round(variance_explained[1] + variance_explained[2], 1), "%\n")
  
  if(group_sep > 20) {
    cat("Strong separation by Diagnosis on PC1\n")
  } else if(group_sep > 10) {
    cat("Moderate separation by Diagnosis\n")
  } else {
    cat("No clear separation by Diagnosis\n")
  }
  
  if(sex_sep > 20) {
    cat("Strong separation by Sex on PC1\n")
  } else if(sex_sep > 10) {
    cat("Moderate separation by Sex\n")
  } else {
    cat("No clear separation by Sex\n")
  }
  
  cor_age <- cor(pca_df$PC1, pca_df$Age)
  cor_rin <- cor(pca_df$PC1, pca_df$RIN)
  cor_tin <- cor(pca_df$PC1, pca_df$TIN)
  
  if(abs(cor_age) > 0.5) {
    cat("Age strongly correlates with PC1 (r =", round(cor_age, 2), ")\n")
    cat("   Consider that Age might be a confounder\n")
  }
  
  if(abs(cor_rin) > 0.5) {
    cat("RIN strongly correlates with PC1 (r =", round(cor_rin, 2), ")\n")
    cat("RNA quality might be driving the variation\n")
  }
  if(abs(cor_tin) > 0.5) {
    cat("TIN strongly correlates with PC1 (r =", round(cor_tin, 2), ")\n")
    cat("TIN quality might be driving the variation\n")
  }
}
interpret_pca(pca_df, variance_explained)

cor(pca_df$RIN, pca_df$TIN)
t.test(RIN ~ Group, data = pca_df)
t.test(TIN ~ Group, data = pca_df)

# STEP 3

resultsNames(dds)
res_interaction <- results(dds, name = "sexf.groupMD")
res_interaction <- res_interaction[order(res_interaction$pvalue), ]
head(res_interaction, 20)
summary(res_interaction)

sig_interaction <- res_interaction[!is.na(res_interaction$padj) & 
                                     res_interaction$padj < 0.05, ]
cat(nrow(sig_interaction))

library(AnnotationDbi)
library(org.Hs.eg.db)

df <- as.data.frame(res_interaction)
df$gene <- rownames(df)
df <- df[df$gene != "ENSG00000099725", ] 

df$gene_clean <- gsub("\\..*$", "", df$gene)
df$symbol <- mapIds(org.Hs.eg.db,
                    keys = df$gene_clean,
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")

df$symbol <- ifelse(is.na(df$symbol), df$gene, df$symbol)

df$significant <- df$pvalue < 0.01 & !is.na(df$pvalue)
sig_df <- df[df$significant == TRUE, ]
top_genes_df <- sig_df[order(sig_df$pvalue), ][1:10, ]

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("gray", "red")) +

  geom_text_repel(data = top_genes_df, aes(label = symbol), 
                  max.overlaps = 20, 
                  size = 3.5, 
                  fontface = "bold",
                  box_padding = 0.4) +
  
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "black", alpha = 0.5) +
  
  labs(x = "log2 Fold Change (Interaction)", 
       y = "-log10(p-value)",
       title = "Genes with Sex-specific response to depression",
       subtitle = "Positive LFC = stronger effect in females (p-value < 0.01)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13))

print(top_genes_df)

ggplot(df, aes(x = log2FoldChange)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "log2 Fold Change (Interaction)",
       y = "Count",
       title = "Distribution of Sex:Diagnosis interaction effects") +
  theme_minimal()

target_gene <- "ENSG00000132744" 
gene_data <- plotCounts(dds, gene = target_gene, intgroup = c("sex", "group"), returnData = TRUE)

ggplot(gene_data, aes(x = sex, y = count, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 2, alpha = 0.6) +
  scale_y_log10() +
  labs(title = paste("Expression profile for ACY3"),
       y = "Normalized Counts (log10)",
       x = "Sex") +
  theme_minimal() +
  scale_fill_manual(values = c("HC" = "grey", "MD" = "orange"))

library(dplyr)
plot_gene_interaction <- function(gene_id) {
  counts_data <- plotCounts(dds, 
                            gene = gene_id, 
                            intgroup = c("group", "sex"), 
                            returnData = TRUE)
  
  ggplot(counts_data, aes(x = group, y = count, color = sex, group = sex)) + 
    geom_point(position = position_jitter(width = 0.1), size = 3, alpha = 0.6) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
    stat_summary(fun = mean, geom = "point", size = 4, shape = 18) +
    scale_y_log10() +
    scale_color_manual(values = c("f" = "#F8766D", "m" = "#00BFC4")) +
    labs(title = paste("Interaction plot for", gene_id),
         x = "Diagnosis Group",
         y = "Normalized Counts (log10)") + 
    theme_minimal() 
}
sig_genes <- df[df$significant, ]

top_genes <- sig_genes$gene[1:min(6, nrow(sig_genes))]
print(top_genes)

plots_list <- list()

for(gene in top_genes) {
  plots_list[[gene]] <- plot_gene_interaction(gene)
}

wrap_plots(plots_list, ncol = 3) + plot_layout(guides = "collect")

# STEP 4

library(RedRibbon)

res_males <- results(dds, name = "group_MD_vs_HC")
res_females <- results(dds, contrast = list(c("group_MD_vs_HC", "sexf.groupMD")))

top_overlap <- function(n) {
  res_m_clean <- res_males[!is.na(res_males$pvalue), ]
  res_f_clean <- res_females[!is.na(res_females$pvalue), ]
  
  male_genes <- rownames(res_m_clean)[order(res_m_clean$pvalue)][1:n]
  female_genes <- rownames(res_f_clean)[order(res_f_clean$pvalue)][1:n]
  
  overlap <- length(intersect(male_genes, female_genes))
  return(overlap)
}

for(n in c(50, 100, 200, 500, 1000)) {
  ov <- top_overlap(n)
  cat("Top", n, "genes overlap:", ov, "(", round(ov/n*100, 1), "%)\n")
}

library(pheatmap)

rld <- rlog(dds, blind = TRUE)
sig_genes <- df[df$significant, ]
top_20_ids <- rownames(sig_genes)[1:20]

mat <- assay(rld)[top_20_ids, ]

rownames(mat) <- gsub("\\..*$", "", rownames(mat))
gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(mat), 
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(mat) <- ifelse(is.na(gene_symbols), rownames(mat), gene_symbols)

anno_col <- data.frame(
  Diagnosis = dds$group,
  Sex = dds$sex
)
rownames(anno_col) <- colnames(mat)

anno_colors <- list(
  Diagnosis = c(HC = "#4DAF4A", MD = "#984EA3"),
  Sex = c(m = "#377EB8", f = "#E41A1C")
)

pheatmap(mat, 
         scale = "row",
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         show_colnames = FALSE,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Hierarchical clustering of top sex-specific genes")

lfc_males <- res_males$log2FoldChange
lfc_females <- res_females$log2FoldChange
valid <- !is.na(lfc_males) & !is.na(lfc_females)

wilcox_p <- wilcox.test(lfc_males[valid], lfc_females[valid], paired = TRUE)$p.value
rep_cor <- cor(lfc_males[valid], lfc_females[valid], method = "spearman")

cat("Wilcoxon test p-value:", format(wilcox_p, scientific = TRUE), "\n")

cor_df <- data.frame(males = lfc_males[valid], females = lfc_females[valid])
cor_df_clean <- cor_df[abs(cor_df$males) < 5 & abs(cor_df$females) < 5, ]

ggplot(cor_df_clean, aes(x = males, y = females)) +
  geom_point(alpha = 0.15, color = "black", size = 0.8) + 
  geom_smooth(method = "lm", color = "blue", linewidth = 1.2) + 
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + 
  labs(x = "log2 Fold Change in Males", 
       y = "log2 Fold Change in Females",
       title = paste0("Global transcriptome correlation (r = ", round(rep_cor_clean, 3), ")"),
       subtitle = "Filtered: extreme outliers (|LFC| >= 5) removed") +
  theme_minimal()

sum(df$significant)

library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)

go_genes <- sig_df$gene

go_enrich <- enrichGO(gene          = go_genes,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENSEMBL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.2)

go_results <- as.data.frame(go_enrich)
head(go_results, 10)

dotplot(go_enrich, showCategory = 15, title = "GO Biological Processes for Sex:Group Interaction") +
  theme_minimal()

library(fgsea)
library(msigdbr)

df_gsea <- df[!is.na(df$log2FoldChange) & !is.na(df$pvalue), ]

df_gsea$gene_clean <- gsub("\\..*", "", df_gsea$gene)

df_gsea$stat <- sign(df_gsea$log2FoldChange) * (-log10(df_gsea$pvalue))

df_gsea <- df_gsea %>% 
  group_by(gene_clean) %>% 
  filter(stat == max(stat)) %>% 
  ungroup()

ranked_genes <- df_gsea$stat
names(ranked_genes) <- df_gsea$gene_clean

ranked_genes <- sort(ranked_genes, decreasing = TRUE)

msig_h <- msigdbr(species = "Homo sapiens", category = "H")

pathways_h <- split(msig_h$ensembl_gene, msig_h$gs_name)

set.seed(42)
fgsea_res <- fgsea(pathways = pathways_h, 
                   stats    = ranked_genes,
                   minSize  = 15, 
                   maxSize  = 500)

fgsea_res_df <- as.data.frame(fgsea_res) %>% arrange(pval)

print(head(fgsea_res_df[, c("pathway", "pval", "padj", "NES", "size")], 10))

top_pathways <- fgsea_res_df %>% head(10)

ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "steelblue"),
                    labels = c("TRUE" = "Enriched in females", "FALSE" = "Enriched in males")) +
  labs(x = "Biological pathway (MSigDB Hallmark)",
       y = "Normalized Enrichment Score (NES)",
       title = "GSEA: Top biological pathways",
       subtitle = "Positive NES = shift towards female response\nNegative NES = shift towards male response",
       fill = "Direction of effect") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "bottom")