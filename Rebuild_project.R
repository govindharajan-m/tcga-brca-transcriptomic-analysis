# Load clinical data
clinical_data <- read.csv("clinical_data.csv")

# Create progression groups
clinical_data$progression_group <- NA

clinical_data$progression_group[
  grepl("^Stage I", clinical_data$ajcc_pathologic_stage)
] <- "Early"

clinical_data$progression_group[
  grepl("^Stage II", clinical_data$ajcc_pathologic_stage)
] <- "Early"

clinical_data$progression_group[
  grepl("^Stage III", clinical_data$ajcc_pathologic_stage)
] <- "Late"

clinical_data$progression_group[
  grepl("^Stage IV", clinical_data$ajcc_pathologic_stage)
] <- "Late"
nrow(clinical_data)

head(clinical_data$ajcc_pathologic_stage)

table(clinical_data$ajcc_pathologic_stage)

length(clinical_data$progression_group)

clinical_data$progression_group <- ifelse(
  grepl("^Stage I|^Stage II", clinical_data$ajcc_pathologic_stage),
  "Early",
  ifelse(
    grepl("^Stage III|^Stage IV", clinical_data$ajcc_pathologic_stage),
    "Late",
    NA
  )
)

table(clinical_data$progression_group, useNA = "ifany")


dim(read.csv("clinical_data.csv"))

dim(read.csv("clinical_data_cleaned.csv"))


library(TCGAbiolinks)

clinical_data <- GDCquery_clinic(
  project = "TCGA-BRCA",
  type = "clinical"
)

dim(clinical_data)

saveRDS(
  clinical_data,
  "clinical_data.rds"
)


clinical_data$progression_group <- ifelse(
  grepl("^Stage I|^Stage II",
        clinical_data$ajcc_pathologic_stage),
  "Early",
  ifelse(
    grepl("^Stage III|^Stage IV",
          clinical_data$ajcc_pathologic_stage),
    "Late",
    NA
  )
)

table(clinical_data$progression_group, useNA = "ifany")


table(clinical_data$ajcc_pathologic_stage)

clinical_data$progression_group <- NA

clinical_data$progression_group[
  grepl("^Stage IV", clinical_data$ajcc_pathologic_stage)
] <- "Late"

clinical_data$progression_group[
  grepl("^Stage III", clinical_data$ajcc_pathologic_stage)
] <- "Late"

clinical_data$progression_group[
  grepl("^Stage II", clinical_data$ajcc_pathologic_stage)
] <- "Early"

clinical_data$progression_group[
  grepl("^Stage I[A-B]?$|^Stage I$", clinical_data$ajcc_pathologic_stage)
] <- "Early"

table(clinical_data$progression_group, useNA = "ifany")

table(clinical_data$ajcc_pathologic_stage)

clinical_data$analysis_group <- NA

clinical_data$analysis_group[
  clinical_data$ajcc_pathologic_stage %in%
    c("Stage I","Stage IA","Stage IB")
] <- "Early"

clinical_data$analysis_group[
  clinical_data$ajcc_pathologic_stage %in%
    c(
      "Stage III",
      "Stage IIIA",
      "Stage IIIB",
      "Stage IIIC",
      "Stage IV"
    )
] <- "Advanced"

table(clinical_data$analysis_group, useNA = "ifany")

head(clinical_data$submitter_id)

table(clinical_data$analysis_group, useNA = "ifany")
early_ids <- clinical_data$submitter_id[
  clinical_data$analysis_group == "Early"
]

advanced_ids <- clinical_data$submitter_id[
  clinical_data$analysis_group == "Advanced"
]

length(early_ids)
length(advanced_ids)

head(early_ids)
head(advanced_ids)

early_ids <- clinical_data$submitter_id[
  !is.na(clinical_data$analysis_group) &
    clinical_data$analysis_group == "Early"
]

advanced_ids <- clinical_data$submitter_id[
  !is.na(clinical_data$analysis_group) &
    clinical_data$analysis_group == "Advanced"
]

length(early_ids)
length(advanced_ids)

sum(is.na(early_ids))
sum(is.na(advanced_ids))

head(early_ids)
head(advanced_ids)


grep("TCGA", list.files(recursive = TRUE), value = TRUE)

length(early_ids)
length(advanced_ids)

sum(is.na(early_ids))
sum(is.na(advanced_ids))

results <- getResults(query.exp)

colnames(results)
  
head(results)

results_subset <- results[
  results$patient_id %in%
    c(early_ids, advanced_ids),
]

saveRDS(clinical_data, "clinical_data.rds")

saveRDS(early_ids, "early_ids.rds")

saveRDS(advanced_ids, "advanced_ids.rds")

colnames(results)

head(results)

results <- getResults(query.exp)

rna_meta <- data.frame(
  patient_id = results$cases.submitter_id,
  sample_id = results$sample.submitter_id,
  file_name = results$file_name,
  stringsAsFactors = FALSE
)

head(rna_meta)

selected_ids <- c(early_ids, advanced_ids)

rna_meta_subset <- rna_meta[
  rna_meta$patient_id %in% selected_ids,
]

dim(rna_meta_subset)
table(rna_meta_subset$patient_id %in% early_ids)
table(rna_meta_subset$patient_id %in% advanced_ids)

length(unique(rna_meta_subset$patient_id))

rna_meta_subset <- rna_meta_subset[
  grepl("-01A$", rna_meta_subset$sample_id),
]

dim(rna_meta_subset)


dim(rna_meta_subset)

length(unique(rna_meta_subset$patient_id))

head(rna_meta_subset)


saveRDS(clinical_data, "clinical_data.rds")
saveRDS(rna_meta_subset, "rna_meta_subset.rds")
saveRDS(early_ids, "early_ids.rds")
saveRDS(advanced_ids, "advanced_ids.rds")

all_files <- list.files(
  recursive = TRUE,
  pattern = "star_gene_counts\\.tsv$",
  full.names = TRUE
)

matched_files <- all_files[
  basename(all_files) %in% rna_meta_subset$file_name
]

test <- read.delim(
  matched_files[1],
  header = TRUE
)

colnames(test)

head(test)

readLines(matched_files[1], n = 15)

test <- read.delim(
  matched_files[1],
  comment.char = "#"
)

head(test)

small_files <- matched_files[1:5]

small_counts <- lapply(small_files, function(f){
  
  x <- read.delim(
    f,
    comment.char = "#"
  )
  
  x[, c("gene_id", "unstranded")]
  
})

head(small_counts[[1]])

nrow(small_counts[[1]])

small_files <- matched_files[1:5]

small_counts <- lapply(small_files, function(f){
  
  x <- read.delim(
    f,
    comment.char = "#"
  )
  
  x[, c("gene_id", "unstranded")]
  
})

head(small_counts[[1]])
nrow(small_counts[[1]])


small_files <- matched_files[1:5]

small_counts <- lapply(small_files, function(f) {
  
  x <- read.delim(
    f,
    comment.char = "#"
  )
  
  x[, c("gene_id", "unstranded")]
  
})

head(small_counts[[1]])

nrow(small_counts[[1]])

test <- read.delim(
  matched_files[1],
  comment.char = "#"
)

head(test)

colnames(test)

test_files <- matched_files[1:20]

count_list <- lapply(test_files, function(f){
  
  x <- read.delim(
    f,
    comment.char = "#"
  )
  
  x <- x[grep("^ENSG", x$gene_id), ]
  
  x$unstranded
  
})

length(count_list)

length(count_list[[1]])

count_matrix <- do.call(
  cbind,
  count_list
)

dim(count_matrix)

gene_table <- read.delim(
  test_files[1],
  comment.char = "#"
)

gene_table <- gene_table[
  grep("^ENSG", gene_table$gene_id),
]

rownames(count_matrix) <- gene_table$gene_id

head(count_matrix)

dim(count_matrix)

head(count_matrix[,1:3])

save.image("breast_cancer_workspace.RData")

gene_table <- read.delim(
  matched_files[1],
  comment.char = "#"
)

gene_table <- gene_table[
  grep("^ENSG", gene_table$gene_id),
]

genes <- gene_table$gene_id

length(genes)

count_list <- lapply(matched_files, function(f){
  
  x <- read.delim(
    f,
    comment.char = "#"
  )
  
  x <- x[
    grep("^ENSG", x$gene_id),
  ]
  
  x$unstranded
  
})

count_matrix <- do.call(
  cbind,
  count_list
)

rownames(count_matrix) <- genes

colnames(count_matrix) <- rna_meta_subset$sample_id

dim(count_matrix)

saveRDS(
  count_matrix,
  "count_matrix.rds"
)

dim(count_matrix)


saveRDS(count_matrix, "count_matrix.rds")

sample_info <- data.frame(
  sample_id = rna_meta_subset$sample_id,
  patient_id = rna_meta_subset$patient_id,
  stringsAsFactors = FALSE
)

list.files(pattern = "count_matrix")

count_matrix <- readRDS("count_matrix.rds")
dim(count_matrix)

list.files()

clinical_data <- readRDS("clinical_data.rds")
count_matrix <- readRDS("count_matrix.rds")
rna_meta_subset <- readRDS("rna_meta_subset.rds")

clinical_data <- readRDS("clinical_data.rds")

count_matrix <- readRDS("count_matrix.rds")

rna_meta_subset <- readRDS("rna_meta_subset.rds")

early_ids <- readRDS("early_ids.rds")

advanced_ids <- readRDS("advanced_ids.rds")

sample_info <- data.frame(
  sample_id = rna_meta_subset$sample_id,
  patient_id = rna_meta_subset$patient_id,
  stringsAsFactors = FALSE
)


group_lookup <- clinical_data[
  !is.na(clinical_data$analysis_group),
  c("submitter_id", "analysis_group")
]

sample_info <- merge(
  sample_info,
  group_lookup,
  by.x = "patient_id",
  by.y = "submitter_id"
)

dim(sample_info)

table(sample_info$analysis_group)

all(sample_info$sample_id %in% colnames(count_matrix))

colnames(sample_info)

table(sample_info$analysis_group)

all(sample_info$sample_id %in% colnames(count_matrix))

colnames(sample_info)
head(sample_info)

group_lookup <- clinical_data[
  !is.na(clinical_data$analysis_group),
  c("submitter_id", "analysis_group")
]

sample_info <- data.frame(
  sample_id = rna_meta_subset$sample_id,
  patient_id = rna_meta_subset$patient_id,
  stringsAsFactors = FALSE
)

sample_info <- merge(
  sample_info,
  group_lookup,
  by.x = "patient_id",
  by.y = "submitter_id"
)
table(sample_info$analysis_group)

nrow(sample_info)

table(sample_info$analysis_group)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

dim(count_matrix)

dim(sample_info)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = sample_info,
  design = ~ analysis_group
)
library(DESeq2)

dds <- dds[
  rowSums(counts(dds)) >= 10,
]

dim(dds)

dds <- DESeq(dds)

saveRDS(dds, "dds.rds")

res_df <- as.data.frame(res)

write.csv(
  res_df,
  "deseq2_results.csv"
)

exists("dds")
exists("res")

dds <- DESeq(dds)

res <- results(dds)

head(res)

summary(res)

dds <- DESeq(dds)

dim(dds)

dds <- dds[
  rowSums(counts(dds)) >= 10,
]

dim(dds)

dds <- DESeq(dds)

res <- results(dds)

summary(res)

res_df <- as.data.frame(res)

write.csv(
  res_df,
  "deseq2_results.csv",
  row.names = TRUE
)

saveRDS(
  res,
  "deseq2_results.rds"
)

summary(res)

res_df <- as.data.frame(res)

sig_genes <- subset(
  res_df,
  padj < 0.05
)

nrow(sig_genes)

res <- results(dds)

summary(res)

head(res)

summary(res)

res_df <- as.data.frame(res)

sig_genes <- subset(
  res_df,
  !is.na(padj) & padj < 0.05
)

nrow(sig_genes)

sig_genes <- sig_genes[
  order(sig_genes$padj),
]
head(sig_genes, 20)

all_files <- list.files(
  recursive = TRUE,
  pattern = "star_gene_counts\\.tsv$",
  full.names = TRUE
)

length(all_files)

matched_files <- all_files[
  basename(all_files) %in% rna_meta_subset$file_name
]

length(matched_files)
head(matched_files)

gene_table <- read.delim(
  matched_files[1],
  comment.char = "#"
)

gene_table <- gene_table[
  grep("^ENSG", gene_table$gene_id),
]

head(gene_table[, c("gene_id", "gene_name")])
exists("res")
exists("res_df")

summary(res)
res_df <- as.data.frame(res)

res_df$gene_id <- rownames(res_df)

res_named <- merge(
  res_df,
  gene_table[, c("gene_id", "gene_name")],
  by = "gene_id"
)

sig_genes <- subset(
  res_named,
  !is.na(padj) & padj < 0.05
)

nrow(sig_genes)

sig_genes <- sig_genes[
  order(sig_genes$padj),
]

head(
  sig_genes[, c(
    "gene_name",
    "log2FoldChange",
    "padj"
  )],
  30
)

write.csv(
  sig_genes,
  "significant_genes.csv",
  row.names = FALSE
)

saveRDS(
  sig_genes,
  "significant_genes.rds"
)

library(ggplot2)

res_named$significant <- ifelse(
  !is.na(res_named$padj) &
    res_named$padj < 0.05,
  "Significant",
  "Not Significant"
)

ggplot(
  res_named,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significant
  )
) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Early vs Advanced Breast Cancer",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

ggsave(
  "volcano_plot.png",
  width = 8,
  height = 6
)

  library(ggplot2)
  
  volcano_data <- subset(
    res_named,
    !is.na(padj)
  )
  
  p <- ggplot(
    volcano_data,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = significant
    )
  ) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(
      title = "Differential Expression: Early vs Advanced Breast Cancer",
      x = "Log2 Fold Change",
      y = "-log10 Adjusted P-value"
    )
  
  p
  
  p <- ggplot(
    volcano_data,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = significant
    )
  ) +
    geom_point(alpha = 0.5) +
    theme_minimal()
  
  print(p)
  
  graphics.off()
  print(p)
  
  file.exists("volcano_plot.png")
  file.info("volcano_plot.png")$size
  
  sig_genes <- sig_genes[
    order(sig_genes$padj),
  ]
  
  head(
    sig_genes[, c(
      "gene_name",
      "log2FoldChange",
      "padj"
    )],
    30
  )
  
  write.csv(
    res_named,
    "all_deseq_results.csv",
    row.names = FALSE
  )
  
  write.csv(
    sig_genes,
    "significant_genes.csv",
    row.names = FALSE
  )
  graphics.off()
  windows()
  print(p)
  
  head(
    sig_genes[, c(
      "gene_name",
      "log2FoldChange",
      "padj"
    )],
    30
  )
  
  up_genes <- subset(
    sig_genes,
    log2FoldChange > 1
  )
  
  down_genes <- subset(
    sig_genes,
    log2FoldChange < -1
  )
  
  nrow(up_genes)
  nrow(down_genes)
  
  top50 <- sig_genes[
    order(sig_genes$padj),
  ][1:50, ]
  
  write.csv(
    top50,
    "top50_genes.csv",
    row.names = FALSE
  )
  
  top50[, c("gene_name","log2FoldChange","padj")]
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  BiocManager::install("enrichplot")
  
  BiocManager::install("clusterProfiler")
  BiocManager::install("org.Hs.eg.db")
  write.csv(
    sig_genes,
    "significant_genes.csv",
    row.names = FALSE
  )
  
  write.csv(
    res_named,
    "all_deseq_results.csv",
    row.names = FALSE
  )
  
  saveRDS(
    dds,
    "dds.rds"
  )
  
  saveRDS(
    res,
    "res.rds"
  )
  
  up_genes <- subset(
    sig_genes,
    log2FoldChange > 1
  )
  
  down_genes <- subset(
    sig_genes,
    log2FoldChange < -1
  )
  
  nrow(up_genes)
  nrow(down_genes)
  
  top50 <- sig_genes[
    order(sig_genes$padj),
  ][1:50, ]
  
  write.csv(
    top50,
    "top50_genes.csv",
    row.names = FALSE
  )
  
  top50[, c("gene_name","log2FoldChange","padj")]
  
  gene_symbols <- unique(sig_genes$gene_name)
  
  gene_ids <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  head(gene_ids)
  
  BiocManager::install("enrichplot")
  
  BiocManager::install("clusterProfiler")
  
  library(enrichplot)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  search()
  
  head(
    sig_genes[order(abs(sig_genes$log2FoldChange), decreasing = TRUE),
              c("gene_name","log2FoldChange","padj")],
    30
  )
  
  table( 
    cut(
      abs(sig_genes$log2FoldChange),
      breaks = c(0,1,2,3,Inf)
    )
  )
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  gene_symbols <- unique(sig_genes$gene_name)
  
  gene_ids <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  head(gene_ids)
  
  BiocManager::install("org.Hs.eg.db")
  exists("org.Hs.eg.db")
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  table(clinical_data$ajcc_pathologic_stage)
  
  install.packages("pheatmap")
  
  top20 <- head(sig_genes[order(sig_genes$padj), ], 20)
  
  top_ids <- top20$gene_id
  
  heatmap_counts <- count_matrix[
    rownames(count_matrix) %in% top_ids,
  ]
  
  pheatmap::pheatmap(
    log2(heatmap_counts + 1),
    show_rownames = TRUE
  )
  
  group_counts <- as.data.frame(
    table(sample_info$analysis_group)
  )
  
  colnames(group_counts) <- c("Group", "Count")
  
  ggplot(group_counts,
         aes(
           x = Group,
           y = Count,
           fill = Group
         )) +
    geom_bar(
      stat = "identity",
      width = 0.7
    ) +
    geom_text(
      aes(label = Count),
      vjust = -0.4,
      size = 5
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Sample Distribution for Differential Expression Analysis",
      x = "",
      y = "Number of Samples"
    ) +
    theme(
      legend.position = "none"
    )
  
  ggsave(
    "early_vs_advanced_samples.png",
    width = 6,
    height = 5,
    dpi = 300
  )
  exists("sample_info")
  exists("res")
  exists("res_named")
  exists("sig_genes")
  
  "analysis_group" %in% colnames(sample_info)
  
  table(sample_info$analysis_group, useNA = "ifany")
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  gene_table <- read.delim(
    matched_files[1],
    comment.char = "#"
  )
  
  gene_table <- gene_table[
    grep("^ENSG", gene_table$gene_id),
  ]
  
  res_named <- merge(
    res_df,
    gene_table[, c("gene_id", "gene_name")],
    by = "gene_id"
  )
  
  sig_genes <- subset(
    res_named,
    !is.na(padj) & padj < 0.05
  )
  
  nrow(sig_genes)
  matched_files <- list.files(
    "GDCdata",
    pattern = "star_gene_counts\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  length(matched_files)
  gene_table <- read.delim(
    matched_files[1],
    comment.char = "#"
  )
  
  gene_table <- gene_table[
    grep("^ENSG", gene_table$gene_id),
  ]
  
  head(gene_table[, c("gene_id", "gene_name")])
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  res_named <- merge(
    res_df,
    gene_table[, c("gene_id", "gene_name")],
    by = "gene_id"
  )
  
  sig_genes <- subset(
    res_named,
    !is.na(padj) & padj < 0.05
  )
  
  nrow(sig_genes)
  
  top20 <- head(
    sig_genes[order(sig_genes$padj), ],
    20
  )
  
  top20[, c("gene_name", "log2FoldChange", "padj")]
  
  write.csv(sig_genes, "significant_genes.csv", row.names = FALSE)
  write.csv(top20, "top20_genes.csv", row.names = FALSE)
  saveRDS(sig_genes, "sig_genes.rds")
  
  top20 <- head(
    sig_genes[order(sig_genes$padj), ],
    20
  )
  
  top20[, c(
    "gene_name",
    "log2FoldChange",
    "padj"
  )]
  
  write.csv(
    top20,
    "top20_genes.csv",
    row.names = FALSE
  )
  
  install.packages("pheatmap")
  library(pheatmap)
  
  top20_ids <- top20$gene_id
  
  heatmap_counts <- count_matrix[
    rownames(count_matrix) %in% top20_ids,
  ]
  
  pheatmap(
    log2(heatmap_counts + 1),
    scale = "row",
    show_rownames = TRUE,
    show_colnames = FALSE
  )
  png(
    "top20_heatmap.png",
    width = 1400,
    height = 1200
  )
  
  pheatmap(
    log2(heatmap_counts + 1),
    scale = "row",
    show_rownames = TRUE,
    show_colnames = FALSE
  )
  
  dev.off()
  
  write.csv(
    sig_genes,
    "significant_genes.csv",
    row.names = FALSE
  )
  sum(sig_genes$log2FoldChange > 1)
  sum(sig_genes$log2FoldChange < -1)
  summary(sig_genes$log2FoldChange)
  
  head(
    top20[, c(
      "gene_name",
      "log2FoldChange",
      "padj"
    )],
    20
  )
  
  sum(sig_genes$log2FoldChange > 1)
  sum(sig_genes$log2FoldChange < -1)
  dev.list()
  while (!is.null(dev.list())) dev.off()
  
  sum(sig_genes$log2FoldChange > 1)
  sum(sig_genes$log2FoldChange < -1)
  
  summary(sig_genes$log2FoldChange)
  write.csv(
    sig_genes,
    "significant_genes.csv",
    row.names = FALSE
  )
  top20[, c(
    "gene_name",
    "log2FoldChange",
    "padj"
  )]
  library(DESeq2)
  resultsNames(dds)  results(dds, contrast = c("group", "Early", "Advanced"))
  results(dds, contrast = c("group", "Advanced", "Early"))
  
  library(DESeq2)
  
  dds <- readRDS("dds.rds")
  
  resultsNames(dds)
  levels(colData(dds)$analysis_group)
  
  dds <- readRDS("dds.rds")
  
  ls()
  file.exists("dds.rds")
  file.exists("res.rds")
  
  library(DESeq2)
  
  resultsNames(dds)
  
  levels(colData(dds)$analysis_group)
  head(sig_genes[
    order(sig_genes$padj),
    c("gene_name","log2FoldChange","padj")
  ], 10)
  
  resultsNames(dds)
  
  head(
    sig_genes[sig_genes$gene_name %in%
                c("CHGA","KRT6A","KRT13","KRT20","CR2"),
              c("gene_name","log2FoldChange","padj")]
  )