clinical_data$progression_group <- ifelse(
  clinical_data$simple_stage %in% c("Stage I", "Stage II"),
  "Early",
  "Late"
)

table(clinical_data$progression_group)
table(clinical_data$simple_stage)
unique(clinical_data$simple_stage)
"simple_stage" %in% names(clinical_data)
table(clinical_data$ajcc_pathologic_stage)
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

table(clinical_data$progression_group)
table(clinical_data$progression_group)
sapply(clinical_data, class)
library(TCGAbiolinks)

query.exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
getResults(query.exp)
nrow(getResults(query.exp))
unique(getResults(query.exp)$data_type)
nrow(getResults(query.exp))
unique(getResults(query.exp)$workflow_type)
colnames(getResults(query.exp))
unique(getResults(query.exp)$workflow_type)
GDCdownload(query.exp)
getwd()

length(list.files(
  "GDCdata/TCGA-BRCA",
  recursive = TRUE
))

list.files()

length(list.files(
  recursive = TRUE,
  pattern = "\\.tsv$"
))

file.info("Tue_Jun__9_16_28_18_2026_0.tar.gz")

length(list.files(recursive = TRUE))
length(list.dirs(recursive = FALSE))
length(list.files(recursive = TRUE))
library(TCGAbiolinks)
exp_data <- GDCprepare(query.exp)
GDCprepare
library(TCGAbiolinks)

query.exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

nrow(getResults(query.exp))

list.files(
  recursive = TRUE,
  pattern = "\\.tsv$",
  full.names = TRUE
)[1:5]

length(
  list.files(
    recursive = TRUE,
    pattern = "star_gene_counts\\.tsv$"
  )
)

library(TCGAbiolinks)

query.exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

exp_data <- GDCprepare(query.exp)
colnames(clinical_data)
head(clinical_data$submitter_id)


length(
  list.files(
    recursive = TRUE,
    pattern = "star_gene_counts\\.tsv$"
  )
)

gc()
table(clinical_data$ajcc_pathologic_stage)
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

table(clinical_data$progression_group)