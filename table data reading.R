ls()
getwd()
list.files()
list.files(recursive = TRUE)
exists("clinical_data")
list.files(pattern = ".csv")
read.csv("clinical_data.csv")
names(clinical_data)
head(clinical_data)
dim(clinical_data)
names(clinical_data)
table(clinical_data$ajcc_pathologic_stage)
clinical_data$simple_stage <- gsub(
  "Stage ([IV]+).*",
  "Stage \\1",
  clinical_data$ajcc_pathologic_stage
)
grep("death|follow", names(clinical_data), value = TRUE)
clinical_data$survival_time <- ifelse(
  is.na(clinical_data$days_to_death),
  clinical_data$days_to_last_follow_up,
  clinical_data$days_to_death
)
table(clinical_data$vital_status)
grep("age", names(clinical_data), value = TRUE)
clinical_data$age_group <- ifelse(
  clinical_data$age_at_index < 40,
  "Young",
  "Older"
)
names(clinical_data)
grep("er|pr|her2|race|ethnicity|subtype", names(clinical_data),
     ignore.case = TRUE,
     value = TRUE)
table(clinical_data$progression_or_recurrence)
colSums(!is.na(clinical_data))
sort(colSums(!is.na(clinical_data)), decreasing = TRUE)
grep("stage", names(clinical_data), ignore.case = TRUE, value = TRUE)

grep("death", names(clinical_data), ignore.case = TRUE, value = TRUE)

grep("follow", names(clinical_data), ignore.case = TRUE, value = TRUE)

grep("age", names(clinical_data), ignore.case = TRUE, value = TRUE)
clinical_data$survival_time <- ifelse(
  is.na(clinical_data$days_to_death),
  clinical_data$days_to_last_follow_up,
  clinical_data$days_to_death
)

clinical_data$event <- ifelse(
  clinical_data$vital_status == "Dead",
  1,
  0
)

summary(clinical_data$survival_time)
table(clinical_data$event)
query.exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)