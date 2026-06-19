library(TCGAbiolinks)
library(dplyr)
library(ggplot2)
ggplot(clinical_data,
       aes(x = vital_status)) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "Patient Survival Status",
    x = "Status",
    y = "Count"
  )
