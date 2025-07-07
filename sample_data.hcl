// Sample Gene Expression HCL Configuration
version = "1.0"
dataset_name = "sample_gene_expression"

gene_data "BRCA1" {
  expression_values = [2.3, 4.1, 1.8, 3.2, 5.0, 2.7, 4.4, 3.9, 2.1, 3.6]
  chromosome = "17"
  start_position = 43044295
  end_position = 43125483
}

gene_data "TP53" {
  expression_values = [3.8, 2.9, 4.2, 3.1, 2.4, 3.7, 4.0, 2.8, 3.3, 3.9]
  chromosome = "17"
  start_position = 7661779
  end_position = 7687538
}

gene_data "EGFR" {
  expression_values = [1.9, 3.4, 2.8, 4.1, 3.6, 2.2, 3.8, 4.3, 2.9, 3.1]
  chromosome = "7"
  start_position = 55019017
  end_position = 55211628
}

gene_data "MYC" {
  expression_values = [4.2, 3.7, 4.8, 3.9, 4.1, 3.5, 4.4, 4.0, 3.8, 4.3]
  chromosome = "8"
  start_position = 127735434
  end_position = 127742951
}

gene_data "PIK3CA" {
  expression_values = [2.1, 2.8, 3.4, 2.9, 3.1, 2.6, 3.0, 2.7, 3.2, 2.8]
  chromosome = "3"
  start_position = 179148114
  end_position = 179240093
}

experimental_conditions {
  condition_1 = "Control_24h"
  condition_2 = "Treatment_24h"
  condition_3 = "Control_48h"
  condition_4 = "Treatment_48h"
  condition_5 = "Control_72h"
  condition_6 = "Treatment_72h"
  condition_7 = "High_dose_24h"
  condition_8 = "High_dose_48h"
  condition_9 = "High_dose_72h"
  condition_10 = "Recovery_24h"
}

metadata {
  study_type = "time_course"
  organism = "Homo sapiens"
  platform = "RNA-seq"
  normalization = "FPKM"
  date_created = "2024-01-15"
}

matrix_config {
  rows = 5
  cols = 10
  missing_value_threshold = 0.1
  log_transformed = true