library(tidyverse)

scr_dir <- dirname("/fslhome/nmella/OpenPBTA-analysis")

setwd(paste0(scr_dir, "/data/release-v13-20200116"))
print(scr_dir)

# dat_rsem_poly1 = readRDS("pbta-gene-expression-rsem-fpkm-collapsed.polya.rds")
# dat_rsem_poly2 = readRDS("pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")
# dat_rsem_poly3 = readRDS("pbta-gene-expression-rsem-fpkm.polya.rds")
# dat_rsem_poly4 = readRDS("pbta-gene-expression-rsem-fpkm.stranded.rds")
# dat_rsem_poly5 = readRDS("pbta-gene-expression-rsem-tpm.polya.rds")
# dat_rsem_poly6 = readRDS("pbta-gene-expression-rsem-tpm.stranded.rds")
dat_kallisto_stranded <- readRDS("pbta-gene-expression-kallisto.stranded.rds")
dat_kallisto_stranded = dat_kallisto_stranded[,2:ncol(dat_kallisto_stranded)]
dat_kallisto_poly <- readRDS("pbta-gene-expression-kallisto.polya.rds")

setwd(paste0(scr_dir, "/data/release-v13-20200116"))
covariate = read_tsv("pbta-histologies.tsv",col_types = cols(molecular_subtype = "c"))

covariate = filter(covariate, experimental_strategy == "RNA-Seq")
covariate = unite(covariate, col = "batch", c("cohort", "seq_center"), sep = "_", remove = TRUE, na.rm = FALSE)
covariate$batch = as_factor(covariate$batch)
covariate = covariate %>% select(1:4, batch, everything())
# covariate = gather(covariate, key = "condition", value = "condition_value", 7:ncol(covariate))
# covariate$condition_value = as_factor(covariate$condition_value)




remove_char = function(a){
  b = NULL
  for (x in a){
    if (!is.na(as.numeric(x))){
      b = c(b, x)
    }
  }
  return(b)
}

categorize = function(a){
  
  b = remove_char(a)
  
  bound1 = quantile(b)[[2]]
  bound2 = quantile(b)[[3]]
  bound3 = quantile(b)[[4]]
  bound4 = quantile(b)[[5]]
  
  i = 1
  for (x in a){
    if (!is.na(as.numeric(x))){
      x = as.integer(x)
      if (x < bound2){
        a[i] = 1
      }else if (x >= bound2 && x < bound3){
        a[i] = 2
      }else if (x >= bound3 & x < bound4){
        a[i] = 3
      }else {
        a[i] = 4
      }
    }
    i = i + 1
  }
  return(a)
}

covariate$age_at_diagnosis_days = categorize(as.numeric(covariate$age_at_diagnosis_days))
covariate = unite(covariate, col = "condition", c("age_at_diagnosis_days", "reported_gender"), sep = "_", remove = TRUE, na.rm = TRUE)
covariate$condition = as_factor(covariate$condition)
covariate = covariate %>% select(1:4, condition, everything())


trans = function(data, Kids_First_Biospecimen_ID, val, starting_col){
  data = data %>%
    gather(Kids_First_Biospecimen_ID, val, starting_col:(ncol(data)-1)) %>%
    spread("gene_name", val)
  return(data)
}

clean = function(data, gene_name){
  data_f = mutate(data, gene_name) 
  data_t = trans(data_f, "id", "gene_expression", 1)
  print(data_t)
  return(data_t)
}

id_batch = cbind("Kids_First_Biospecimen_ID" = covariate$Kids_First_Biospecimen_ID, "batch" = covariate$batch, "condition" = covariate$condition,"condition_value" = covariate$condtion_value)
colnames(id_batch) <- c("Kids_First_Biospecimen_ID", "batch", "condition")
id_batch = as_tibble(id_batch)


dat_kallisto_poly = mutate(dat_kallisto_poly, "batch" = rep("polya", times = nrow(dat_kallisto_poly)))
dat_kallisto_stranded = mutate(dat_kallisto_stranded, "batch" = rep("stranded", times = nrow(dat_kallisto_stranded)))

data = full_join(dat_kallisto_poly, dat_kallisto_stranded)



# dat_rsem_poly1_t = clean(dat_rsem_poly1, "gene_name" = rownames(dat_rsem_poly1))
# dat_rsem_poly1_t = inner_join(id_batch, dat_rsem_poly1_t, by = "Kids_First_Biospecimen_ID")


# poly1batch = as_factor(dat_rsem_poly1_t$batch)
# levels(poly1batch) <- 1:length(levels(poly1batch))
# poly1condition = as_factor(dat_rsem_poly1_t$condition)
# levels(poly1condition) <- 1:length(levels(poly1condition))
# poly1matrix = t(dat_rsem_poly1_t[4:ncol(dat_rsem_poly1_t)])
# i = 1
# poly1condition = as.integer(poly1condition)
# for (x in poly1condition){
#   if (x > mean(poly1condition)){
#      poly1condition[i] = 2
#   }
#   else {
#     poly1condition[i] = 1
#   }
#   i = i + 1
# }



# dat_rsem_poly2_t = clean(dat_rsem_poly2)
# dat_rsem_poly3_t = clean(dat_rsem_poly3)
# dat_rsem_poly4_t = clean(dat_rsem_poly4)
# dat_rsem_poly5_t = clean(dat_rsem_poly5)
# dat_rsem_poly6_t = clean(dat_rsem_poly6)

library(dplyr)
library(tidyr)


df = dat_kallisto_stranded[,2:ncol(dat_kallisto_stranded)]

df_t = df %>%
  gather(Kids_First_Biospecimen_ID, gene_expression, 2:ncol(df)) %>%
  spread(gene_id, gene_expression)

df = dat_kallisto_poly[,2:ncol(dat_kallisto_poly)]

df_t2 = df %>%
  gather(Kids_First_Biospecimen_ID, gene_expression, 2:ncol(df)) %>%
  spread(gene_id, gene_expression)

df_t = df_t %>%
  mutate("batch" = rep(1, times = nrow(df_t)))

df_t2 = df_t2 %>%
  mutate("batch" = rep(2, times = nrow(df_t2)))

data = full_join(df_t, df_t2)

data = data %>% select(1:1, batch, everything())

id_batch2 = rename(id_batch, "condition_1" = batch) %>% rename("condition_2" = condition)

fdata = inner_join(id_batch2, data)

fdata = arrange(fdata, condition_1)




poly1batch = as_factor(fdata$batch)
levels(poly1batch) <- 1:length(levels(poly1batch))
poly1condition = as_factor(fdata$condition_1)
levels(poly1condition) <- 1:length(levels(poly1condition))
poly1matrix = t(fdata[5:ncol(fdata)])

poly1batch = as.numeric(poly1batch)
poly1condition = as.numeric(poly1condition)




library(BatchQC)


batchQC(poly1matrix, batch=poly1batch, condition=poly1condition,
        report_file="batchqc_report.html", report_dir=".",
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE, batchqc_output=TRUE)


