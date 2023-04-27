# S. Spielman for ALSF CCDL & Jo Lynne Rokita for D3b, 2022
#
# Makes a pdf panel of forest plot of survival analysis on MB samples 
#  with immune cell fractions and PDL-1 expression predictors


library(survival) # needed to parse model output
library(tidyverse)
library(patchwork) # for this forest plot export, patchwork is better than cowplot for lining up p-values

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig5", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Zenodo CSV output directory and file path
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
fig5d_csv <- file.path(zenodo_tables_dir, "figure-5d-data.csv")


# Directory with result data to input for plot
input_dir <- file.path(root_dir, "analyses", "survival-analysis", "results", "immune")


## Input and output files -------------------------------------
# note that $model contains the model, and $data is the data that went into the model
survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_per_Medulloblastoma_terms_quantiseq.RDS"
  )
)

forest_pdf <- file.path(output_dir, "forest_MB_quantiseq_CD274.pdf")

## Make plot --------------------------------------------------

# NOTE:
# The `Macrophage_M1` lower and upper bounds are INFINITE, with an HR = 0 and p = 0.99.
# This coefficient will therefore NOT be included in the plot, so it must be discussed in text.

ref_term <- "Tumor resection: Biopsy (ref)"

# Set up ordering and labels for y-axis
term_order <- rev(c("CD274",
                "B_cell", 
                "Macrophage_M2",  
                "Monocyte", 
                "Myeloid_dendritic_cell",
                "Neutrophil",          
                "NK_cell",
                "T_cell_CD4_non_regulatory",
                "T_cell_CD8",               
                "T_cell_regulatory_Tregs",
                "extent_of_tumor_resectionGross/Near total resection",
                "extent_of_tumor_resectionPartial resection",
                ref_term))

term_labels <- rev(c("CD274 expression (FPKM)",
                "B cell", 
                "Macrophage M2",  
                "Monocyte", 
                "Myeloid dendritic cell",
                "Neutrophil",          
                "NK cell",
                "CD4+ T cell (non-regulatory)",
                "CD8+ T cell",               
                "Regulatory T cell (Treg)",
                "Tumor resection: Total",
                "Tumor resection: Partial",
                ref_term))


# Get n and event info from glance outpout
survival_n <- broom::glance(survival_result$model) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result$model) %>%
  # add reference
  add_row(term = ref_term, estimate = 0) %>%
  # remove the unbounded Macrophage_M1:
  filter(term != "Macrophage_M1") %>%
  mutate(estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high), 
         # significance indicator column for filling points.
         # Note T/F these are strings for type compatibility with "REF"
         significant = case_when(p.value <= 0.05 ~ "TRUE", 
                                 p.value > 0.05 ~ "FALSE", 
                                 is.na(p.value) ~ "REF"),
         # y-axis factor re-labeling
         term = factor(term, 
                       levels = term_order,
                       labels = term_labels)
         )



# Forest plot of the model
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.25
  ) + 
  geom_point(
    size = 4.5,
    stroke = 0.25,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "REF" = "gray", 
               "TRUE" = "black"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3, 
    size = 0.25
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold", size = rel(1.1)),
    # thinner axes, ticks for compilation
    axis.line = element_line(size = rel(0.25)),
    axis.ticks = element_line(size = rel(0.25)),
    plot.margin = margin(0,0,0,0)
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()

# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  # remove CI for ref
  mutate(value = ifelse(grepl("ref", term), NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0, size = 4) +
  labs(
    # hack!
    subtitle = paste0("                            ",
                      "HR (95% CI)                             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get covered
    plot.margin = margin(6, 0, 36, -26, unit = "pt"),
    plot.subtitle = element_text(face = "bold", size = rel(1.1))
  ) 


forest_panels <- forest_plot + labels_panel + plot_layout(widths = c(1, 1.2))

# Export plot
ggsave(forest_pdf, forest_panels, width = 12, height = 4)


# Export CSV for Zenodo upload
survival_result$data %>%
  # order columns and only keep columns that are in the model:
  # OS_status and OS_years
  # ~ <quantiseq terms>+CD274+extent_of_tumor_resection"
  dplyr::select(Kids_First_Biospecimen_ID, 
                OS_years, 
                OS_status,
                extent_of_tumor_resection,
                CD274,
                B_cell:last_col()) %>%
  # recode OS_status back to LIVING/DECEASED
  dplyr::mutate(OS_status = ifelse(
    OS_status == 1, 
    "LIVING", 
    "DECEASED")
  ) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_csv(fig5d_csv)
