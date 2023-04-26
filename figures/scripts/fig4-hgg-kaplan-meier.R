# S. Spielman for ALSF CCDL 2022-3
#
# Makes a pdf panel for the HGG Kaplan-Meier survival analysis for Figure 4

library(tidyverse)
library(survival)
library(patchwork)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig4", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}



# Zenodo CSV output directory and file path
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
fig4h_csv <- file.path(zenodo_tables_dir, "figure-4h-data.csv")

# Input directory
input_dir <- file.path(root_dir, "analyses", "survival-analysis", "results", "subtypes")



## Input and output files ------------------------------------------------
km_result <- read_rds(
  file.path(input_dir,
            "logrank_hgg_subtypes.RDS"
  )
)
km_output_pdf <- file.path(output_dir, "km_hgg_panel.pdf")


# Establish color-blind friendly colors
okabe_palette <- colorblindr::palette_OkabeIto[1:5]


# We need to separately calculate the P-value for including it in the figure without scientific notation
diff_obj <- survdiff(survival::Surv(OS_days, OS_status) ~ molecular_subtype,  
                     km_result$original_data)
diff_pvalue <- 1 - pchisq(diff_obj$chisq, length(diff_obj$n) - 1)
diff_pvalue_formatted <- format(
  signif(diff_pvalue, 2),
  scientific = FALSE)

km_plot <- survminer::ggsurvplot(fit = km_result$model,
                       data = km_result$original_data,
                       palette = okabe_palette,
                       pval = paste0("P = ", diff_pvalue_formatted),
                       pval.coord = c(400, 0.95),
                       risk.table = TRUE,
                       xlim = c(0, 4000),
                       break.time.by = 500,
                       ggtheme = ggpubr::theme_pubr(),
                       legend = "right",
                       xlab = "Time (days)",
                       legend.title = "Molecular subtype",
                       legend.labs = c("DMG, H3 K28", 
                                       "DMG, H3 K28, TP53 activated", 
                                       "DMG, H3 K28, TP53 lost",
                                       "HGG, H3 wildtype",
                                       "HGG, H3 wildtype, TP53 lost"),
                       risk.table.y.text.col = TRUE,
                       risk.table.y.text = FALSE,
                       # table font size
                       fontsize = 3.75
) 
# add grid to plot separately; this allows us to use theme_pubr() above and still have a grid
km_plot_graph <- km_plot$plot + 
  cowplot::background_grid() + 
  # slightly larger legend text
  theme(legend.text = element_text(size = rel(1)))
km_plot_table <- km_plot$table

# Re-combine plot and table
km_final <- km_plot_graph/km_plot_table +
  plot_layout(widths = c(2, 1), 
              heights = c(2, 0.75))


# Save
ggsave(km_output_pdf, km_final, width = 8, height = 4)

# Export CSV for Zenodo upload, from the `$original_data` field already 
#  in the model object
km_result$original_data %>%
  # ensure OS_status is living/deceased, not 1/2:
  # https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/aa929753fca0019294571ec813e6ae7224b1d8b8/analyses/survival-analysis/util/survival_models.R#L97
  dplyr::mutate(OS_status = ifelse(
    OS_status == 1, 
    "LIVING", 
    "DECEASED"
  )) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_csv(fig4h_csv)



