# S. Spielman for ALSF CCDL 2022
#
# Makes a pdf panel for the HGG Kaplan-Meier survival analysis for Figure TBD

library(tidyverse)
library(survival)
library(patchwork)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# TBD!!!!!!!!!!!!!!!!!!
# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs")#, "fig4", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Input directory
input_dir <- file.path(root_dir, "analyses", "survival-analysis", "results", "subtypes")



## Input and output files ------------------------------------------------
km_result <- read_rds(
  file.path(input_dir,
            "logrank_hgg_subtypes.RDS"
  )
)
km_output_pdf <- file.path(output_dir, "km_hgg.pdf")


# re-factor the km data so legend order matches line order
df <- km_result$original_data
df$molecular_subtype <- factor(df$molecular_subtype, 
                                 levels = c("HGG, H3 wildtype",
                                            "HGG, H3 wildtype, TP53 loss",
                                            "DMG, H3 K28",
                                            "DMG, H3 K28, TP53 activated",
                                            "DMG, H3 K28, TP53 loss"))

# Rebuild the model with this reordered df
fit_factored <- survival::survfit(
  formula(survival::Surv(time = OS_days, event = OS_status) ~ molecular_subtype),
  data = df
)

# Establish color-blind friendly colors
okabe_palette <- colorblindr::palette_OkabeIto[1:5]

km_plot <- survminer::ggsurvplot(fit = fit_factored,
                       data = df,
                       palette = okabe_palette,
                       pval = "P = 0.0002", # P=2e-4, which will look better I think spelled out?
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
                       risk.table.y.text = FALSE
) 
# add grid to plot separately; this allows us to use theme_pubr() above and still have a grid
km_plot_graph <- km_plot$plot + cowplot::background_grid()
km_plot_table <- km_plot$table

# Re-combine plot and table
km_final <- km_plot_graph/km_plot_table +
  plot_layout(widths = c(2, 1), 
              heights = c(2, 0.75))


# Save
ggsave(km_output_pdf, km_final, width = 12, height = 6)

