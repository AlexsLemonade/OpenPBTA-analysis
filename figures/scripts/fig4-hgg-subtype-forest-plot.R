# S. Spielman for ALSF CCDL & Jo Lynne Rokita for D3b, 2022-3
#
# Makes a pdf panel of forest plot of survival analysis on HGG samples 
#  with molecular subtype as predictors


library(survival) # needed to parse model output
library(tidyverse)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig4", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Zenodo CSV output directory and file path
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
fig4g_csv <- file.path(zenodo_tables_dir, "figure-4g-data.csv")


# Directory with result data to input for plot
input_dir <- file.path(root_dir, "analyses", "survival-analysis", "results", "subtypes")


## Input and output files -------------------------------------
# here, $model has the model and $data has the associated data
survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_subtype.RDS"
  )
)

forest_pdf <- file.path(output_dir, "forest_hgg_subtypes.pdf")

## Make plot --------------------------------------------------
ref_term <- "molecular_subtypeDMG, H3 K28"

# Set up ordering and labels for y-axis
term_order <- rev(c("molecular_subtypeHGG, H3 wildtype",
                    "molecular_subtypeHGG, H3 wildtype, TP53 loss",
                    "molecular_subtypeDMG, H3 K28, TP53 activated",
                    "molecular_subtypeDMG, H3 K28, TP53 loss",
                    ref_term))

term_labels <- rev(c("HGG - H3 wildtype",
                     "HGG - H3 wildtype, TP53 lost",
                     "DMG - H3 K28, TP53 activated",
                     "DMG - H3 K28, TP53 lost",
                     "DMG - H3 K28 (ref)"))


# Get n and event info from glance output
survival_n <- broom::glance(survival_result$model) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result$model) %>%
  # Add DMG, H3 K28 as reference
  add_row(term = ref_term, 
          estimate = 0) %>% 
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
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
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
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
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
    plot.subtitle = element_text(face = "bold"), 
    # thinner axes, ticks for compilation
    axis.line = element_line(size = rel(0.25)),
    axis.ticks = element_line(size = rel(0.25))
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
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
  #remove values for reference
  mutate(value = ifelse(term == "DMG - H3 K28 (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                ",
                      "HR (95% CI)         P-value")
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
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(forest_pdf, forest_panels, width = 10, height = 3)


# Export CSV for Zenodo upload
survival_result$data %>%
  # order columns and only keep columns that are in the model:
  # OS_status and OS_years
  # ~ molecular_subtype
  dplyr::select(Kids_First_Biospecimen_ID, 
                OS_years, 
                OS_status,
                molecular_subtype) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_csv(fig4g_csv)

