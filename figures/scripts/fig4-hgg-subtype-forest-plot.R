# S. Spielman for ALSF CCDL 2022
#
# Makes a pdf panel of forest plot of survival analysis on MB samples 
#  with immune cell fractions and PDL-1 expression predictors


library(survival) # needed to parse model output
library(tidyverse)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig4", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Directory with result data to input for plot
input_dir <- file.path(root_dir, "analyses", "survival-analysis", "results")


## Input and output files -------------------------------------
survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_subtype.RDS"
  )
)

forest_pdf <- file.path(output_dir, "forest_hgg_subtypes.pdf")

## Make plot --------------------------------------------------

# Set up ordering and labels for y-axis
term_order <- rev(c("molecular_subtypeHGG, H3 wildtype",
                    "molecular_subtypeHGG, H3 wildtype, TP53 loss",
                    "molecular_subtypeDMG, H3 K28",
                    "molecular_subtypeDMG, H3 K28, TP53 activated",
                    "molecular_subtypeDMG, H3 K28, TP53 loss"
                    ))

term_labels <- rev(c("HGG - H3 wildtype",
                     "HGG - H3 wildtype, TP53 loss",
                     "DMG - H3 K28 (reference)",
                     "DMG - H3 K28, TP53 activated",
                     "DMG - H3 K28, TP53 loss"))


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add DMG, H3 K28 as reference
  add_row(term = "molecular_subtypeDMG, H3 K28", 
          estimate = 0) %>% 
  mutate(estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high), 
         # significance indicator column for filling points.
         significant = p.value <= 0.05,
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
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("white", "black"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
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
    plot.subtitle = element_text(face = "bold")
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
  #remove values for reference
  mutate(value = ifelse(term == "DMG - H3 K28 (reference)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                             ",
                      "HR (95% CI)                              P-value")
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
ggsave(forest_pdf, forest_panels, width = 15, height = 4)



