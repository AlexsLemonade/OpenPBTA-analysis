# Mutational Landscape Figure 
#
# 2020 
# C. Savonen for ALSF - CCDL
#
# Purpose  Run all steps needed for mutational-landscape Figure.  

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

######## Declare relative file paths of modules used for this figure ###########
snv_callers_dir <- file.path(root_dir, 
                             "analyses", 
                             "snv-callers")
mut_sig_dir <- file.path(root_dir, 
                             "analyses", 
                             "mutational-signatures")
tmb_dir <- file.path(root_dir, 
                     "analyses",
                     "tmb-compare-tcga")

# Import specialized functions from mutational-signatures
source(file.path(mut_sig_dir, "util", "mut_sig_functions.R"))
source(file.path(tmb_dir, "util", "tmb-plot-function.R"))

# Consensus script file paths
pbta_consensus_script <- file.path(snv_callers_dir, "run_caller_consensus_analysis-pbta.sh")
tcga_consensus_script <- file.path(snv_callers_dir, "run_caller_consensus_analysis-tcga.sh")

# Plot file paths for this figure
tmb_cdf_plot_file <- file.path(root_dir, 
                               "analyses", 
                               "tmb-compare-tcga", 
                               "plots", 
                               "tmb-cdf-pbta-tcga.png")

mut_sig_plot_nature_file <- file.path(root_dir, 
                                 "analyses", 
                                 "mutational-signatures", 
                                 "plots", 
                                 "nature", 
                                 "bubble_matrix_nature_mutation_sig.png")

mut_sig_plot_cosmic_file <- file.path(root_dir, 
                                 "analyses", 
                                 "mutational-signatures", 
                                 "plots", 
                                 "cosmic", 
                                 "bubble_matrix_cosmic_mutation_sig.png")

########################### Run the analyses needed ############################
# Run both SNV caller consensus scripts
system(paste("bash", pbta_consensus_script))
system(paste("bash", tcga_consensus_script))

# Run mutational signatures analysis
rmarkdown::render('analyses/mutational-signatures/mutational_signatures.Rmd', 
                  clean = TRUE)

###################### Re-run the individual plots #############################
# TODO: update if the this tmb consensus file gets updated in a future data release
tmb_pbta <- data.table::fread(file.path(snv_callers_dir,
                                        "results",
                                        "consensus", 
                                        "pbta-snv-mutation-tmb-coding.tsv")) %>% 
  # This variable is weird when binding but we don't need it for the plot so we'll just remove it. 
  dplyr::select(-genome_size) %>% 
  dplyr::filter(experimental_strategy != "Panel")

# TODO: update if this tmb consensus file gets updated in a future data release
tmb_tcga <- data.table::fread(file.path(snv_callers_dir,
                                        "results",
                                        "consensus",
                                        "tcga-snv-mutation-tmb-coding.tsv")) %>% 
  dplyr::select(-genome_size)

# Make PBTA TMB plot
pbta_plot <- tmb_cdf_plot(tmb_pbta, plot_title = "PBTA", colour = "#3BC8A2") +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 12), 
    plot.margin = ggplot2::unit(c(1, 1, -2, 1), "cm")
  )

# Make TCGA plot
tcga_plot <- tmb_cdf_plot(tmb_tcga, plot_title = "TCGA (Adult)", colour = "#630882") +
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    strip.text.x = ggplot2::element_text(size = 11), 
    plot.margin = ggplot2::unit(c(1, 1, -4, 0), "cm")
  ) 

# Set up cosmic signature plots
nature_sigs_df <- readr::read_tsv(file.path(mut_sig_dir, 
                                            "results",
                                            "nature_signatures_results.tsv"))

cosmic_sigs_df <- readr::read_tsv(file.path(mut_sig_dir, 
                                            "results",
                                            "nature_signatures_results.tsv"))

mut_sig_plot_cosmic <- bubble_matrix_plot(cosmic_sigs_df, 
                                          label = "COSMIC") + 
  ggplot2::theme(legend.position = "none",
                 axis.text = ggplot2::element_text(size = 10), 
                 plot.margin = ggplot2::unit(c(.5, 1, .5, .5), "cm")) 

mut_sig_plot_nature <- bubble_matrix_plot(nature_sigs_df, 
                                          label = "Alexandrov et al, 2013") +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 10), 
                 plot.margin = ggplot2::unit(c(.5, .5, .5, .5), "cm")) 

# Arrange plots in grid
ggpubr::ggarrange(ggpubr::ggarrange(pbta_plot, 
                                    tcga_plot, 
                                    ncol = 2, 
                                    labels = c("A", ""), 
                                    widths = c(2, 1), 
                                    font.label = list(size = 22)),
                  ggpubr::ggarrange(mut_sig_plot_cosmic, 
                                    mut_sig_plot_nature, 
                                    ncol = 2, 
                                    labels = c("B", "C"), 
                                    widths = c(1.6, 2), 
                                    font.label = list(size = 22)),
                  nrow = 2, 
                  heights = c(1.5, 2))

# Save to PNG
ggplot2::ggsave(file.path("figures", "fig2-mutational-landscapes.png"), 
                width = 16.5, height = 15,
                units = "in")
