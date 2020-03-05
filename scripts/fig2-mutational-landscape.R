# Mutational Landscape Figure Run
#
# 2020 
# C. Savonen for ALSF - CCDL
#
# Purpose  Run all steps needed for mutational-landscape Figure.  

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

########## Declare relative file paths
snv_callers_dir <- file.path(root_dir, 
                             "analyses", 
                             "snv-callers")

# 
pbta_consensus_script <- file.path(snv_callers_dir, "run_caller_consensus_analysis-pbta.sh")
tcga_consensus_script <- file.path(snv_callers_dir, "run_caller_consensus_analysis-tcga.sh")

# Plot file paths for this figure
tmb_cdf_plot_file <- file.path(root_dir, "analyses", "tmb-compare-tcga", "plots", "tmb-cdf-pbta-tcga.png")
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


# Run both SNV caller consensus scripts
system(paste("bash", pbta_consensus_script))
system(paste("bash", tcga_consensus_script))

# Run mutational signatures analysis
rmarkdown::render('analyses/mutational-signatures/mutational_signatures.Rmd', 
                  clean = TRUE)

# Run tmb-compare analysis 
rmarkdown::render('analyses/tmb-compare-tcga/compare-tmb.Rmd', 
                  clean = TRUE)

# Read in each file
tmb_cdf_plot <- ggplot2::ggplot()+
  ggpubr::background_image(
    png::readPNG(tmb_cdf_plot_file))

mut_sig_plot_cosmic <- ggplot2::ggplot()+
  ggpubr::background_image(
    png::readPNG(mut_sig_plot_cosmic_file))

mut_sig_plot_nature <- ggplot2::ggplot()+
  ggpubr::background_image(
    png::readPNG(mut_sig_plot_nature_file))

# Arrange plots in grid
ggpubr::ggarrange(tmb_cdf_plot,
                  ggpubr::ggarrange(mut_sig_plot_cosmic, 
                                    mut_sig_plot_nature, 
                                    ncol = 2, labels = c("B", "C")),
                  nrow = 2, 
                  labels = "A")

# Save to PNG
ggplot2::ggsave(file.path("figures", "fig2-mutational-landscapes.png"), 
                width = 17, height = 15,
                units = "in")
