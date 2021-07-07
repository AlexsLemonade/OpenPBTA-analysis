# JN Taroni for ALSF CCDL 2021
#
# In this script we compile pathology diagnosis and pathology free text 
# diagnosis terms/strings used as part of inclusion or exclusion criteria for
# LGAT subtyping 
#
# USAGE: Rscript --vanilla 00-LGAT-select-pathology-dx.R

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# JSON file containing the terms/string
output_file <- file.path(root_dir,
                         "analyses",
                         "molecular-subtyping-LGAT",
                         "lgat-subset",
                         "lgat_subtyping_path_dx_strings.json")


# Inclusion criteria from https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/753#issuecomment-697008356
include_path_dx <- stringr::str_to_lower(
  c(
    "Low-grade glioma/astrocytoma",
    "Ganglioglioma"
  ))

# Exclusion criterion from https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/753#issuecomment-697008356
exclude_path_dx <- stringr::str_to_lower(
  c(
    "Dysembryoplastic neuroepithelial tumor"
  ))

# Update:Recode criteria on the basis of pathology_free_text_diagnosis  
# We were removing these as per https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/995
# but we want to keep these with recode subtypes as GNT now
recode_path_free_text <- stringr::str_to_lower(
  c(
    "desmoplastic infantile astrocytoma",
    "glioneuronal"  # This also covers the more specific cases (e.g., rosette forming glioneuronal tumor)
  ))

# Create a list with the strings we'll use for inclusion/exclusion
terms_list <- list(include_path_dx = include_path_dx,
                   exclude_path_dx = exclude_path_dx,
                   recode_path_free_text = recode_path_free_text)

# Write to file
writeLines(jsonlite::prettify(jsonlite::toJSON(terms_list)), output_file)
