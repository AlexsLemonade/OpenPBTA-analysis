# In this script we will be gathering pathology diagnosis
# and pahology free text diagnosis terms to select HGG
# samples for downstream HGG subtyping analysis and save 
# the json file in hgg-subset folder

output_dir <- "hgg-subset"
output_file <- file.path(output_dir, 
                         "hgg_subtyping_path_dx_strings.json")

# The `pathology_diagnosis` fields for HGG
# as we identified in 00-v17-HGG-select-pathology-dx.Rmd are:
exact_path_dx<- c(
  "High-grade glioma/astrocytoma (WHO grade III/IV)",
  "Brainstem glioma- Diffuse intrinsic pontine glioma"
)

# Gliomatosis Cerebri can be high grade glioma or low grade 
# glioma so we will add an inclusion criteria for v18 release 
# to only keep `Gliomatosis Cerebri` samples if pathology_free_text_diagnosis
# as `anaplastic gliomatosis cerebri (who grade 4)`
gliomatosis_path_free_text_exact = "anaplastic gliomatosis cerebri (who grade 4)"

# Create a list with the strings we'll use for inclusion.
terms_list <- list(exact_path_dx = exact_path_dx,
                   gliomatosis_path_free_text_exact = gliomatosis_path_free_text_exact)


#Save this list as JSON.
writeLines(jsonlite::prettify(jsonlite::toJSON(terms_list)), output_file)

