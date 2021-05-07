# In this script we will be gathering pathology diagnosis
# and pathology free text diagnosis terms to select embryonal
# samples for downstream embryonal subtyping analysis and save 
# the json file in subset-files folder

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

output_file <- file.path(root_dir,
                         "analyses",
                         "molecular-subtyping-embryonal",
                         "subset-files",
                         "embryonal_subtyping_path_dx_strings.json")

# The `pathology_diagnosis` fields for emmbryonal
# as we identified in 00-v17-embryonal-select-pathology-dx.Rmd are:
path_dx_terms<- c(
  "Supratentorial or Spinal Cord PNET",
  "Embryonal Tumor with Multilayered Rosettes"
)


free_text_dx_terms <- c(
  "embryonal tumor with multilayer rosettes, ros (who grade iv)",
  "embryonal tumor, nos, congenital type",
  "ependymoblastoma",
  "medulloepithelioma",
  "medullooepithelioma" # It was noted that this misspelling was included in the `pathology_free_text_diagnosis` column of the histologies file so we will add it to the JSON file [per this comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/788#issuecomment-700717514)
)

terms_list <- list(include_path_dx = path_dx_terms,
                   include_free_text = free_text_dx_terms)

#Save this list as JSON.
writeLines(jsonlite::prettify(jsonlite::toJSON(terms_list)), output_file)
