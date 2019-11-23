# Make full SQL join table
#
# Josh Shapiro and Candace Savonen 
# CCDL for ALSF 2019
#
# Because DBI cannot do a full join, this is how we combine all the callers' data
# into one table. 

# Join Strelka with Lancet data
# Union of two left joins
union_strelka_lancet <- dplyr::union_all(
  dplyr::left_join(lancet, strelka,
                   by = join_cols,
                   suffix = c("_lancet", "_strelka")) %>%
    # Want to keep the VAFs from each caller
    dplyr::select(join_cols, "VAF_lancet", "VAF_strelka"),
  dplyr::left_join(strelka, lancet,
                   by = join_cols,
                   suffix = c("_strelka", "_lancet")) %>%
    dplyr::select(join_cols, "VAF_lancet", "VAF_strelka") %>% 
    # VAFs for lancet will already be included in the previous left_join
    dplyr::filter(is.na(VAF_lancet))
) 

# Join Mutect with VarDict data
union_mutect_vardict <- dplyr::union_all(
  dplyr::left_join(vardict, mutect,
                   by = join_cols,
                   suffix = c("_vardict", "_mutect")) %>%
    dplyr::select(join_cols, "VAF_vardict", "VAF_mutect"),
  dplyr::left_join(mutect, vardict,
                   by = join_cols,
                   suffix = c("_mutect", "_vardict")) %>%
    dplyr::select(join_cols, "VAF_vardict", "VAF_mutect") %>% 
    dplyr::filter(is.na(VAF_vardict))
) 

# Join both combinations made above into one big data set
all_caller <- dplyr::union_all(
  dplyr::left_join(union_strelka_lancet, union_mutect_vardict,
                   by = join_cols) %>%
    dplyr::select(join_cols, "VAF_lancet", "VAF_strelka", "VAF_vardict", "VAF_mutect"),
  dplyr::left_join(union_mutect_vardict, union_strelka_lancet,
                   by = join_cols) %>%
    dplyr::select(join_cols, "VAF_lancet", "VAF_strelka", "VAF_vardict", "VAF_mutect") %>% 
    dplyr::filter(is.na(VAF_strelka) && is.na(VAF_lancet))
) 
