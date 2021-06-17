`%>%` <- dplyr::`%>%`


csv_files <- list.files(path = "./", pattern = "*.csv")

tibble::tibble(filename=csv_files) %>%
  dplyr::mutate(file_contents=purrr::map(filename,
                           ~ readr::read_csv(file.path("./", .)))) %>%
  tidyr::unnest(file_contents) %>%
  dplyr::mutate(signature = stringr::str_replace(filename, ".csv", "")) %>%
  dplyr::select(-filename) %>%
  tidyr::spread(type, probability) -> dat
  
dat %>%
  dplyr::select(-signature) %>% 
  as.matrix() -> refsig_cns


rownames(refsig_cns) <- dat$signature

readr::write_rds(refsig_cns, "refsig_cns_matrix.RDS")
