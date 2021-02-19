# functions takes in a maf format dataframe and a variant_type to return an upset plot
getUpset<-function(maf,variant_type="SNP"){
  # Add up how many callers identify eac mutation
  caller_mat <- maf %>% 
    dplyr::filter(Variant_Type %in% variant_type) %>%  
    dplyr::select("Chromosome",
                  "Start_Position",
                  "Tumor_Sample_Barcode", 
                  "Variant_Classification",
                  "caller") %>%
    unique() %>%
    reshape2::dcast(Chromosome+Start_Position+Tumor_Sample_Barcode+Variant_Classification ~ caller,fun.aggregate = length) %>%
    dplyr::select(-c("Chromosome","Start_Position", "Tumor_Sample_Barcode", "Variant_Classification"))
  
  UpSetR::upset(data = caller_mat ,
                order.by = "freq",main.bar.color = "gray")
}



