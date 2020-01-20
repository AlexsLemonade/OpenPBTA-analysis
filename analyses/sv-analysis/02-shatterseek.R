# Yang Yang 2020
# This script is for analyzing chromothripsis, using modified shaterseek code.
# 
# Input files:
# 1. independent-specimens.wgs.primary-plus.tsv
# this file is used to choose independent specimens
# 2. pbta-cnv-cnvkit.seg or pbta-cnv-controlfreec.tsv.gz
# CNV file is needed in shatterseek
# 3. BSID_withoutYandM.vcf
# Those files are SV files, which generated from the "01-process-sv-file.R" script
# 
# Output files:
# 1. PBTA_shatterseek_call_minsize6.csv
# 2. PBTA_shatterseek_call_minsize3.csv
# 3. PBTA_merge_sv.csv
# above three files are generated from shatterseek with min.size is 6 or 3, and combine them.
# 4. PBTA_shatterseek_call_minsize6_test.csv
# 5. PBTA_shatterseek_call_minsize3_test.csv
# 6. PBTA_shatterseek_call_merge.csv
# above three files are generated from shatterseek test
# 7. PBTA_chromothripsis_newlink_region.csv
# this file is about the chromothripsis region and chromothripsis-link region of each sample.
# this file will be used in future analysis.

## ===================== Load Packages =====================
library(devtools)

# withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true"),remotes::install_github('parklab/ShatterSeek'))
library(ShatterSeek)
library(readr)
library(plyr)


## ===================== Root Directory =====================
# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))


## =====================  Create A Subdirectory to Hold All The Output Files ===================== 
# save output files in "scratch/sv-shatterseek" and "analyses/sv-analysis/tables/sv-shatterseek"
output_directory <- file.path(root_dir, "scratch","sv-shatterseek")
output_directory_visible <- file.path(root_dir, "analyses","sv-analysis","tables","sv-shatterseek")

# if output directories don't esist, creat them
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}
if (!dir.exists(output_directory_visible)) {
  dir.create(output_directory_visible, recursive = TRUE)
}

## ===================== Load  Independent Specimen List =====================
independent_specimen_list <- read.table(file.path(root_dir,"data","independent-specimens.wgs.primary-plus.tsv"),header = TRUE,sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load CNV File =====================
# since there are 2 cnv files and 1 consensus file, we read  all of them
cnvkit <- read_tsv(file.path(root_dir,"data","pbta-cnv-cnvkit.seg.gz"))
controlfreec <- read_tsv(file.path(root_dir,"data","pbta-cnv-controlfreec.tsv.gz"))
cnvconsensus <- read_tsv(file.path(root_dir,"analyses","copy_number_consensus_call","results","pbta-cnv-consensus.seg"))
#  unify the column names and chr format of the two cnv files
names(controlfreec)[names(controlfreec)=="chr"] <- "chrom"
names(controlfreec)[names(controlfreec)=="start"] <- "loc.start"
names(controlfreec)[names(controlfreec)=="end"] <- "loc.end"
names(controlfreec)[names(controlfreec)=="copy.number"] <- "copy.num"
controlfreec$chrom <- paste0("chr",controlfreec$chrom)
# we use cnvconsensus as cnv file
cnv_file <- cnvkit
# choose independent specimens and remove all chrY
# shatterseek can't recognize chrY, so remove them
cnv_analysis <-  cnv_file[cnv_file$ID %in% bioid & cnv_file$chrom != "chrY",]
# remove "chr", because shatterseek can't recognize it
cnv_analysis$chrom <- gsub("chr","",cnv_analysis$chrom)


## ===================== Run shatterseek in min.size = 6 or 3, combine and ouput the result =====================
# ss6, ss3 and sv are used to merge data from different samples
ss6 <- data.frame()
ss3 <- data.frame()
sv <- data.frame()

# main process
for (i in bioid) {
  # read sv file one by one
  sv_shatterseek <- read.table(file.path("scratch","sv-vcf",paste(i,"_withoutYandM.tsv",sep="")),sep="\t",header=TRUE)
  print(head(sv_shatterseek))
  # sv_shatterseek_original is a file with chrY and ChrM, which will be used later
  sv_shatterseek_original <- read.table(file.path("scratch","sv-vcf",paste(i,".tsv",sep="")),sep="\t",header=TRUE)
  print(head(sv_shatterseek_original))
  # sv_shatterseek_original[,"sample"] <- i
  # 
  # # merge all sv_shatterseek_original, will be used later
  # sv <- rbind(sv,sv_shatterseek_original)
  # 
  # #  read cnv one by one
  # cnv_shatterseek <-  cnv_analysis[cnv_analysis$ID == i,]
  # 
  # # build sv and cnv data frame
  # SV_data <-
  #   SVs(
  #     chrom1 = as.character(sv_shatterseek$chrom1),
  #     pos1 = as.numeric(sv_shatterseek$pos1),
  #     chrom2 = as.character(sv_shatterseek$chrom2),
  #     pos2 = as.numeric(sv_shatterseek$pos2),
  #     SVtype = as.character(sv_shatterseek$SVtype),
  #     strand1 = as.character(sv_shatterseek$strand1),
  #     strand2 = as.character(sv_shatterseek$strand2)
  #   )
  # CN_data <-
  #   CNVsegs(
  #     chrom = as.character(cnv_shatterseek$chrom),
  #     start = cnv_shatterseek$loc.start,
  #     end = cnv_shatterseek$loc.end,
  #     total_cn = cnv_shatterseek$copy.num
  #   )
  # 
  # # run shatterseek in min.size = 6
  # chromothripsis6 <- shatterseek(SV.sample=SV_data,seg.sample=CN_data,min.Size=6)
  # chrss6 <- chromothripsis6@chromSummary
  # chrss6$sample <- i
  # ss6 <- rbind(ss6,chrss6)
  # 
  # # run shatterseek in min.size = 3
  # chromothripsis3 <- shatterseek(SV.sample=SV_data,seg.sample=CN_data,min.Size=3)
  # chrss3 <- chromothripsis3@chromSummary
  # chrss3$sample <- i
  # ss3 <- rbind(ss3,chrss3)
}

# # save output
# write.csv(ss6,file.path(output_directory,"PBTA_shatterseek_call_minsize6.csv"),row.names = F)
# write.csv(ss3,file.path(output_directory,"PBTA_shatterseek_call_minsize3.csv"),row.names = F)
# write.csv(sv,file.path(output_directory,"PBTA_merge_sv.csv"),row.names=F)
# 
# ## ===================== Test and Filter The Combined Results, Then Ouput =====================
# # read files
# parkcall6 <- ss6
# parkcall3 <- ss3
# # parkcall6 <- read.csv(file.path(output_directory,"PBTA_shatterseek_call_minsize6.csv"),header = TRUE)
# # parkcall3 <- read.csv(file.path(output_directory,"PBTA_shatterseek_call_minsize3.csv"),header = TRUE)
# 
# # test min.size = 6
# ft <- parkcall6$pval_fragment_joins
# ent <- parkcall6$chr_breakpoint_enrichment
# ext <- parkcall6$pval_exp_chr
# exct <- parkcall6$pval_exp_cluster
# 
# # fdr test
# parkcall6$ft_fdr <- p.adjust(ft, method = "fdr")
# parkcall6$ent_fdr <- p.adjust(ent, method = "fdr")
# parkcall6$ext_fdr <- p.adjust(ext, method = "fdr")
# parkcall6$exct_fdr <- p.adjust(exct, method = "fdr")
# 
# # sum all intra sv, including DEL, DUP and INV
# parkcall6$total_intra <- parkcall6$number_DEL + parkcall6$number_DUP + parkcall6$number_h2hINV + parkcall6$number_t2tINV
# 
# # test: when total_intra > 5 and ft_fdr > 0.2 and (ent_fdr <= 0.2 or exct_fdr <= 0.2), call6 will be chrss, else call6 will be no
# parkcall6$call6 <- ifelse(((parkcall6$total_intra > 5)&(parkcall6$ft_fdr > 0.2)&(parkcall6$ent_fdr <= 0.2|parkcall6$exct_fdr <= 0.2)),"chrss","no")
# 
# # save parkcall6
# write.csv(parkcall6,file.path(output_directory,"PBTA_shatterseek_call_minsize6_test.csv"),row.names = F)
# 
# # test min.size = 3
# ft <- parkcall3$pval_fragment_joins
# ent <- parkcall3$chr_breakpoint_enrichment
# ext <- parkcall3$pval_exp_chr
# exct <- parkcall3$pval_exp_cluster
# 
# # fdr test
# parkcall3$ft_fdr <- p.adjust(ft, method = "fdr")
# parkcall3$ent_fdr <- p.adjust(ent, method = "fdr")
# parkcall3$ext_fdr <- p.adjust(ext, method = "fdr")
# parkcall3$exct_fdr <- p.adjust(exct, method = "fdr")
# 
# # sum all intra sv, including DEL, DUP and INV
# parkcall3$total_intra <- parkcall3$number_DEL + parkcall3$number_DUP + parkcall3$number_h2hINV + parkcall3$number_t2tINV
# 
# # test: when total_intra > 2 and TRA >3 and ft_fdr > 0.2, call3 will be chrss, else call3 will be no
# parkcall3$call3 <- ifelse(((parkcall3$total_intra > 2)&(parkcall3$number_TRA > 3)&(parkcall3$ft_fdr > 0.2)),"chrss","no")
# # save parkcall3
# write.csv(parkcall3,file.path(output_directory,"PBTA_shatterseek_call_minsize3_test.csv"),row.names = F)
# 
# # combine call3 result to call6, combine call6 result to call3
# parkcall6$call3 <- parkcall3$call3
# parkcall3$call6 <- parkcall6$call6
# # keep col "chrom","start","end","sample","call3" and "call6"
# keeps <- c("chrom","start","end","sample","call3","call6")
# parkcall6 <- parkcall6[keeps]
# parkcall3 <- parkcall3[keeps]
# 
# #  mergecall is a final test file
# mergecall <- rbind(parkcall6[parkcall6$call6 == "chrss"&parkcall6$call3 == "no",],parkcall3[parkcall3$call3 == "chrss",])
# # change colname
# colnames(mergecall) <- c("chr","start","end","sample","call3","call6")
# # save mergecall
# write.csv(mergecall,file.path(output_directory,"PBTA_shatterseek_call_merge.csv"),row.names=F)
# 
# ## ===================== Call Linking Groups ===================== 
# # read chromothripsis table
# chrss <- mergecall
# # chrss <- read.csv(file.path(output_directory,"PBTA_shatterseek_call_merge.csv"),header=T,stringsAsFactors = F)
# 
# # assign an unique id for each chrss chr, will be used afterwards
# #  format_id is an unique tag
# chrss$format_id <- paste(chrss$sample,chrss$chr,sep = "_")
# chrss$size <- chrss$end-chrss$start
# chrss$chrss_status <- "chrss"
# # chrss_chr is the chrome of original chr
# chrss$chrss_chr <- chrss$chr
# 
# # read sv file, will be used to figure out link-regions
# svlist <- read.csv(file.path(output_directory,"PBTA_merge_sv.csv"),header=T,stringsAsFactors = F)
# svlist$start1 <- svlist$pos1-1
# svlist$end1 <- svlist$pos1
# svlist$start2 <- svlist$pos2-1
# svlist$end2 <- svlist$pos2
# svlist$size <- svlist$end2-svlist$end1
# 
# # if linked chr is chrss, allow expand to 1sd
# # ss is a data.frame used to save chrss-link region and their chrss_chr
# ss <- data.frame()
# format_id <- chrss$format_id
# 
# # chrsslink function
# chrsslink <- function(chrss,formatid=format_id){
#   for (i in 1:nrow(chrss)){ 
#     # [inter1, inter2] is the chrss range
#     inter1 <- max(chrss$start[i]-10000,0)
#     inter2 <- chrss$end[i]+10000
#     # select i'sv from svlist
#     sv <- svlist[svlist$sample==chrss$sample[i],]
#     # svtra is a list of translocation, either side falls chrss range
#     svtra <- sv[((sv$chrom1==chrss$chr[i]&sv$pos1>=inter1&sv$pos1<=inter2)|(sv$chrom2==chrss$chr[i]&sv$pos2>=inter1&sv$pos2<=inter2))&(sv$SVtype=="TRA"),]
#     # if svstr is not empty
#     if (length(svtra$chrom1)>0){
#       # new  cols "chr" "tra_pos1" are  from another chrome 
#       svtra$chr <- ifelse(svtra$chrom1==chrss$chr[i],svtra$chrom2,svtra$chrom1)
#       svtra$tra_pos1 <- ifelse(svtra$chrom1==chrss$chr[i],svtra$pos2,svtra$pos1)
#       # new  col "chrss_chr" is from chrss chrome
#       svtra$chrss_chr <- ifelse(svtra$chrom1==chrss$chr[i],svtra$chrom1,svtra$chrom2)
#       # new  col "n" is the sum  by  chr
#       svtra <- ddply(svtra, .(chr), transform, n = length(chr))
#       # svtramin is the min position by chr
#       svtramin <- aggregate(svtra$tra_pos1, by = list(svtra$chr), min)
#       colnames(svtramin) <- c("chr","start")
#       # svtramax is the max position by chr
#       svtramax <- aggregate(svtra$tra_pos1+1, by = list(svtra$chr), max)
#       colnames(svtramax) <- c("chr","end")
#       #  merge svtramin and svtramax into svtra2
#       svtra2 <- merge(svtramin,svtramax,by.x=("chr"),by.y=("chr"))
#       # merge svtra2 and svtra into svtra3 
#       svtra <- svtra[!duplicated(svtra[,c('chr')]),]
#       svtra3 <- merge(x = svtra2, y = svtra[ , c("chr","chrss_chr", "sample","n")], by = "chr", all.x=TRUE)
#       ss <- rbind(ss,svtra3)
#     }
#   }
#   ss$format_id <- paste(ss$sample,ss$chr,sep="_")
#   idlist <- unique(ss$format_id)
#   
#   # linkss is used to save data
#   linkss <- data.frame(matrix(ncol = 7, nrow = 1))
#   colnames(linkss) <- colnames(ss)
#   # linkss2 collects linkss
#   linkss2 <- data.frame()
#   
#   #  judge if a link-chrss is a chrss region or link region
#   for (id in idlist){
#     # formatid=format_id
#     ss_id <- ss[ss$format_id==id,]
#     # requires at least two chrss linked tra, so sum(ss_id$n)>1
#     if(sum(ss_id$n)>1){
#       # if the chromosome is not an original chrss
#       if(!(id %in% formatid)){
#         # caculate the region for the link region
#         linkss$chr <- unique(ss_id$chr)
#         linkss$start <- min(ss_id$start,ss_id$end)
#         linkss$end <- max(ss_id$start,ss_id$end)
#         linkss$chrss_chr <- paste(c(chrss[chrss$format_id==id,]$chrss_chr,paste(c(ss_id$chrss_chr),collapse="_",sep="_")),collapse="_",sep="_")
#         linkss$sample <- unique(ss_id$sample)
#         linkss$format_id <- unique(ss_id$format_id)
#         linkss$chrss_status <- "no"
#         linkss2 <- rbind(linkss2,linkss)
#       }
#       # if the chromosome is an original chrss
#       if((id %in% formatid)){
#         # caculate the region for the link region
#         linkss$chr <- unique(ss_id$chr)
#         # assign start and end based on TRA position
#         linkss$start <- min(ss_id$start,ss_id$end)
#         linkss$end <- max(ss_id$start,ss_id$end)
#         linkss$chrss_chr <- paste(c(chrss[chrss$format_id==id,]$chrss_chr,paste(c(ss_id$chrss_chr),collapse="_",sep="_")),collapse="_",sep="_")
#         linkss$sample <- unique(ss_id$sample)
#         linkss$format_id <- unique(ss_id$format_id)
#         linkss$chrss_status <- "chrss"
#         linkss2=rbind(linkss2,linkss)
#       }
#     }
#   }
#   
#   # assign start and end based on chrss cluster
#   linkss2 <- merge(linkss2,chrss[,c("format_id","start","end")],by.x=("format_id"),by.y=("format_id"),all.x=T)
#   linkss2$start <- ifelse(!is.na(linkss2$start.y),ifelse(linkss2$start.x<linkss2$start.y,linkss2$start.x,linkss2$start.y),linkss2$start.x)
#   linkss2$end <- ifelse(!is.na(linkss2$end.y),ifelse(linkss2$end.x>linkss2$end.y,linkss2$end.x,linkss2$end.y),linkss2$end.x)
#   keeps <- c("chr","start","end","chrss_chr","sample","format_id","chrss_status")
#   linkss2 <- linkss2[keeps]
#   # merge with original chrss
#   chrss2 <- chrss[keeps]
#   linkss2 <- rbind(linkss2,chrss2)
#   idlist <- unique(linkss2$format_id)
#   
#   linkss <- data.frame(matrix(ncol = 7, nrow = 1))
#   colnames(linkss) <- colnames(linkss2)
#   linkss3 <- data.frame()
#   for (id in idlist){
#     ss_id <- linkss2[linkss2$format_id==id,]
#     # caculate the region for the link region
#     linkss$chr <- unique(ss_id$chr)
#     linkss$start <- min(ss_id$start,ss_id$end)
#     linkss$end <- max(ss_id$start,ss_id$end)
#     linkss$chrss_chr <- paste(c(ss_id$chr,ss_id$chrss_chr),collapse="_",sep="_")
#     linkss$sample <- unique(ss_id$sample)
#     linkss$format_id <- unique(ss_id$format_id)
#     linkss$chrss_status <- unique(linkss2[linkss2$format_id==id,]$chrss_status)
#     linkss3 <- rbind(linkss3,linkss)
#   }
#   return(linkss3)   
# }
# 
# 
# # link region analysis
# chrsslink1 <- chrsslink(chrss)
# chrsslink2 <- chrsslink(chrsslink1)
# 
# # selfrenew until identical
# round <- 0
# while (!identical(chrsslink1,chrsslink2)){
#   chrsslink1 <- chrsslink(chrsslink2)
#   chrsslink2 <- chrsslink(chrsslink1)
#   chrsslink1$chrss_chr <- sapply(strsplit(chrsslink1$chrss_chr, "_", fixed = TRUE), function(x) 
#     paste(unique(x), collapse = "_"))
#   chrsslink2$chrss_chr <- sapply(strsplit(chrsslink2$chrss_chr, "_", fixed = TRUE), function(x) 
#     paste(unique(x), collapse = "_"))
#   round <- round+1
#   print(round)
#   print(nrow(chrsslink1)-nrow(chrsslink2))
# }
# 
# # split and unique value in chrss_chr
# chrsslink1$chr <- gsub("X","23",chrsslink1$chr)
# chrsslink1$chrss_chr <- gsub("X","23",chrsslink1$chrss_chr)
# chrsslink1$chr <- gsub("Y","24",chrsslink1$chr)
# chrsslink1$chrss_chr <- gsub("Y","24",chrsslink1$chrss_chr)
# chrsslink1$chr <- gsub("M","25",chrsslink1$chr)
# chrsslink1$chrss_chr <- gsub("M","25",chrsslink1$chrss_chr)
# chrsslink1$link_chromosome <- 0
# chrsslink1$chrss_chr <- strsplit(as.character(chrsslink1$chrss_chr), "_", fixed = TRUE)
# for (j in 1:length(chrsslink1$sample)){
#   id <- chrsslink1$sample[j]
#   chr <- chrsslink1$chr[j]
#   chrssid <- chrsslink1[chrsslink1$sample==id,]
#   # sort chrss_chr 
#   link_chrss <- as.character(sort(as.numeric(unique(unlist(c(chrssid[chrssid$chr==chr,]$chrss_chr))))))
#   # find intersect rows with link_chrss
#   link_tra <- as.character(sort(as.numeric(unique(unlist(c(chrssid$chrss_chr[unlist(lapply(chrssid$chrss_chr, function(x,y=link_chrss) length(intersect(y, x))>0))]))))))
#   round <- 0
#   while (!identical(link_chrss,link_tra)){
#     link_chrss = link_tra=as.character(sort(as.numeric(unique(unlist(c(chrssid$chrss_chr[unlist(lapply(chrssid$chrss_chr, function(x,y=link_tra) length(intersect(y, x))>0))]))))))
#     link_tra = link_tra=as.character(sort(as.numeric(unique(unlist(c(chrssid$chrss_chr[unlist(lapply(chrssid$chrss_chr, function(x,y=link_chrss) length(intersect(y, x))>0))]))))))
#     round <- round+1
#     print(round)}
#   link_chr <- paste(link_chrss,collapse="_")
#   chrsslink1$link_chromosome[j] <- link_chr
# }
# chrsslink1$chrss_chr <- sapply(chrsslink1$chrss_chr, function(x) 
#   paste(unique(x), collapse = "_"))
# chrsslink1$chr <- gsub("23","X",chrsslink1$chr)
# chrsslink1$chr <- gsub("24","Y",chrsslink1$chr)
# chrsslink1$chrss_chr <- gsub("23","X",chrsslink1$chrss_chr)
# chrsslink1$link_chromosome <- gsub("23","X",chrsslink1$link_chromosome)
# chrsslink1$chrss_chr <- gsub("24","Y",chrsslink1$chrss_chr)
# chrsslink1$link_chromosome <- gsub("24","Y",chrsslink1$link_chromosome)
# chrsslink1$chrss_chr <- gsub("25","M",chrsslink1$chrss_chr)
# chrsslink1$link_chromosome <- gsub("25","M",chrsslink1$link_chromosome)
# 
# # add some column for group plotting
# chrsslink1$cluster_id <- paste(chrsslink1$sample,chrsslink1$link_chromosome,sep="_")
# chrsslink1$chrss_status <- gsub("no","link",chrsslink1$chrss_status)
# chrsslink1$circos_plot <- paste(chrsslink1$sample,".png",sep="")
# chrsslink1$cnv <- paste(chrsslink1$sample,".txt",sep="")
# 
# write.csv(chrsslink1,file.path(output_directory_visible,"PBTA_chromothripsis_newlink_region.csv"),row.names = F)
# write.csv(chrsslink1,file.path(output_directory,"PBTA_chromothripsis_newlink_region.csv"),row.names = F)
