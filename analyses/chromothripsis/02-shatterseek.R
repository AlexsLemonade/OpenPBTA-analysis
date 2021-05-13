# Yang Yang 2020
# This script is for analyzing chromothripsis, using modified ShatterSeek code.
#
# Input files:
# 1. independent-specimens.wgs.primary-plus.tsv
# this file is used to choose independent specimens
# 2. analyses/copy_number_consensus_call/results/pbta-cnv-consensus.seg.gz
# CNV file is needed in shatterseek
# 3. BSID_withoutYandM.tsv and BSID.tsv
# Those files are SV files, which generated from the "01-process-sv-file.R" script
#
# Output files:
# PBTA_chromothripsis_newlink_region.csv
# this file is about the chromothripsis region and chromothripsis-link region of each sample.
# this file will be used in future analysis.
# chrss_status: chrss means this region is an original chromothripsis region, link means this region is an chromothripsis-link region
# link_chromosome: all the chromosomes are linked together, they are a chromothripsis cluster

## ===================== Load Packages =====================
library(devtools)
library(ShatterSeek)
library(readr)
library(plyr)
library(dplyr)


## ===================== Root Directory =====================
# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))


## =====================  Create A Subdirectory to Hold All The Output Files =====================
# save output files in "analyses/sv-analysis/tables/sv-shatterseek"
output_directory <- file.path(root_dir, "analyses","sv-analysis","result")


# if output directories don't esist, creat them
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}


## ===================== Load  Independent Specimen List =====================
independent_specimen_list <- read.table(file.path(root_dir,"data","independent-specimens.wgs.primary-plus.tsv"),header = TRUE,sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load CNV File =====================
# read cnv  consensus file
cnvconsensus <- read_tsv(file.path(root_dir, "data", "pbta-cnv-consensus.seg.gz"))
# For shatterseek can not work with NA copy number, remove rows with NA copy number
cnvconsensus <- cnvconsensus[is.na(cnvconsensus$copy.num)==F,]

# choose independent specimens and remove all chrY
# shatterseek can't recognize chrY, so remove them
cnv_analysis <-  cnvconsensus[cnvconsensus$ID %in% bioid & cnvconsensus$chrom != "chrY",]
# remove "chr", because shatterseek can't recognize it
cnv_analysis$chrom <- gsub("chr","",cnv_analysis$chrom)


## ===================== Run shatterseek in min.size = 6 or 3, combine and ouput the result =====================
# ss6, ss3 and sv are used to merge data from different samples
ss6 <- data.frame()
ss3 <- data.frame()
svlist <- data.frame()

# main process
for (i in bioid) {
  # read sv file one by one
  sv_shatterseek <- read.table(file.path("scratch","sv-vcf",paste(i,"_withoutYandM.tsv",sep="")),sep="\t",header=TRUE)

  # sv_shatterseek_original is a file with chrY and ChrM, which will be used later
  sv_shatterseek_original <- read.table(file.path("scratch","sv-vcf",paste(i,".tsv",sep="")),sep="\t",header=TRUE)

  #  read cnv one by one
  cnv_shatterseek <-  cnv_analysis[cnv_analysis$ID == i,]

  # For CI
  # if we meet empty sv file, jump into next loop
  nsv1 <- nrow(sv_shatterseek)
  nsv2  <-  nrow(sv_shatterseek_original)
  # if nsv1 or nsv2 is 0, jump into next loop
  if (nsv1 == 0 | nsv2 == 0) {
    print(paste0(i," has an empty sv file"))
    next;
  }


  # add sample id
  sv_shatterseek_original$sample <- i
  # merge all sv_shatterseek_original, will be used later
  svlist <- rbind(svlist,sv_shatterseek_original)


  # build sv and cnv data frame
  SV_data <-
    SVs(
      chrom1 = as.character(sv_shatterseek$chrom1),
      pos1 = as.numeric(sv_shatterseek$pos1),
      chrom2 = as.character(sv_shatterseek$chrom2),
      pos2 = as.numeric(sv_shatterseek$pos2),
      SVtype = as.character(sv_shatterseek$SVtype),
      strand1 = as.character(sv_shatterseek$strand1),
      strand2 = as.character(sv_shatterseek$strand2)
    )
  CN_data <-
    CNVsegs(
      chrom = as.character(cnv_shatterseek$chrom),
      start = cnv_shatterseek$loc.start,
      end = cnv_shatterseek$loc.end,
      total_cn = cnv_shatterseek$copy.num
    )

  # run shatterseek in min.size = 6
  chromothripsis6 <- shatterseek(SV.sample=SV_data,seg.sample=CN_data,min.Size=6)
  chrss6 <- chromothripsis6@chromSummary
  chrss6$sample <- i
  ss6 <- rbind(ss6,chrss6)

  # run shatterseek in min.size = 3
  chromothripsis3 <- shatterseek(SV.sample=SV_data,seg.sample=CN_data,min.Size=3)
  chrss3 <- chromothripsis3@chromSummary
  chrss3$sample <- i
  ss3 <- rbind(ss3,chrss3)
}


## ===================== Test and Filter The Combined Results, Then Ouput =====================
# fdr function
add_shatterseek_fdr <- function(shatterseek_result) {
  shatterseek_result <- shatterseek_result %>%
    dplyr::mutate(
      ft_fdr = p.adjust(pval_fragment_joins, method = "fdr"),
      ent_fdr = p.adjust(chr_breakpoint_enrichment, method = "fdr"),
      ext_fdr = p.adjust(pval_exp_chr, method = "fdr"),
      exct_fdr = p.adjust(pval_exp_cluster, method = "fdr")
    )
  return(shatterseek_result)
}

# fdr test in min.size=6
ss6 <- add_shatterseek_fdr(ss6)
# sum all intra sv, including DEL, DUP and INV
ss6$total_intra <- ss6$number_DEL + ss6$number_DUP + ss6$number_h2hINV + ss6$number_t2tINV
# test: when total_intra > 5 and ft_fdr > 0.2 and (ent_fdr <= 0.2 or exct_fdr <= 0.2), call6 will be chrss, else call6 will be no
ss6$call6 <- ifelse(((ss6$total_intra > 5)&(ss6$ft_fdr > 0.2)&(ss6$ent_fdr <= 0.2|ss6$exct_fdr <= 0.2)),"chrss","no")

# fdr test in min.size=3
ss3 <- add_shatterseek_fdr(ss3)
# sum all intra sv, including DEL, DUP and INV
ss3$total_intra <- ss3$number_DEL + ss3$number_DUP + ss3$number_h2hINV + ss3$number_t2tINV
# test: when total_intra > 2 and TRA >3 and ft_fdr > 0.2, call3 will be chrss, else call3 will be no
ss3$call3 <- ifelse(((ss3$total_intra > 2)&(ss3$number_TRA > 3)&(ss3$ft_fdr > 0.2)),"chrss","no")

# combine call3 result to call6, combine call6 result to call3
ss6 <- ss3 %>%
  # pick the relevant columns from ss3 to join to ss6
  select(chrom, `sample`, call3) %>%
  # join by chromosome and sample
  inner_join(ss6, by = c("chrom", "sample"))
ss3 <- ss6 %>%
  # pick the relevant columns from ss6 to join to ss3
  select(chrom, `sample`, call6) %>%
  # join by chromosome and sample
  inner_join(ss3, by = c("chrom", "sample"))


# keep col "chrom","start","end","sample","call3" and "call6"
keeps <- c("chrom","start","end","sample","call3","call6")
ss6 <- ss6[keeps]
ss3 <- ss3[keeps]

#  chrss is a final test file
chrss <- rbind(ss6[ss6$call6 == "chrss"&ss6$call3 == "no",],ss3[ss3$call3 == "chrss",])
# change colname
colnames(chrss) <- c("chr","start","end","sample","call3","call6")


## ===================== Call Linking Groups =====================
# chrsslink function
# this function is for finding chromothripsis-link regions
# parameter chrss is the output of section "Test and Filter The Combined Results"
# parameter formatid is the colcunm "formatid" of chrss, is an unique tag
# the output includes all orginal chromothripsis region and link region
chrsslink <- function(chrss,formatid=format_id){
  for (i in 1:nrow(chrss)){
    # [inter1, inter2] is the chrss range
    inter1 <- max(chrss$start[i]-10000,0)
    inter2 <- chrss$end[i]+10000
    # select i'sv from svlist
    sv <- svlist[svlist$sample==chrss$sample[i],]
    # chrsschri is the chr of this line
    chrsschri <- as.character(chrss$chr[i])
    # svtra is a list of translocation, either side falls chrss range
    svtra <- sv[((sv$chrom1==chrsschri&sv$pos1>=inter1&sv$pos1<=inter2)|(sv$chrom2==chrsschri&sv$pos2>=inter1&sv$pos2<=inter2))&(sv$SVtype=="TRA"),]
    # if svstr is not empty
    if (length(svtra$chrom1)>0){
      # new  cols "chr" "tra_pos1" are  from another chrome
      svtra$chr <- ifelse(svtra$chrom1==chrsschri,svtra$chrom2,svtra$chrom1)
      svtra$tra_pos1 <- ifelse(svtra$chrom1==chrsschri,svtra$pos2,svtra$pos1)
      # new  col "chrss_chr" is from chrss chrome
      svtra$chrss_chr <- ifelse(svtra$chrom1==chrsschri,svtra$chrom1,svtra$chrom2)
      # new  col "n" is the sum  by  chr
      svtra <- ddply(svtra, .(chr), transform, n = length(chr))
      # svtramin is the min position by chr
      svtramin <- aggregate(svtra$tra_pos1, by = list(svtra$chr), min)
      colnames(svtramin) <- c("chr","start")
      # svtramax is the max position by chr
      svtramax <- aggregate(svtra$tra_pos1+1, by = list(svtra$chr), max)
      colnames(svtramax) <- c("chr","end")
      #  merge svtramin and svtramax into svtra2
      svtra2 <- merge(svtramin,svtramax,by.x=("chr"),by.y=("chr"))
      # merge svtra2 and svtra into svtra3
      svtra <- svtra[!duplicated(svtra[,c('chr')]),]
      svtra3 <- merge(x = svtra2, y = svtra[ , c("chr","chrss_chr", "sample","n")], by = "chr", all.x=TRUE)
      ss <- rbind(ss,svtra3)
    }
  }
  ss$format_id <- paste(ss$sample,ss$chr,sep="_")
  idlist <- unique(ss$format_id)

  # linkss is used to save data
  linkss <- data.frame(matrix(ncol = 7, nrow = 1))
  colnames(linkss) <- colnames(ss)
  # linkss2 collects linkss
  linkss2 <- data.frame()

  #  judge if a link-chrss is a chrss region or link region
  for (id in idlist){
    # formatid=format_id
    ss_id <- ss[ss$format_id==id,]
    # requires at least two chrss linked tra, so sum(ss_id$n)>1
    if(sum(ss_id$n)>1){
      # if the chromosome is not an original chrss
      if(!(id %in% formatid)){
        # caculate the region for the link region
        linkss$chr <- unique(ss_id$chr)
        linkss$start <- min(ss_id$start,ss_id$end)
        linkss$end <- max(ss_id$start,ss_id$end)
        linkss$chrss_chr <- paste(c(chrss[chrss$format_id==id,]$chrss_chr,paste(c(ss_id$chrss_chr),collapse="_",sep="_")),collapse="_",sep="_")
        linkss$sample <- unique(ss_id$sample)
        linkss$format_id <- unique(ss_id$format_id)
        linkss$chrss_status <- "no"
        linkss2 <- rbind(linkss2,linkss)
      }
      # if the chromosome is an original chrss
      if((id %in% formatid)){
        # caculate the region for the link region
        linkss$chr <- unique(ss_id$chr)
        # assign start and end based on TRA position
        linkss$start <- min(ss_id$start,ss_id$end)
        linkss$end <- max(ss_id$start,ss_id$end)
        linkss$chrss_chr <- paste(c(chrss[chrss$format_id==id,]$chrss_chr,paste(c(ss_id$chrss_chr),collapse="_",sep="_")),collapse="_",sep="_")
        linkss$sample <- unique(ss_id$sample)
        linkss$format_id <- unique(ss_id$format_id)
        linkss$chrss_status <- "chrss"
        linkss2=rbind(linkss2,linkss)
      }
    }
  }

  # assign start and end based on chrss cluster
  linkss2 <- merge(linkss2,chrss[,c("format_id","start","end")],by.x=("format_id"),by.y=("format_id"),all.x=T)
  linkss2$start <- ifelse(!is.na(linkss2$start.y),ifelse(linkss2$start.x<linkss2$start.y,linkss2$start.x,linkss2$start.y),linkss2$start.x)
  linkss2$end <- ifelse(!is.na(linkss2$end.y),ifelse(linkss2$end.x>linkss2$end.y,linkss2$end.x,linkss2$end.y),linkss2$end.x)
  keeps <- c("chr","start","end","chrss_chr","sample","format_id","chrss_status")
  linkss2 <- linkss2[keeps]
  # merge with original chrss
  chrss2 <- chrss[keeps]
  linkss2 <- rbind(linkss2,chrss2)
  idlist <- unique(linkss2$format_id)

  linkss <- data.frame(matrix(ncol = 7, nrow = 1))
  colnames(linkss) <- colnames(linkss2)
  linkss3 <- data.frame()
  for (id in idlist){
    ss_id <- linkss2[linkss2$format_id==id,]
    # caculate the region for the link region
    linkss$chr <- unique(ss_id$chr)
    linkss$start <- min(ss_id$start,ss_id$end)
    linkss$end <- max(ss_id$start,ss_id$end)
    linkss$chrss_chr <- paste(c(ss_id$chr,ss_id$chrss_chr),collapse="_",sep="_")
    linkss$sample <- unique(ss_id$sample)
    linkss$format_id <- unique(ss_id$format_id)
    linkss$chrss_status <- unique(linkss2[linkss2$format_id==id,]$chrss_status)
    linkss3 <- rbind(linkss3,linkss)
  }
  return(linkss3)
}


# assign an unique id for each chrss chr, will be used afterwards
#  format_id is an unique tag
chrss$format_id <- paste(chrss$sample,chrss$chr,sep = "_")
chrss$size <- chrss$end-chrss$start
chrss$chrss_status <- "chrss"
# chrss_chr is the chrome of original chr
chrss$chrss_chr <- chrss$chr

# read sv file, will be used to figure out link-regions
svlist$start1 <- svlist$pos1-1
svlist$end1 <- svlist$pos1
svlist$start2 <- svlist$pos2-1
svlist$end2 <- svlist$pos2
svlist$size <- svlist$end2-svlist$end1

# if linked chr is chrss, allow expand to 1sd
# ss is a data.frame used to save chrss-link region and their chrss_chr
ss <- data.frame()
format_id <- chrss$format_id

# link region analysis
# chrsslink(chrss) helps chrss to find the linked region —— chrsslink1
# chrsslink(chrsslink1) helps chrsslink1 to find the linked region —— chrsslink2
# we will compare chrsslink1 and chrsslink2 in the next lines (selfrenew)
# if they are the same, means we have found all linked regions
# if they are not the same, we should keep finding more linked regions until chrsslink1 and chrsslink2 are identical
chrsslink1 <- chrsslink(chrss)
chrsslink2 <- chrsslink(chrsslink1)

# selfrenew until identical
round <- 0
while (!identical(chrsslink1,chrsslink2)){
  chrsslink1 <- chrsslink(chrsslink2)
  chrsslink2 <- chrsslink(chrsslink1)
  chrsslink1$chrss_chr <- sapply(strsplit(chrsslink1$chrss_chr, "_", fixed = TRUE), function(x)
    paste(unique(x), collapse = "_"))
  chrsslink2$chrss_chr <- sapply(strsplit(chrsslink2$chrss_chr, "_", fixed = TRUE), function(x)
    paste(unique(x), collapse = "_"))
  # selfrenew time
  round <- round+1
  print(paste("Round:", round))
  # print how many different lines between chrsslink1 and chrsslink2
  # just to see the pace of selfrenew, in case of endless-loop
  print(nrow(chrsslink1)-nrow(chrsslink2))
}

# split and unique value in chrss_chr
chrsslink1$chr <- gsub("X","23",chrsslink1$chr)
chrsslink1$chrss_chr <- gsub("X","23",chrsslink1$chrss_chr)
chrsslink1$chrss_chr <- gsub("Y","24",chrsslink1$chrss_chr)
chrsslink1$chrss_chr <- gsub("M","25",chrsslink1$chrss_chr)
chrsslink1$link_chromosome <- 0
chrsslink1$chrss_chr <- strsplit(as.character(chrsslink1$chrss_chr), "_", fixed = TRUE)
for (j in 1:length(chrsslink1$sample)){
  id <- chrsslink1$sample[j]
  chr <- chrsslink1$chr[j]
  chrssid <- chrsslink1[chrsslink1$sample==id,]
  # sort chrss_chr
  link_chrss <- as.character(sort(as.numeric(unique(unlist(c(chrssid[chrssid$chr==chr,]$chrss_chr))))))
  # find intersect rows with link_chrss
  link_tra <- as.character(sort(as.numeric(unique(unlist(c(chrssid$chrss_chr[unlist(lapply(chrssid$chrss_chr, function(x,y=link_chrss) length(intersect(y, x))>0))]))))))
  round <- 0
  while (!identical(link_chrss,link_tra)){
    link_chrss = link_tra=as.character(sort(as.numeric(unique(unlist(c(chrssid$chrss_chr[unlist(lapply(chrssid$chrss_chr, function(x,y=link_tra) length(intersect(y, x))>0))]))))))
    link_tra = link_tra=as.character(sort(as.numeric(unique(unlist(c(chrssid$chrss_chr[unlist(lapply(chrssid$chrss_chr, function(x,y=link_chrss) length(intersect(y, x))>0))]))))))
    round <- round+1
    paste("Round:", round)}
  link_chr <- paste(link_chrss,collapse="_")
  chrsslink1$link_chromosome[j] <- link_chr
}
chrsslink1$chrss_chr <- sapply(chrsslink1$chrss_chr, function(x)
  paste(unique(x), collapse = "_"))
chrsslink1$chr <- gsub("23","X",chrsslink1$chr)
chrsslink1$chrss_chr <- gsub("23","X",chrsslink1$chrss_chr)
chrsslink1$link_chromosome <- gsub("23","X",chrsslink1$link_chromosome)
chrsslink1$chrss_chr <- gsub("24","Y",chrsslink1$chrss_chr)
chrsslink1$link_chromosome <- gsub("24","Y",chrsslink1$link_chromosome)
chrsslink1$chrss_chr <- gsub("25","M",chrsslink1$chrss_chr)
chrsslink1$link_chromosome <- gsub("25","M",chrsslink1$link_chromosome)

# add some column for group plotting
chrsslink1$cluster_id <- paste(chrsslink1$sample,chrsslink1$link_chromosome,sep="_")
chrsslink1$chrss_status <- gsub("no","link",chrsslink1$chrss_status)


write.csv(chrsslink1,file.path(output_directory,"PBTA_chromothripsis_newlink_region.csv"),row.names = F)
