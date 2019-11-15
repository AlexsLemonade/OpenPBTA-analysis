# this script is  for plotting SV distribution figures
# 
# Input files:
# 1.independent-specimens.wgs.primary-plus.tsv
# this file is used to choose independent specimens
# 2. BSID.vcf
# this is the sv file generated from script "01-process-sv-file.R"
# 3. pbta-histologies.tsv
# this is the patient information file
# 4. PBTA_chromothripsis_newlink_region.csv
# this file was generated from script "02-shatterseek.R"
# 
# Output files:
# 1. figure1: SV proportion
# 2. figure2: SV number across tumors
# 3. figure3: SV chromothripsis percentage

library(ggplot2)
library(reshape2)
library(ggthemes)
library("Rmisc")
library("plyr")
library("cowplot")
library(magick)
library(grid)
library(randomcoloR)

## ===================== Load  Independent specimen list =====================
independent_specimen_list <- read.table("independent-specimens.wgs.primary-plus.tsv",header = TRUE,sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Generate data.frames =====================
stats <- data.frame(row.names =bioid )
stats$DUP <- 0
stats$DEL <- 0
stats$INV <- 0
stats$TRA <- 0
stats$SUM <- 0
stats$LOGSUM <- 0
# read sv files one by one
for (i in 1:length(bioid)) {
  id <- bioid[i]
  file <- read.delim(paste("D:/Project/Signature/PBTA/V7/sv/",id,".vcf",sep = ""))
  stats$DUP[i] <- nrow(file[file$SVtype == "DUP",])
  stats$DEL[i] <- nrow(file[file$SVtype == "DEL",])
  stats$INV[i] <- nrow(file[(file$SVtype == "h2hINV" |file$SVtype == "t2tINV"),])
  stats$TRA[i] <- nrow(file[file$SVtype == "TRA",])
  stats$SUM[i] <- sum(stats[i,c(1,2,3,4)])
  stats$LOGSUM[i] <-sum(log10(stats[i,1]),log10(stats[i,2]),log10(stats[i,3]),log10(stats[i,4]))
}

#  melt stats
stats_melt <- melt(stats[,c(1,2,3,4)],variable.name="Type",value.name="Number")
stats_melt$Name <- rep(bioid,4)
stats_melt$Log <- log10(stats_melt$Number)
stats_melt$Sum <- rep(stats$SUM,4)
stats_melt$Logsum <- rep(stats$LOGSUM,4)



## ===================== Imfort Patient Information and Order =====================
pbta_histologies <- read.table("pbta-histologies.tsv",header = TRUE,sep = "\t")
pbta_chrss <- read.csv("PBTA_chromothripsis_newlink_region.csv",header = TRUE)
info <- pbta_histologies[pbta_histologies$Kids_First_Biospecimen_ID %in% bioid,]

# Chrss col can recode chromothripsis sample
info$Chrss <- 0
info[info$Kids_First_Biospecimen_ID %in% pbta_chrss$sample,]$Chrss <- 1

# stats_melt_info
stats_melt_info<-merge(stats_melt,info,by.x = "Name",by.y = "Kids_First_Biospecimen_ID")
# SampleNumber is the total number of specific tumor
stats_melt_info$SampleNumber <- 0
# Disease type come from disease_type_new, and change diseases whose total number < 10 to "Others"
stats_melt_info$disease_type_update <- "NA"
for (i in unique(stats_melt_info$disease_type_new)) {
  stats_melt_info[stats_melt_info$disease_type_new == i,]$SampleNumber = summary(stats_melt_info$disease_type_new)[i]/4
  if (stats_melt_info[stats_melt_info$disease_type_new == i,]$SampleNumber < 10) {
    stats_melt_info[stats_melt_info$disease_type_new == i,]$disease_type_update = "Others"
  }  else {
    stats_melt_info[stats_melt_info$disease_type_new == i,]$disease_type_update = i
  }
}

# stats_info
stats_info <- stats
stats_info$Name <- row.names(stats)
stats_info<-merge(stats_info,info,by.x = "Name",by.y = "Kids_First_Biospecimen_ID")
stats_info$SampleNumber <- 0
stats_info$disease_type_update <- "NA"
for (i in unique(stats_info$disease_type_new)) {
  stats_info[stats_info$disease_type_new == i,]$SampleNumber = summary(stats_info$disease_type_new)[i]
  if (stats_info[stats_info$disease_type_new == i,]$SampleNumber < 10) {
    stats_info[stats_info$disease_type_new == i,]$disease_type_update = "Others"
  }  else {
    stats_info[stats_info$disease_type_new == i,]$disease_type_update = i
  }
}
stats_info$SampleNumberUpdate = 0
for (i in unique(stats_info$disease_type_update)) {
  stats_info[stats_info$disease_type_update == i,]$SampleNumberUpdate = summary(stats_info$disease_type_update)[i]
}


# order by and sum
order <- order(stats_melt_info$SampleNumber,stats_melt_info$disease_type_update,stats_melt_info$Sum,decreasing = TRUE)
order2 <- order(stats_info$SampleNumber,stats_info$disease_type_update,stats_info$SUM,decreasing = TRUE)

# do order
stats_melt_info <- stats_melt_info[order,]
stats_info <- stats_info[order2,]

## ===================== PLOT =====================
# 1. plot "SV proportion"
stats_melt_info$Name <- factor(stats_melt_info$Name,levels = unique(stats_melt_info$Name))
stats_info$disease_type_update <- factor(stats_info$disease_type_update,levels=unique(stats_info$disease_type_update))
# samplenumber is used for annotation
samplenumber <- as.vector(table(stats_info$disease_type_update))
max = c()
ini =  0
for (i in samplenumber) {
  ini =  ini  + i
  max =  c(max,ini)
}
# chrssloc is used for annotation
chrssloc = as.vector(stats_info$Chrss)
max2 =c()
ini = 0
for (i in chrssloc) {
  ini =  ini  + 1
  if (i == 1) {
    max2 =  c(max2,ini)
  }
}

ggplot(stats_melt_info,
       aes(Name, Number, fill = Type),
       main = "par(las=2)") +
  geom_bar(stat = "identity",
           position = "fill",
           width = 1) +
  ggtitle("SV proportion_v7") +
  labs(y = "Percentage",  x = "") +
  scale_fill_wsj("colors6", "")  +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.ticks.x = element_blank()
  ) +
  guides(fill = guide_legend(title = NULL)) +
  annotate(
    "rect",
    xmin = -1,
    xmax = 797,
    ymin = -0.01,
    ymax = -0.05,
    fill = "cyan4"
  ) +
  annotate(
    "rect",
    xmin = max,
    xmax = max,
    ymin = -0.01,
    ymax = -0.05,
    colour = "white"
  ) +
  annotate(
    "rect",
    xmin = max2,
    xmax = max2,
    ymin = -0.06,
    ymax = -0.1,
    colour = "dodgerblue4"
  ) +
  expand_limits(y = c(-0.1, 1))

ggsave("SV proportion_v7.tiff",dpi = 1000,width = 20, height = 10, units = "cm")


# 2. plot "SV number across tumors"
# order by sample number SampleNumber and sum
order3 <- order(stats_info$SampleNumberUpdate,stats_info$SUM,decreasing = TRUE)
# do order
stats_info <- stats_info[order3,]
# put "Others" at the end
levels <- as.vector(unique(stats_info$disease_type_update))
levels <- levels[-which(levels=="Others")]
levels <- c(levels,"Others")

stats_info$disease_type_update <- factor(stats_info$disease_type_update,levels = levels)
stats_info$Name <- factor(stats_info$Name,levels = unique(stats_info$Name))

p <- ggplot(stats_info,aes(x=Name,y=SUM,color=disease_type_update))
pp <- p+geom_point(size=1)+facet_grid(. ~ disease_type_update,scales="free_x")  # to plot multi panels in a figure, scales="free_x":x axis has different scale, but y axis has fixed scale
ppp <- pp+theme_bw() + scale_y_log10()+scale_fill_manual("legend_title") +labs(col="Tumor Type")
pppp <- ppp+theme(axis.text.x=element_blank(),axis.text.y=element_text(size=10), # cancel x text
               axis.title.x=element_blank(),axis.title.y=element_text(size=10),
               axis.ticks.length=unit(0,"cm"), # cancel ticks
               panel.background = element_rect(fill = "gray94",
                                               colour = "gray94",
                                               size = 0.5, linetype = "solid"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                               colour = "gray94"), 
               panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                               colour = "gray94"),
               strip.background = element_blank(),
               strip.text.x = element_blank(),
               legend.title=element_text(size=10), 
               legend.text=element_text(size=9),
               legend.position="bottom",axis.line=element_line(colour="black"), # cancel legend, change x & y axis lines to black
               axis.line.x=element_line(size=2),axis.line.y=element_line(size=1),line=element_line(size=1.5))#line size of x is 2, line size of y is 1, scale line size is 1.5
# generate 16 colors
palette <- distinctColorPalette(16)
pppp + scale_color_manual(values = palette)

ggsave("SV number across tumors.tiff",dpi = 2000,width = 35, height = 20, units = "cm")


# 3. Plot "SV chromothripsis percentage"
# chrss_summary is a new data.frame, we save chromothripsis  mean value in it
chrss_summary <- data.frame(stats_info[,"Name"],stats_info[,"disease_type_update"],stats_info[,"Chrss"])
names(chrss_summary) <- c("Name","disease_type_update","Chrss")
# caculate mean value
chrssmean <-  tapply(chrss_summary$Chrss,chrss_summary$disease_type_update,mean)
chrssmean  <- as.data.frame(chrssmean)
chrssmean$Name  <- row.names(chrssmean)
# order chrssmean
chrssmean <- chrssmean[order(chrssmean$chrssmean,decreasing = TRUE),]
chrssmean$Name <- factor(chrssmean$Name, levels=factor(unique(chrssmean$Name))) 

ggplot(chrssmean,aes(x=Name,y=chrssmean,fill =  chrssmean))+geom_col()+theme(axis.text.x = element_text(angle = 90))+labs(y="Chromothripsis percentage",x="Tumor Type") +scale_fill_gradient(low = "blue", high = "red")

ggsave("SV chromothripsis percentage.tiff",dpi = 1000,width = 10, height = 20, units = "cm")
