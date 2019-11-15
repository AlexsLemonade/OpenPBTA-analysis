# this script is for generating chromothripsis and cnv figures
# we will use those figures to classify various chromothripsis mechanism
#
# Input files:
# 1.PBTA_chromothripsis_newlink_region.csv
# this file was generated from "02-shatterseek.R" script
# 2.independent-specimens.wgs.primary-plus.tsv
# this file is used to choose independent specimens
# 3. pbta-cnv-cnvkit.seg
# we will plot cnv 
#
# Output files:
# jpeg figures for each chromothripsis event


library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(cowplot)
library(ggforce)

`%!in%` = Negate(`%in%`)

## ===================== Load Shatterseek_link file ===================== 
list <- read.csv("D:/Project/Shatterseek/PBTA_v7/PBTA_chromothripsis_newlink_region.csv")
# read file
list=list[order(list$cluster_id,as.numeric(list$chr)),]
# unique cluster
listuni <-  list[!duplicated(list[,c('cluster_id')]),]

## ===================== Load  Independent specimen list =====================
independent_specimen_list <- read.table("D:/Data/PBTA/release-v7-20191031/independent-specimens.wgs.primary-plus.tsv",header = TRUE,sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load cnv file =====================
cnv_file <- read.table("D:/Data/PBTA/release-v7-20191031/pbta-cnv-cnvkit.seg/pbta-cnv-cnvkit.seg",header = TRUE,sep = "\t")
# choose independent specimens
cnv_analysis <-  cnv_file[cnv_file$ID %in% bioid & cnv_file$chrom != "chrY",]
# remove "chr", because shatterseek can't recognize it
cnv_analysis$chrom <- gsub("chr","",cnv_analysis$chrom)


for (j in 1:length(listuni$sample)){ 
  seq <- unlist(strsplit(as.character(listuni$link_chromosome[j]), split="_"))
  seq <- gsub("X","23",seq)
  seq <- gsub("Y","24",seq)
  seq <- gsub("M","25",seq)
  seq <- as.numeric(seq)
  sample <- as.character(listuni$sample[j])
  
  #  region  is specific  cluster in list
  region <- list[list$cluster_id == listuni$cluster_id[j],]
  region$chr <- gsub("X","23",region$chr)
  region$chr <- gsub("Y","24",region$chr)
  region$chr <- gsub("M","25",region$chr)
  
  # read cnv file
  cnv <- cnv_analysis[cnv_analysis$ID == sample,]
  # XY to 23 24
  cnv$chrom <- gsub("X","23",cnv$chrom)
  cnv$chrom <- gsub("Y","24",cnv$chrom)
  cnv$chrom <- gsub("M","25",cnv$chrom)
  # filter cnv in seq
  cnv <- cnv[cnv$chrom %in% seq,]
  
  # read sv file
  sv.data <- read.delim(paste("D:/Project/Signature/PBTA/V7/sv/",sample,".vcf",sep = ""))
  
  # XY to 23 24
  sv.data$chrom1 <- gsub("X","23",sv.data$chrom1)
  sv.data$chrom2 <- gsub("X","23",sv.data$chrom2)
  sv.data$chrom1 <- gsub("Y","24",sv.data$chrom1)
  sv.data$chrom2 <- gsub("Y","24",sv.data$chrom2)
  sv.data$chrom1 <- gsub("M","25",sv.data$chrom1)
  sv.data$chrom2 <- gsub("M","25",sv.data$chrom2)
  # filter sv in seq
  sv.data <- sv.data[sv.data$chrom1 %in% seq | sv.data$chrom2 %in% seq,]
  # sv position
  sv.data$position1 <- ifelse(sv.data$chrom1 %!in% seq,sv.data$pos2-1,sv.data$pos1)
  sv.data$position2 <- ifelse(sv.data$chrom2 %!in% seq,sv.data$pos1+1,sv.data$pos2)
  sv.data$chrom2 <- ifelse(sv.data$chrom2 %!in% seq,sv.data$chrom1,sv.data$chrom2)
  sv.data$chrom1 <- ifelse(sv.data$chrom1 %!in% seq,sv.data$chrom2,sv.data$chrom1)
  sv.data <- sv.data[!duplicated(sv.data[,c('chrom1','position1','chrom2','position2','SVtype')]),]
  # 23 24 to XY
  sv.data$chrom1 <- gsub("23","X",sv.data$chrom1)
  sv.data$chrom2 <- gsub("23","X",sv.data$chrom2)
  sv.data$chrom1 <- gsub("24","Y",sv.data$chrom1)
  sv.data$chrom2 <- gsub("24","Y",sv.data$chrom2)
  sv.data$chrom1 <- gsub("25","M",sv.data$chrom1)
  sv.data$chrom2 <- gsub("25","M",sv.data$chrom2)
  
  # Get the sv chromosomes...
  chr1 <- paste0("chr", sv.data$chrom1)
  chr2 <- paste0("chr", sv.data$chrom2)
  
  # ...positions
  pos1 <- sv.data$position1
  pos2 <- sv.data$position2
  
  # ... and the events type
  type <- sv.data$SVtype
  
  # genome length
  genome.length <- sum(seqlengths(Hsapiens)[seq])
  
  # chrs start x coordinate
  chrs_fake_starts <- vector("list", length(seq))
  chrs_fake_starts  <- setNames(chrs_fake_starts,  names(Hsapiens)[seq] )
  chrs_fake_starts[[names(Hsapiens)[seq[1]]]] <- 0
  length_sum <- 0
  if (length(seq)>1){
    for ( i in 2:length(seq) ) { 
      length_sum <- length_sum + as.numeric(seqlengths(Hsapiens)[[seq[i-1]]])+50000000
      chrs_fake_starts[[names(Hsapiens)[seq[i]]]] <- length_sum
    }
  }
  
  
  # assign 500k to cnv and region to make the segement visible
  if (length(seq) < 3){
    for (i in 1:length(seq)) {
      cnv[cnv$chrom == seq[i],]$loc.start <- cnv[cnv$chrom == seq[i],]$loc.start + chrs_fake_starts[[i]]
      cnv[cnv$chrom == seq[i],]$loc.end <- cnv[cnv$chrom == seq[i],]$loc.end + chrs_fake_starts[[i]] + 500000
      region[region$chr==seq[i],]$start <- region[region$chr==seq[i],]$start + chrs_fake_starts[[i]]
      region[region$chr==seq[i],]$end <- region[region$chr==seq[i],]$end + chrs_fake_starts[[i]] + 500000
    }
  } else {
    for (i in 1:length(seq)) {
      cnv[cnv$chrom==seq[i],]$loc.start <- cnv[cnv$chrom==seq[i],]$loc.start + chrs_fake_starts[[i]]
      cnv[cnv$chrom==seq[i],]$loc.end <- cnv[cnv$chrom==seq[i],]$loc.end + chrs_fake_starts[[i]]+500000
      region[region$chr==seq[i],]$start=region[region$chr==seq[i],]$start+chrs_fake_starts[[i]]
      region[region$chr==seq[i],]$end=region[region$chr==seq[i],]$end+chrs_fake_starts[[i]]+500000
    }
  }
  
  
  #set lable coordinate and boundary coordinate
  chrs_fake_label.pos <- vector("list", length(seq))
  chrs_fake_label.pos <- setNames(chrs_fake_label.pos,  names(Hsapiens)[seq] )
  chrs_boundary.pos <- vector("list", length(seq))
  chrs_boundary.pos <- setNames(chrs_fake_label.pos,  names(Hsapiens)[seq] )
  for ( i in 1:length(chrs_fake_starts) ) {
    chrs_fake_label.pos[[names(Hsapiens)[seq[i]]]] <- seqlengths(Hsapiens)[[seq[i]]]/2 + chrs_fake_starts[[names(Hsapiens)[[seq[i]]]]]
    chrs_boundary.pos[[names(Hsapiens)[seq[i]]]] <- seqlengths(Hsapiens)[[seq[i]]] + chrs_fake_starts[[names(Hsapiens)[[seq[i]]]]]
  }
  
  # set fake coordinate
  pos1_fake <- vector("list", nrow(sv.data))
  pos2_fake <- vector("list", nrow(sv.data))
  for ( i in 1:nrow(sv.data) ) {
    pos1_fake[[i]] <- chrs_fake_starts[[chr1[i]]] + pos1[i]  #if chr not a chrss_chr,  value will be interger
    pos2_fake[[i]] <- chrs_fake_starts[[chr2[i]]] + pos2[i]  #if chr not a chrss_chr,  value will be interger
  }
  
  #set bezier  coordinate
  beziers.height <- runif(nrow(sv.data), 0.5, 1)
  beziers.mid <- (unlist(pos2_fake) + unlist(pos1_fake))/2  
  
  # set baziers dataframe
  beziers <- data.frame(
    x <- c(rbind( unlist(pos1_fake), beziers.mid, unlist(pos2_fake) )),
    y <- c(rbind( 0.2, beziers.height, 0.2 ) ),
    group <- rep( paste( chr1, pos1, chr2, pos2,type, sep="_" ), each=3),
    svtype <- rep( type, each=3)
  )
  
  # give tra beziers a higher y
  beziers_intra <- beziers[beziers$svtype!='TRA',]
  beziers_inter <- beziers[beziers$svtype=='TRA',]
  beziers_inter$y <- beziers_inter$y+0.5
  
  # 23 24 to XY
  cnv <- na.omit(cnv)
  cnv$chrom <- gsub("23","X",cnv$chrom)
  cnv$chrom <- gsub("24","Y",cnv$chrom)
  cnv$chrom <- gsub("25","M",cnv$chrom)
  
  # set boundary and boundary_even
  boundary  <- as.data.frame(c(unlist(chrs_boundary.pos),unlist(chrs_fake_starts)))
  colnames(boundary) <- "position"
  boundary <- as.data.frame(sort(boundary$position,index.return=FALSE))
  colnames(boundary) <- "position"
  boundary_even <-  as.data.frame(boundary[seq(2,nrow(boundary),2),])
  colnames(boundary_even) <- "position"
  ymax <- max(c(cnv$copy.num,2))
  
  #set title coordinate
  title_position <- tail(unlist(chrs_boundary.pos),n=1)/2
  
  # set plot range
  plot_range <- ifelse(length(seq)>1,(tail(unlist(chrs_boundary.pos),n=1)+50000000),(tail(unlist(chrs_boundary.pos),n=1)+5000000))
  
  #plot cnv
  if (length(seq) > 1) {
    if (nrow(cnv) > 0 ){
      p2  <-  ggplot(cnv) + geom_segment(
        aes(
          x = (loc.start),
          y = copy.num,
          xend = (loc.end),
          yend = copy.num
        ),
        size = 5,
        colour = "black"
      )  + geom_segment(
        data = boundary,
        aes(
          x = position ,
          xend = position,
          y = 0,
          yend = as.numeric(ymax) # boundary
        ),
        linetype = 2,
        colour = 'black',
        size = 0.5
      ) + scale_x_continuous(
        expand = c(0, 0),
        breaks = NULL,
        limits = c(0, plot_range)
      ) + scale_y_continuous(
        trans = "sqrt",
        breaks = function(x)
          unique(c(floor(pretty(seq(
            0, (max(x) + 1) * 1.1
          ))),0,1,2,3))
      ) + ylab("CNV") + xlab("") + theme_bw() + theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.grid.major=element_line(colour='grey'),
        panel.grid.minor=element_blank()
      ) } else {
        p2  <-  ggplot(cnv) + geom_segment(
          aes(
            x = (loc.start),
            y = copy.num,
            xend = (loc.end),
            yend = copy.num
          ),
          size = 5,
          colour = "black"
        )  + geom_segment(
          data = boundary,
          aes(
            x = position ,
            xend = position,
            y = 0,
            yend = as.numeric(ymax) # boundary
          ),
          linetype = 2,
          colour = 'black',
          size = 0.5
        ) + scale_x_continuous(
          expand = c(0, 0),
          breaks = NULL,
          limits = c(0, plot_range)
        ) + scale_y_continuous(
          trans = "sqrt",
          breaks = function(x)
            unique(c(floor(pretty(seq(
              0, (max(x) + 1) * 1.1
            ))),0,1,2,3))
        ) + ylab("CNV") + xlab("") + theme_bw() + theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12),
          panel.grid.major=element_line(colour='grey'),
          panel.grid.minor=element_blank()
        )
      }
  } else {
    if (nrow(cnv) > 0) {
      p2  <-  ggplot(cnv) + geom_segment(
        aes(
          x = (loc.start),
          y = copy.num,
          xend = (loc.end),
          yend = copy.num
        ),
        size = 5,
        colour = "black"
      ) + geom_segment(
        data = boundary,
        aes(
          x = position ,
          xend = position,
          y = 0,
          yend = as.numeric(ymax) # boundary
        ),
        linetype = 2,
        colour = 'black',
        size = 0.5
      ) + scale_x_continuous(
        expand = c(0, 0),
        breaks = NULL,
        limits = c(0, plot_range)
      ) + scale_y_continuous(
        trans = "sqrt",
        breaks = function(x)
          unique(c(floor(pretty(seq(
            0, (max(x) + 1) * 1.1
          ))),0,1,2,3))
      ) + ylab("CNV") + xlab("") + theme_bw() + theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.grid.major=element_line(colour='grey'),
        panel.grid.minor=element_blank()
      )} else {
        p2  <-  ggplot(cnv) + geom_segment(
          aes(
            x = (loc.start),
            y = copy.num,
            xend = (loc.end),
            yend = copy.num
          ),
          size = 5,
          colour = "black"
        ) + geom_segment(
          data = boundary,
          aes(
            x = position ,
            xend = position,
            y = 0,
            yend = as.numeric(ymax) # boundary
          ),
          linetype = 2,
          colour = 'black',
          size = 0.5
        ) + scale_x_continuous(
          expand = c(0, 0),
          breaks = NULL,
          limits = c(0, plot_range)
        ) + scale_y_continuous(
          trans = "sqrt",
          breaks = function(x)
            unique(c(floor(pretty(seq(
              0, (max(x) + 1) * 1.1
            ))),0,1,2,3))
        ) + ylab("CNV") + xlab("") + theme_bw() + theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12),
          panel.grid.major=element_line(colour='grey'),
          panel.grid.minor=element_blank()
        )
      }
  }
  
  # set label size  and title size
  label_size <- 7-log2(length(seq))
  title_size <-  10-log2(length(seq))
  
  # 23 24 to XY
  region$chr <-  gsub("23","X",region$chr)
  region$chr <- gsub("24","Y",region$chr)
  region$chr <- gsub("25","M",region$chr)
  
  # give region
  region$fake_start <- as.numeric(chrs_fake_starts)
  region$fake_end <- as.numeric(chrs_boundary.pos)
  
  # define  link  region and chrss region
  link_region <- region[region$chrss_status=="link",]
  chrss_region <- region[region$chrss_status=="chrss",]
  
  # plot sv
  if (length(link_region$chr)==0) {
    print("1")
    print(listuni$cluster_id[j])
    p1  <-  ggplot() + geom_bezier(
      aes(
        x = x,
        y = y,
        group = group,
        color = svtype
      ),
      data = beziers_inter,
      show.legend = TRUE,
      size = 0.5
    ) + geom_bezier(
      aes(
        x = x,
        y = y,
        group = group,
        color = svtype
      ),
      data = beziers_intra,
      show.legend = TRUE,
      size = 0.5
    ) + scale_x_continuous(expand = c(0, 0), limits = c(0, plot_range)) + xlab("") +
      # Remove default axes labels and grey backgroud
      theme(
        axis.title.x = element_text(size = title_size),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # ...and the grey backgroud
        panel.background = element_rect(fill = NA),
        # ...change the legend parameters
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_blank(),
        legend.position = c(0.15, 0.8),
        legend.direction = "horizontal"
      ) +
      # Set the axes limits
      scale_y_continuous(limits = c(-0.2, 1.5)) +
      # Add chromosomes boundaries
      geom_segment(
        aes(
          x = c(unlist(chrs_boundary.pos), unlist(chrs_fake_starts)) ,
          xend = c(unlist(chrs_boundary.pos), unlist(chrs_fake_starts)),
          y = 0,
          yend = 0.2
        ),
        linetype = 2,
        colour = 'black',
        size = 0.5
      ) + geom_segment(
        data = chrss_region,
        aes(
          x = start,
          y = 0,
          xend = end,
          yend = 0
        ),
        color = "red",
        size = 5
      ) +
      # Add chromosomes labels
      annotate(
        geom = 'text',
        label = names(chrs_fake_label.pos),
        x = unlist(chrs_fake_label.pos),
        y = 0.1,
        size = label_size
      ) + scale_color_manual(
        "SV type",
        values = c(
          "DEL" = "#428BCA",
          "DUP" = "#d9534f",
          "h2hINV" = "#FFC425",
          "t2tINV" = "#5cb85c",
          "TRA" = "#aa96dc"
        )
      ) + annotate(
        geom = 'text',
        label = listuni$cluster_id[j],
        x = title_position,
        y = 1.45,
        size = title_size
      )
    ##############################################################################################
  } else if (length(chrss_region$chr) == 0) {
    print("2")
    print(listuni$cluster_id[j])
    p1 <- ggplot() + geom_bezier(
      aes(
        x = x,
        y = y,
        group = group,
        color = svtype
      ),
      data = beziers_inter,
      show.legend = TRUE,
      size = 0.5
    ) + geom_bezier(
      aes(
        x = x,
        y = y,
        group = group,
        color = svtype
      ),
      data = beziers_intra,
      show.legend = TRUE,
      size = 0.5
    ) + scale_x_continuous(expand = c(0, 0), limits = c(0, plot_range)) + xlab("") +
      # Remove default axes labels and grey backgroud
      theme(
        axis.title.x = element_text(size = title_size),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # ...and the grey backgroud
        panel.background = element_rect(fill = NA),
        # ...change the legend parameters
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_blank(),
        legend.position = c(0.15, 0.8),
        legend.direction = "horizontal"
      ) +
      # Set the axes limits
      scale_y_continuous(limits = c(-0.2, 1.5)) +
      # Add chromosomes boundaries
      geom_segment(
        aes(
          x = c(unlist(chrs_boundary.pos), unlist(chrs_fake_starts)) ,
          xend = c(unlist(chrs_boundary.pos), unlist(chrs_fake_starts)),
          y = 0,
          yend = 0.2
        ),
        linetype = 2,
        colour = 'black',
        size = 0.5
      ) + geom_segment(
        data = link_region,
        aes(
          x = start,
          y = 0,
          xend = end,
          yend = 0
        ),
        color = "blue",
        size = 5
      ) +
      # Add chromosomes labels
      annotate(
        geom = 'text',
        label = names(chrs_fake_label.pos),
        x = unlist(chrs_fake_label.pos),
        y = 0.1,
        size = label_size
      ) + scale_color_manual(
        "SV type",
        values = c(
          "DEL" = "#428BCA",
          "DUP" = "#d9534f",
          "h2hINV" = "#FFC425",
          "t2tINV" = "#5cb85c",
          "TRA" = "#aa96dc"
        )
      ) + annotate(
        geom = 'text',
        label = listuni$cluster_id[j],
        x = title_position,
        y = 1.45,
        size = title_size
      )
    ################################################################################################
  } else {
    print("3")
    print(listuni$cluster_id[j])
    p1 = ggplot() + geom_bezier(
      aes(
        x = x,
        y = y,
        group = group,
        color = svtype
      ),
      data = beziers_inter,
      show.legend = TRUE,
      size = 0.7
    ) + geom_bezier(
      aes(
        x = x,
        y = y,
        group = group,
        color = svtype
      ),
      data = beziers_intra,
      show.legend = TRUE,
      size = 0.7
    ) + scale_x_continuous(expand = c(0, 0), limits = c(0, plot_range)) + xlab("") +
      # Remove default axes labels and grey backgroud
      theme(
        axis.title.x = element_text(size = title_size),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # ...and the grey backgroud
        panel.background = element_rect(fill = NA),
        # ...change the legend parameters
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_blank(),
        legend.position = c(0.15, 0.8),
        legend.direction = "horizontal"
      ) +
      # Set the axes limits
      scale_y_continuous(limits = c(-0.2, 1.5)) +
      # Add chromosomes boundaries
      geom_segment(
        aes(
          x = c(unlist(chrs_boundary.pos), unlist(chrs_fake_starts)) ,
          xend = c(unlist(chrs_boundary.pos), unlist(chrs_fake_starts)),
          y = 0,
          yend = 0.2
        ),
        linetype = 2,
        colour = 'black',
        size = 0.5
      ) + geom_segment(
        data = link_region,
        aes(
          x = start,
          y = 0,
          xend = end,
          yend = 0
        ),
        color = "blue",
        size = 5
      ) + geom_segment(
        data = chrss_region,
        aes(
          x = start,
          y = 0,
          xend = end,
          yend = 0
        ),
        color = "red",
        size = 5
      ) +
      # Add chromosomes labels
      annotate(geom = 'text', label = names(chrs_fake_label.pos), x = unlist(chrs_fake_label.pos), y = 0.1, size = label_size)+scale_color_manual("SV type",values = c("DEL"="#428BCA","DUP"="#d9534f","h2hINV"="#FFC425","t2tINV"="#5cb85c","TRA"="#aa96dc"))+annotate(geom = 'text', label = listuni$cluster_id[j], x = title_position, y = 1.45, size = title_size)
  }
  
  # combine p1 and p2
  p3 <- plot_grid(p1,p2,ncol=1, align="v",rel_heights =c(2,1))
  
  # output figure
  jpegname <- paste0(listuni$cluster_id[j],"_cnvkit",".jpeg")
  jpeg(jpegname,width=10000,height=6000,res = 600)
  print(p3)
  dev.off()
}
