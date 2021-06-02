#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

## read data and rename df header
dat <- read.table(args[1], sep = "\t")
names(dat) <- c('sample','mutation','project','caller','tumor_type')

## generate combined boxplot as well as plots for each caller
boxplot_all <- ggplot(
    dat,
    aes(reorder(tumor_type, mutation, FUN=median), 
    mutation, color=project
)) + geom_boxplot() + ylim(0, 200) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
facet_grid(rows=vars(caller))

# boxplot_lancet <- ggplot(
#     subset(dat, caller=="lancet"), 
#     aes(reorder(tumor_type, mutation, FUN=median), 
#     mutation, color=project
# )) + geom_boxplot() + ylim(0, 200) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# boxplot_mutect2 <- ggplot(
#     subset(dat, caller=="mutect2"), 
#     aes(reorder(tumor_type, mutation, FUN = median), 
#     mutation, color=project
# )) + geom_boxplot() + ylim(0, 200) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# boxplot_strelka2 <- ggplot(
#     subset(dat, caller=="strelka2"), 
#     aes(reorder(tumor_type, mutation, FUN = median), 
#     mutation, color=project
# )) + geom_boxplot() + ylim(0, 200) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


## output plots in png format
ggsave(
    paste(args[2], "boxplot-all.png", sep="/"),
    boxplot_all, device="png",
    width=10, height=5
)
# ggsave(
#     paste(args[2], "boxplot-lancet.png", sep="/"),
#     boxplot_lancet, device="png",
#     width=10, height=4
# )
# ggsave(
#     paste(args[2], "boxplot-mutect2.png", sep="/"),
#     boxplot_mutect2, device="png",
#     width=10, height=4
# )
# ggsave(
#     paste(args[2], "boxplot-strelka2.png", sep="/"),
#     boxplot_strelka2, device="png",
#     width=10, height=4
# )