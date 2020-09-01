library('ggpubr')
stringsAsFactors=FALSE
library(stringr)
library(gridBase)
library(gridGraphics)
library(optparse)
library(forcats) ### for fct_reorder()

###################################################################################################

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 

A = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAPolya_counts.txt")   ### Variable representing Telomerase activity-prediction for Polya_counts
				 
B = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_counts.txt")  ### Variable representing Telomerase activity-prediction for Stranded_counts
				 
C = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAPolya_FPKM.txt")      ### Variable representing Telomerase activity-prediction for PolyA_FPKM data
				 
D = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_FPKM.txt")   ### Variable representing Telomerase activity-prediction for Stranded_FPKM data


palette_dir <- file.path(root_dir, "figures", "palettes")
histology_palette <- read_tsv(file.path(palette_dir,
                                        "histology_color_palette.tsv"))				 
colnames(histology_palette)[1]='short_histology'

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")
		 
telomerase_png <- file.path(output_dir, "Telomerase_Activities.png")

########################################## Figure A and B #########################################################

				 
PTBA_GE_PolyA_TMScores1 = read.table(A,sep='\t',head=T)
colnames(PTBA_GE_PolyA_TMScores1)[colnames(PTBA_GE_PolyA_TMScores1)=="NormEXTENDScores"]='NormEXTENDScores_PolyACounts'

PTBA_GE_Standard_TMScores1 = read.table(B,sep='\t',head=T)
colnames(PTBA_GE_Standard_TMScores1)[colnames(PTBA_GE_Standard_TMScores1)=="NormEXTENDScores"]='NormEXTENDScores_StrandedCounts'

PTBA_GE_PolyA_TMScores2 = read.table(C,sep='\t',head=T)
colnames(PTBA_GE_PolyA_TMScores2)[colnames(PTBA_GE_PolyA_TMScores2)=="NormEXTENDScores"]='NormEXTENDScores_PolyA_FPKM'

PTBA_GE_Standard_TMScores2 = read.table(D,sep='\t',head=T)
colnames(PTBA_GE_Standard_TMScores2)[colnames(PTBA_GE_Standard_TMScores2)=="NormEXTENDScores"]='NormEXTENDScores_Stranded_FPKM'

PTBA_GE_PolyA_TMScores = merge(PTBA_GE_PolyA_TMScores1,PTBA_GE_PolyA_TMScores2,by='SampleID')
PTBA_GE_Standard_TMScores = merge(PTBA_GE_Standard_TMScores1,PTBA_GE_Standard_TMScores2,by='SampleID')

########################################## Figure C and D #########################################################

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 

A = file.path(root_dir, "data", "pbta-histologies.tsv")   ### Variable representing clinical data 
				 
B = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_FPKM.txt")   ### Variable representing Telomerase activity-prediction for Stranded_FPKM data



PTBA_Histology = read.table(A,sep='\t',head=T)    ## Reading the clinical data

colnames(PTBA_Histology)[colnames(PTBA_Histology)=="Kids_First_Biospecimen_ID"]='SampleID'   ## Renaming "Kids_First_Biospecimen_ID" as SampleID for comparison purpose

TMScores = read.table(B,sep='\t',head=T)  ## Reading Stranded FPKM telomerase scores

PTBA_GE_Standard_Histology = merge(PTBA_Histology,TMScores,by='SampleID')   ### Merging Clinical data with the Telomerase scores

Stranded_Histology = PTBA_GE_Standard_Histology
Stranded_Histology = Stranded_Histology[-which(Stranded_Histology$short_histology == "Other"),]   ### Removing the tumors with catagory labelled as "Others"
Frequency = data.frame(table(Stranded_Histology$short_histology))  ### Counting the number of cases for all histologies to avoid less number further
colnames(Frequency)=c('Variables','Freq')
Frequency = Frequency[which(Frequency$Freq == 1),]
Stranded_Histology = Stranded_Histology[-which(Stranded_Histology$short_histology %in% Frequency$Variables),]     ### Removing the tumors with only one case in histologies 

Shorthist = merge(histology_palette,Stranded_Histology,by='short_histology')

########################################## Figure E #########################################################

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 

A = file.path(root_dir, "data", "pbta-histologies.tsv")   ### Variable representing clinical data 
				 
B = file.path(root_dir, "analyses", "telomerase-activity-prediction",
                 "results", "TelomeraseScores_PTBAStranded_FPKM.txt")   ### Variable representing Telomerase activity-prediction for Stranded_FPKM data


PTBA_Histology = read.table(A,sep='\t',head=T)    ## Reading the clinical data

colnames(PTBA_Histology)[colnames(PTBA_Histology)=="Kids_First_Biospecimen_ID"]='SampleID'   ## Renaming "Kids_First_Biospecimen_ID" as SampleID for comparison purpose

PTBA_GE_Standard_TMScores = read.table(B,sep='\t',head=T)  ## Reading Stranded FPKM telomerase scores

PTBA_GE_Standard_Histology = merge(PTBA_Histology,PTBA_GE_Standard_TMScores,by='SampleID')   ### Merging Clinical data with the Telomerase scores

PTBA_GE_Standard_Histology = PTBA_GE_Standard_Histology[which(PTBA_GE_Standard_Histology$short_histology == "Medulloblastoma"),]   ### Select tumors with catagory labelled as "Medulloblastoma"

stat.test <- data.frame(compare_means(
  NormEXTENDScores ~ molecular_subtype, data = PTBA_GE_Standard_Histology,
  method = "t.test"
)) 

combinations = nrow(stat.test)

statistics = stat.test%>%
  mutate(y.position = seq(1,by=0.04,length.out=combinations))




############# Saving PNG format 

png(telomerase_png,width = 6, height = 6, units = "in", res = 1200)

theme_set(theme_classic() +
          theme(plot.title = element_text(size=10, face="bold"),axis.text.x=element_text(angle=25,size=6,vjust=1,hjust=1),axis.text.y=element_text(size=7), axis.title.x = element_text(size=0), axis.title.y = element_text(size=8),
          legend.position = "none",
          legend.key.size= unit(0.3,"cm"),
          legend.key.width = unit(0.3,"cm"),
          legend.title = element_text(size=7),
          legend.text =element_text(size=6)
        )
)


P1 = ggscatter(PTBA_GE_PolyA_TMScores, x = "NormEXTENDScores_PolyACounts", y = "NormEXTENDScores_PolyA_FPKM",color = "red",size = 0.2,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)

		   
P2 = ggscatter(PTBA_GE_Standard_TMScores, x = "NormEXTENDScores_StrandedCounts", y = "NormEXTENDScores_Stranded_FPKM",color = "red",size = 0.2,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)



P3 = ggplot(Shorthist, aes(x=fct_reorder(short_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size=0.2,notch=FALSE,outlier.size = 0,outlier.shape=NA,aes(color=hex_codes,fill=hex_codes),alpha=0.4)+ geom_jitter(shape=16, cex=0.1,aes(color=hex_codes))



P4 = ggplot(Stranded_Histology, aes(x=fct_reorder(broad_histology,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size=0.2,notch=FALSE,outlier.size = 0,outlier.shape=NA,color="maroon",fill="pink",alpha=0.4)+ geom_jitter(shape=16, cex=0.1,color="maroon")


P5 = ggplot(PTBA_GE_Standard_Histology, aes(x=fct_reorder(molecular_subtype,NormEXTENDScores,.desc =TRUE),y=NormEXTENDScores))+geom_boxplot(size=0.2,notch=FALSE,outlier.size = 0,outlier.shape=NA,color="maroon",fill="pink",alpha=0.4)+ geom_jitter(shape=16, width = 0.1,size=0.2,color="maroon")+stat_pvalue_manual(
    data = statistics, label = "p.adj",size=1.7,
    xmin = "group1", xmax = "group2",tip.length = 0.003,
    y.position = "y.position"
    )

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 9, ncol =3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =6,font.ytickslab =6,font.x = 6,font.y=6,xlab="PolyA Counts EXTEND Scores",ylab="PolyA FPKM EXTEND Scores",title="A",font.title=7),vp = define_region(row = 1:3, col = 1))
print(ggpar(P2,font.xtickslab =6,font.ytickslab =6,font.x = 6,font.y=6,xlab="Stranded Counts EXTEND Scores ",ylab="Stranded FPKM EXTEND Scores",title="B",font.title=7),vp = define_region(row = 1:3, col = 2))
print(ggpar(P5,font.xtickslab =c(5,"black"),font.ytickslab =6,font.x = 6,font.y=6,font.legend=6,xlab="Medulloblastoma Subgroups",ylab="EXTEND Scores",title="E",font.title=7),vp = define_region(row = 1:3, col = 3))
print(ggpar(P3,font.xtickslab =c(5,"black"),font.ytickslab =6,font.x = 6,font.y=6,ylab="EXTEND Scores",xlab = "Tumor Short Histology",title="C",font.title=7),vp = define_region(row = 4:6, col = 1:3))
print(ggpar(P4,font.xtickslab =c(5,"black"),font.ytickslab =6,font.x = 6,font.y=6,ylab="EXTEND Scores",xlab = "Tumor Broad Histology",title="D",font.title=7),vp = define_region(row = 7:9, col = 1:3))

dev.off()

