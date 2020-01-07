## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----include=FALSE------------------------------------------------------------
require(taRifx)
require(data.table)
require(tidyr)
require(circlize)
require(GenomicRanges)
require(D3GB)

## ----message=FALSE------------------------------------------------------------
library(svcnvplus)
dim(nbl_segdat)
head(nbl_segdat)

dim(nbl_svdat)
head(nbl_svdat)

## -----------------------------------------------------------------------------
segdf <- validate.seg(nbl_segdat)
head(segdf)

## -----------------------------------------------------------------------------
svdf <- validate.sv(nbl_svdat)
head(svdf)

## ----message=FALSE, echo = FALSE----------------------------------------------
segdf_clean <- clean.cnv.artifact(segdf, verbose=FALSE,n.reps = 4)  # remove likely artifacts from segmentation data

## ----plot1,  fig.width=9, fig.height=4, fig.align='center', fig.cap = "Genome wide CNV frequencies", message=FALSE----
cnv_freq <- cnv.freq.plot(segdf)  # plot cnv frequencies

## ----message=FALSE------------------------------------------------------------
head(cnv_freq$freqsum)  # data.frame contains every genomic bin 


## ----message=FALSE------------------------------------------------------------
gene_cnv <- gene.cnv(segdf_clean,genome.v = "hg19",fill.gaps = TRUE,verbose=FALSE)

## -----------------------------------------------------------------------------
dim(gene_cnv$cnvmat)

## ----message=FALSE------------------------------------------------------------
pct_change <- pct.genome.changed(segdf)

## ----plot2,  fig.width=9, fig.height=5, fig.align='center', fig.cap = "Common breakpoints by sample",  message=FALSE----
common.breaks <- match.variant.breaks(segdf,svdf,maxgap=50000,low.cov = NULL,verbose=FALSE)

restab <- data.frame(common.breaks$restab)[order(common.breaks$restab$total.sv),]
m2 <- sprintf("%.1f",100*mean(restab$matched.sv/restab$total.sv))
barplot(rbind(restab$matched.sv,restab$total.sv - restab$matched.sv),
        border=NA,las=2,xlab="",horiz=FALSE,cex.main=.7,cex.names=.4, names=rownames(restab))
legend("top",paste("SV breaks matched by CNV breaks\n","Average = ",m2,"%",sep=""),bty='n')
grid(ny=NULL,nx=NA)

## ----plot3, fig.width=5, fig.height=5, fig.align='center', fig.cap = "SV versus CNV breakpoint burden",  message=FALSE----
sv_breaks_df  <- sv.breaks(svdf)  # define breakpoints from SV data
sv_burden <- table(sv_breaks_df$sample)
seg_breaks_df  <- seg.breaks(segdf,fc.pct = 0.2,verbose=FALSE)  # define breakpoints from seg data based on certain CNV change cutoff
seg_burden <- table(seg_breaks_df$sample)
common_samples <- intersect(names(sv_burden),names(seg_burden))
dat <- cbind(sv_burden[common_samples],seg_burden[common_samples])
plot(dat,xlab="SV burden",ylab="CNV breakpoint burden")
legend("topright",paste("cor=",cor(dat)[1,2], sep=""))

## ----message=FALSE------------------------------------------------------------
segdf <- validate.seg(segdat_lung_ccle)
svdf <- validate.sv(svdat_lung_ccle)

## ----message=FALSE------------------------------------------------------------
shatt_lung_cnv <- shattered.regions.cnv(segdf, fc.pct = 0.2, clean.brk = 4, window.size = 10,
                                        slide.size = 2,num.breaks = 8, num.sd = 5,  
                                        dist.iqm.cut = 150000,verbose=FALSE)

shatt_lung_cnv$regions.summary$A549_LUNG

## ----message=FALSE------------------------------------------------------------
shatt_lung <- shattered.regions(segdf, svdf, fc.pct = 0.2,  min.num.probes = 5, clean.brk = 8,
                                window.size = 10,slide.size = 2, num.seg.breaks = 6, 
                                num.seg.sd = 5, num.sv.breaks = 6, num.sv.sd = 5, 
                                num.common.breaks = 3, num.common.sd = 0, interleaved.cut = 0.33,
                                dist.iqm.cut = 100000,verbose=FALSE)
shatt_lung$regions.summary$SCLC21H_LUNG

## ----plot4, fig.width=5, fig.height=5, fig.align='center', fig.cap = "Circos plot representing c LUNG cancer cell line with chromothripsis",  message=FALSE----
circ.chromo.plot(shatt_lung,sample.id = "SCLC21H_LUNG")

## ----message=FALSE------------------------------------------------------------
fdr.test <- freq.p.test(shatt_lung_cnv$high.density.regions.hc, method="bonferroni", p.cut = 0.05)

## ----plot5,  fig.width=9, fig.height=4, fig.align='center', fig.cap = "Recurrently shattered regions map", message=FALSE----
recurrent.regions.plot(shatt_lung_cnv, fdr = fdr.test$freq.cut)

## ----message=FALSE------------------------------------------------------------
# obtain genomic bins within above the FDR cutoff
freq.matrix <- apply(shatt_lung_cnv$high.density.regions.hc,2,sum)
textRegions <- names(which(freq.matrix >= fdr.test$freq.cut))
hitRegions <- remove.factors((data.frame(do.call(rbind,strsplit(textRegions," ")))))
hitRegions[,2] <- as.numeric(hitRegions[,2])
hitRegions[,3] <- as.numeric(hitRegions[,3])
colnames(hitRegions) <- c("chr","start","end")
rownames(hitRegions) <-textRegions

# collapes contiguous bins into unique regions
bins2remove <- c()
for(i in 2:nrow(hitRegions)){ 
  if(hitRegions[i,"chr"] == hitRegions[i-1,"chr"] ){
    if(hitRegions[i,"start"] < (hitRegions[i-1,"end"])){
      hitRegions[i,"start"] <- hitRegions[i-1,"start"]
      bins2remove <- c(bins2remove,textRegions[i-1])
    }
  }
}
hitRegionsPost<- hitRegions[setdiff(rownames(hitRegions),bins2remove),]

require(GenomicRanges)
hitRegions_gr <- with(hitRegions, GRanges(chr, IRanges(start=start, end=end)))
hitRegionsPost_gr <- with(hitRegionsPost, GRanges(chr, IRanges(start=start, end=end)))
hits <-GenomicAlignments::findOverlaps(hitRegionsPost_gr,hitRegions_gr)

regList <- list()
for(i in unique(queryHits(hits))) regList[[paste(hitRegionsPost[i,],collapse=" ") ]] <- textRegions[subjectHits(hits)[which(queryHits(hits) == i)]]
# obtain the genomic bins with maximum number of samples
regListPeak <- lapply(regList, function(x) names(which(freq.matrix[x] == max(freq.matrix[x]))))
# collect samples with shattered region in the peaks 
regListPeakSamples <- lapply(regListPeak, function(x) names(which(apply(cbind(shatt_lung_cnv$high.density.regions.hc[,x]),1,sum) > 0)))

## ----message=FALSE------------------------------------------------------------

results_cnv <- cnv.gene.breaks(nbl_segdat, fc.pct = 0.2, genome.v="hg19",clean.brk = 8,upstr = 50000,verbose=FALSE)

## ----plot6,  fig.width=9, fig.height=4, fig.align='center', fig.cap = "Recurrently altered genes with overlapping CNV breakpoints", message=FALSE----
barplot(sort(unlist(lapply(results_cnv$geneSamples,length)),decreasing=T)[1:20],las=2)

## ----plot7,  fig.width=9, fig.height=4, fig.align='center', fig.cap = "Recurrently altered genes with overlapping CNV breakpoints", message=FALSE----
barplot(sort(unlist(lapply(results_cnv$upstreamSamples,length)),decreasing=T)[1:20],las=2)

