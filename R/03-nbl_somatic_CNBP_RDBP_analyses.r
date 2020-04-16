rm(list = ls(all.names = TRUE))
### look at the density of SV events in specific chromosomic regions
library(devtools)  # not sure if needed to load https
library(beeswarm)
library("BSgenome.Hsapiens.UCSC.hg19")
library(GenomicRanges)
library(GenomeInfoDb)
class(Hsapiens)
library(devtools)  # not sure if needed to load https
library(data.table)
library(tidyr)
require(taRifx)  # contains remove.factors

source('R/my_stat_functions.r')
source('R/01-nbl_somatic_SV_FUNCTIONS.R')


# load gene/exon annotations and generate Ranges objects
load("data/ucsc_hg19_refseq_genes_exon_df_Oct31_2018.rda",verbose=T)
features_tab_GR = with(genes_tab[,c("seqnames","start","end","strand")], GRanges(seqnames, IRanges(start=start, end=end),Rle(strand)))
mcols(features_tab_GR)$name2 <-rownames(genes_tab)


# load segmentation data for 914 samples from the NBL SNP array dataset
load("data/cnv_segmentation_SNP_081318.rda",verbose=T)
# load segmentation data for 135 CGI WGS samples + 884 tumors from ALL, AML, OS and WT
load("data/cnv_segmentation_CGI_081318.rda",verbose=T)

# load low coverage regions (subtelomeric and pericentromeric defined regions)
low_cov <- read.delim("data/CGI_low_coverage_regions.txt",as.is=T,header=F)
colnames(low_cov) <- c("chr","start","stop")
low_cov_GR = with(low_cov, GRanges(chr, IRanges(start=start, end=stop)))
overlap_gr <- GenomicAlignments::findOverlaps(features_tab_GR,low_cov_GR,ignore.strand=TRUE)
mapped_genes <- rownames(genes_tab[setdiff(1:nrow(genes_tab),queryHits(overlap_gr)),])

## obtain amplifications and deep deletions from segmentation datasets 
# reformat segmentation data
segment_cgi2<-segment_cgi
segment_snp2<-segment_snp
segment_4tum2 <- segment_4tum
segment_snp2[,"Chromosome"] <- paste("chr",segment_snp2[,"Chromosome"],sep="")
segment_cgi2[,"Chromosome"] <- paste("chr",segment_cgi2[,"Chromosome"],sep="")
segment_4tum2[,"Chromosome"] <- paste("chr",segment_4tum2[,"Chromosome"],sep="")
colnames(segment_4tum2) <- colnames(segment_cgi2)

# run analysis (FUNCTION=readDepthCopynum)
results_CN <- readDepthCopynum(segment_cgi2, genes_tab[mapped_genes,], ampl_cut = 2, ddel_cut = -1.9)
results_CN_snp<- readDepthCopynum(segment_snp2, genes_tab[mapped_genes,], ampl_cut = 1.2, ddel_cut = -1.2)
results_CN_4tum<- readDepthCopynum(segment_4tum2, genes_tab[mapped_genes,], ampl_cut = 1.2, ddel_cut = -1.2)

### Identification of recurrently altered altered genes using CNV breakpoint localization 
# a graphic description of this analysis is represented in Supplementary Figure S11B


upstr=100000 
dnstr=25000
promoter=1000
offset=100
copynumsize = 2000000
feature_tab<- genes_tab
exons_tab <- exons_tab
segdat <- segment_cgi
breaks_cgi <- readDepthBreaks(segdat, cutoff=0.304 ,segsize=10000,lowcov=low_cov)
breaks<- breaks_cgi[which(breaks_cgi$stop -breaks_cgi$start < 20000),]


breaks_redund_left <-  unite(breaks, newcol, c(chr,start), remove=FALSE,sep=":")$newcol
breaks_redund_right <-  unite(breaks, newcol, c(chr,stop), remove=FALSE,sep=":")$newcol
breaks_redund <- intersect(which(! breaks_redund_left %in% names(which(sort(table(breaks_redund_left)) > 3))),
	which(! breaks_redund_right %in% names(which(sort(table(breaks_redund_right)) > 3))))
breaks <- breaks[breaks_redund,]

rownames(breaks) <- unite(breaks, newcol, c(sample, chr,start,stop), remove=FALSE,sep=":")$newcol
results_BP <- enrichBP(breaks,feature_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100)

sort(unlist(lapply(results_BP$proximalSamples,length)),decreasing=T)[1:20]
sort(unlist(lapply(results_BP$breakSamples,length)),decreasing=T)[1:20]
sort(unlist(lapply(merge2lists(results_BP$breakSamples,results_BP$proximalSamples),length)),decreasing=T)[1:20]

upstr=100000 
dnstr=25000
promoter=1000
offset=100
copynumsize = 2000000
feature_tab<- genes_tab
exons_tab <- exons_tab
segdat <- segment_4tum2
breaks_cgi <- readDepthBreaks(segdat, cutoff=0.304 ,segsize=10000,lowcov=low_cov)
breaks<- breaks_cgi[which(breaks_cgi$stop -breaks_cgi$start < 20000),]

breaks_redund_left <-  unite(breaks, newcol, c(chr,start), remove=FALSE,sep=":")$newcol
breaks_redund_right <-  unite(breaks, newcol, c(chr,stop), remove=FALSE,sep=":")$newcol
breaks_redund <- intersect(which(! breaks_redund_left %in% names(which(sort(table(breaks_redund_left)) > 3))),
                           which(! breaks_redund_right %in% names(which(sort(table(breaks_redund_right)) > 3))))
breaks <- breaks[breaks_redund,]

rownames(breaks) <- unite(breaks, newcol, c(sample, chr,start,stop), remove=FALSE,sep=":")$newcol
results_BP_4tum  <- enrichBP(breaks,feature_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100)

sort(unlist(lapply(results_BP_4tum$proximalSamples,length)),decreasing=T)[1:20]
sort(unlist(lapply(results_BP_4tum$breakSamples,length)),decreasing=T)[1:20]
sort(unlist(lapply(merge2lists(results_BP_4tum$breakSamples,results_BP_4tum$proximalSamples),length)),decreasing=T)[1:20]


upstr=100000 
dnstr=25000
promoter=1000
offset=100
copynumsize = 2000000
feature_tab<- genes_tab
exons_tab <- exons_tab
segdat <- segment_snp
breaks_snp <- readDepthBreaks(segdat, cutoff=0.152 ,segsize=10000,lowcov=low_cov)
breaks<- breaks_snp[which(breaks_snp$stop -breaks_snp$start < 20000),]
breaks_redund_left <-  unite(breaks, newcol, c(chr,start), remove=FALSE,sep=":")$newcol
breaks_redund_right <-  unite(breaks, newcol, c(chr,stop), remove=FALSE,sep=":")$newcol
breaks_redund<- intersect(which(! breaks_redund_left %in% names(which(sort(table(breaks_redund_left)) > 3))),
	which(! breaks_redund_right %in% names(which(sort(table(breaks_redund_right)) > 3))))
breaks <- breaks[breaks_redund,]
rownames(breaks) <- unite(breaks, newcol, c(sample, chr,start,stop), remove=FALSE,sep=":")$newcol
results_BP_snp <- enrichBP(breaks,feature_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100)

sort(unlist(lapply(results_BP_snp$breakSamples,length)),decreasing=T)[1:20]
sort(unlist(lapply(results_BP_snp$proximalSamples,length)),decreasing=T)[1:20]

## save results for 
save(results_BP, results_BP_snp,results_BP_4tum, genes_tab, exons_tab, results_CN,results_CN_4tum, results_CN_snp, low_cov_GR,
	file="data/BP_analysis_Nov15_19.rda")


### break snp seg data into groupsfor GISTIC analysis ###
#########################################################

load("data/clinical_COG_20181129_plus_TARGET_20180331_plus_SCA_20180619_pheno.rda",verbose=T)

## segmentation for HR-NA
seg_noa <- segment_snp[which(segment_snp$Sample %in% intersect(hr_noa,segment_snp$Sample)),]
write.table(seg_noa,file="GISTIC_SNP_GROUPS/SNP_singleTumor_422_HR_NA.txt",sep="\t",quote=F,row.names=F)

## segmentation for MNA
seg_amp <- segment_snp[which(segment_snp$Sample %in% intersect(hr_amp,segment_snp$Sample)),]
write.table(seg_amp,file="GISTIC_SNP_GROUPS/SNP_singleTumor_422_HR_MNA.txt",sep="\t",quote=F,row.names=F)

## segmentation for LOWINT
seg_lint <- segment_snp[which(segment_snp$Sample %in% intersect(c(lowrisk,intrisk) ,segment_snp$Sample)),]
write.table(seg_lint,file="GISTIC_SNP_GROUPS/SNP_singleTumor_422_HR_LOWINT.txt",sep="\t",quote=F,row.names=F)



