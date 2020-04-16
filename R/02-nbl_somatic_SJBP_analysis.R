require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(devtools)
require(GenomicFeatures)
library(data.table)
require(taRifx)  # contains remove.factors

setwd("~/Box Sync/git/NB_structural_variants/")

source('R/my_stat_functions.r')
source('R/01-nbl_somatic_SV_FUNCTIONS.R')

## load all collapsed data by type of SVs
#setwd("~/Box Sync/My_CHOP/SV_paper_V2/rdata/")


# database of genomic variants
DGV <- read.delim("data/GRCh37_hg19_variants_2016-05-15-small.txt",as.is=TRUE)

# table for cosmic genes
#cosmic_genes <- read.delim(paste(init,"/diskin_lab/Gonzalo/SV_nbl/Census_allSun Mar 12 20-20-26 2017.csv",sep=""),sep=",")


# Load all CGI SV calls
load("data/allSomatic_cgi_SV_data_V2.rda",verbose=TRUE)
save(junctions_all,file="data/allSomatic_cgi_SV_data_V2.rda",compress=TRUE)

junall <- junctions_all[which(!duplicated(substr(names(junctions_all),0,16)))]
somatic_junctions_all_df <-  as.data.frame(rbindlist(junall))
rownames(somatic_junctions_all_df) <- unlist(sapply(names(junall),function(i) paste(i,rownames(junall[[i]]),sep=".") ))

# Obtain genomic features (genes and exons) derived from UCSC known genes (no need to run; just load the result file below)
#########################################
# refseq <- read.delim("~/Box Sync/My_CHOP/reference_files/ucsc_refseq_0ct31_2018.txt",as.is=T)
# chromosomes <- paste("chr",c(1:22,"X"),sep="")
# refseq<-refseq[which(refseq$chrom %in% chromosomes),]
# refseq<-refseq[which(!duplicated(refseq$name)),]

# refseq_data <- list()
# refseq_gene <- list()
# for(gene_id in unique(refseq$name2)){
#	refseq_g <- refseq[which(refseq$name2 == gene_id),]
#	#strand <- names(sort(table(refseq_g$strand),decreasing=T)[1])
#	strand <- refseq_g$strand[1]
#	refseq_g<-refseq_g[which(refseq_g$strand == strand),]
#	start <- min(refseq_g$txStart)
#	end <- max(refseq_g$txEnd)
#	seqnames <- unique(refseq_g$chrom)
#	width <- end-start
#	refseq_gene[[gene_id]] <- c(paste(seqnames,collapse=":"),start,end,width,paste(strand,collapse=":"),gene_id)
#	if(length(strand) > 1) message(paste(gene_id,strand))
#	if(length(seqnames) > 1) message(paste(gene_id,seqnames))
#	exStart <- as.numeric(unlist(strsplit(refseq_g$exonStarts,",")))
#	exEnd <- as.numeric(unlist(strsplit(refseq_g$exonEnds,",")))
#	exonPos <- do.call(rbind,strsplit(unique(apply(cbind(exStart[order((exEnd+exStart)/2)],exEnd[order((exEnd+exStart)/2)]),1,paste,collapse="_")),"_"))
#	refseq_data[[gene_id]] <- cbind(rep(gene_id,nrow(exonPos)),rep(seqnames,nrow(exonPos)),exonPos)
#	}
# genes_tab <- data.frame(do.call(rbind,refseq_gene))
# colnames(genes_tab) <- c("seqnames","start","end","width","strand","gene_id")
# genes_tab[,"start"] <-  as.numeric(as.character(genes_tab[,"start"]))
# genes_tab[,"end"] <-  as.numeric(as.character(genes_tab[,"end"]))
# genes_tab[,"width"] <-  as.numeric(as.character(genes_tab[,"width"]))
# genes_tab[,"seqnames"] <-  as.character(genes_tab[,"seqnames"])
# genes_tab[,"strand"] <-  as.character(genes_tab[,"strand"])
# genes_tab[,"gene_id"] <-  as.character(genes_tab[,"gene_id"])


# remove_genes <- unique(c(rownames(genes_tab)[grep(":",genes_tab$seqnames)],grep("-AS",rownames(genes_tab),value=T),grep("^MIR",rownames(genes_tab),value=T)))
# genes_tab <- genes_tab[setdiff(rownames(genes_tab),remove_genes),]
# exons_tab <- data.frame(do.call(rbind,refseq_data[rownames(genes_tab)]))
# colnames(exons_tab) <- c("gene_id","seqnames","start","end")
# exons_tab[,"start"] <-  as.numeric(as.character(exons_tab[,"start"]))
# exons_tab[,"end"] <-  as.numeric(as.character(exons_tab[,"end"]))
# exons_tab[,"seqnames"] <-  as.character(exons_tab[,"seqnames"])
# exons_tab[,"gene_id"] <-  as.character(exons_tab[,"gene_id"])
# width <- exons_tab$end-exons_tab$start
# exons_tab<-data.frame(exons_tab,width)

# save(exons_tab,genes_tab,file="data/ucsc_hg19_refseq_genes_exon_df_Oct31_2018.rda")

### load the data.franes generated with the code commented above
load("data/ucsc_hg19_refseq_genes_exon_df_Oct31_2018.rda",verbose=T)


### Filterin artifacts from CGI Junction and Event data
#########################################

# add patient id column to SV data.frame in lowC junction
TARGET.USI<-as.character(substr(rownames(somatic_junctions_all_df),0,16))
somatic_junctions_all_f0 <- cbind(TARGET.USI,somatic_junctions_all_df)
somatic_junctions_all_f0 <- somatic_junctions_all_f0[which(!somatic_junctions_all_f0$Type == "artifact"),]
somatic_junctions_all_f0 <- somatic_junctions_all_f0[which(somatic_junctions_all_f0$FrequencyInBaselineGenomeSet == "0" ),]
savesvs <- c(intersect(which(somatic_junctions_all_f0$Type %in% c("complex","interchromosomal","probable-inversion") ),which(somatic_junctions_all_f0$DiscordantMatePairAlignments >4)),which(!somatic_junctions_all_f0$Type %in% c("complex","interchromosomal","probable-inversion") ))
somatic_junctions_all_f0 <- somatic_junctions_all_f0[savesvs,]

##################################################################################
# FILTERING COMMON VARIANTS FROM THE SOMATIC DATASET
# create genomicRanges objects for each type of junction (deletion,inversion,etc...) in order to find overlaps with the database of genomic variants

intrachr_all <- which(unlist(lapply(apply(cbind(as.character(somatic_junctions_all_f0$LeftChr),as.character(somatic_junctions_all_f0$RightChr)),1,unique),length)) == 1)
all_complex <- intersect(which(somatic_junctions_all_f0$Type == "complex"),intrachr_all)
all_loss <- intersect(which(somatic_junctions_all_f0$Type == "deletion"),intrachr_all)
all_gain <- intersect(which(somatic_junctions_all_f0$Type == "tandem-duplication"),intrachr_all)
all_inv <- intersect(which(somatic_junctions_all_f0$Type %in% c("inversion","probable-inversion")),intrachr_all)
events_all2 <- data.frame(as.character(somatic_junctions_all_f0$LeftChr),somatic_junctions_all_f0$LeftPosition,somatic_junctions_all_f0$RightPosition + somatic_junctions_all_f0$RightLength)
rownames(events_all2) <- rownames(somatic_junctions_all_f0)
colnames(events_all2) <- c("chr","start","end")
events_all_complex <- events_all2[all_complex,]
events_all_loss <-  events_all2[all_loss,]
events_all_gain <-  events_all2[all_gain,]
events_all_inv <-  events_all2[all_inv,]
events_all_complex_GR =  with(events_all_complex, GRanges(chr, IRanges(start=start, end=end))) 
events_all_loss_GR =  with(events_all_loss, GRanges(chr, IRanges(start=start, end=end))) 
events_all_gain_GR =  with(events_all_gain, GRanges(chr, IRanges(start=start, end=end))) 
events_all_inv_GR =  with(events_all_inv, GRanges(chr, IRanges(start=start, end=end))) 


### Load the DGV and create genomic Ranges objects for each type of events
DGV2 <- DGV[which(DGV$samplesize > 1),]
loss_type <- c("deletion","loss","gain+loss") 
gain_type <- c("duplication","gain","gain+loss","tandem duplication")
inv_type <- c("inversion")
other_sv <- c("insertion","complex","novel sequence insertion","sequence alteration")
DGV_loss <-  DGV2[which(DGV2$variantsubtype %in% loss_type),]
DGV_gain <-  DGV2[which(DGV2$variantsubtype %in% gain_type),]
DGV_inv <-  DGV2[which(DGV2$variantsubtype %in% inv_type),]
DGV_oth <-  DGV2[which(DGV2$variantsubtype %in% other_sv),]
DGV2_loss <-  cbind(paste("chr",DGV_loss$chr,sep=""),DGV_loss[,c("start","end")])
DGV2_gain <-  cbind(paste("chr",DGV_gain$chr,sep=""),DGV_gain[,c("start","end")])
DGV2_inv <-  cbind(paste("chr",DGV_inv$chr,sep=""),DGV_inv[,c("start","end")])
DGV2_oth <-  cbind(paste("chr",DGV_oth$chr,sep=""),DGV_oth[,c("start","end")])
colnames(DGV2_loss) <-colnames(DGV2_gain) <-colnames(DGV2_inv) <-colnames(DGV2_oth) <- c("chr","start","end")
dgv_loss_ranges = with(DGV2_loss, GRanges(chr, IRanges(start=start, end=end))) 
dgv_gain_ranges = with(DGV2_gain, GRanges(chr, IRanges(start=start, end=end))) 
dgv_inv_ranges = with(DGV2_inv, GRanges(chr, IRanges(start=start, end=end))) 
dgv_compl_ranges = with(DGV2_oth, GRanges(chr, IRanges(start=start, end=end))) 


# find overlaps betwee events and junctions and the DGV (HC data)
hits_allJunct_dgv = GenomicAlignments::findOverlaps(dgv_loss_ranges,events_all_loss_GR)
overlaps_all <- pintersect(dgv_loss_ranges[queryHits(hits_allJunct_dgv),], events_all_loss_GR[subjectHits(hits_allJunct_dgv),])
percentOverlapA_all <- width(overlaps_all) / width(dgv_loss_ranges[queryHits(hits_allJunct_dgv),])
percentOverlapB_all <- width(overlaps_all) / width(events_all_loss_GR[subjectHits(hits_allJunct_dgv)])
remove_all_loss <- rownames(events_all_loss[unique(subjectHits(hits_allJunct_dgv[intersect(which(percentOverlapA_all > 0.5),which(percentOverlapB_all > 0.5))])),])

hits_allJunct_dgv = GenomicAlignments::findOverlaps(dgv_gain_ranges,events_all_gain_GR)
overlaps_all <- pintersect(dgv_gain_ranges[queryHits(hits_allJunct_dgv),], events_all_gain_GR[subjectHits(hits_allJunct_dgv),])
percentOverlapA_all <- width(overlaps_all) / width(dgv_gain_ranges[queryHits(hits_allJunct_dgv),])
percentOverlapB_all <- width(overlaps_all) / width(events_all_gain_GR[subjectHits(hits_allJunct_dgv)])
remove_all_gain <- rownames(events_all_gain[unique(subjectHits(hits_allJunct_dgv[intersect(which(percentOverlapA_all > 0.5),which(percentOverlapB_all > 0.5))])),])

hits_allJunct_dgv = GenomicAlignments::findOverlaps(dgv_inv_ranges,events_all_inv_GR)
overlaps_all <- pintersect(dgv_inv_ranges[queryHits(hits_allJunct_dgv),], events_all_inv_GR[subjectHits(hits_allJunct_dgv),])
percentOverlapA_all <- width(overlaps_all) / width(dgv_inv_ranges[queryHits(hits_allJunct_dgv),])
percentOverlapB_all <- width(overlaps_all) / width(events_all_inv_GR[subjectHits(hits_allJunct_dgv)])
remove_all_inv <- rownames(events_all_inv[unique(subjectHits(hits_allJunct_dgv[intersect(which(percentOverlapA_all > 0.5),which(percentOverlapB_all > 0.5))])),])

hits_allJunct_dgv = GenomicAlignments::findOverlaps(dgv_compl_ranges,events_all_complex_GR)
overlaps_all <- pintersect(dgv_compl_ranges[queryHits(hits_allJunct_dgv),], events_all_complex_GR[subjectHits(hits_allJunct_dgv),])
percentOverlapA_all <- width(overlaps_all) / width(dgv_compl_ranges[queryHits(hits_allJunct_dgv),])
percentOverlapB_all <- width(overlaps_all) / width(events_all_complex_GR[subjectHits(hits_allJunct_dgv)])
remove_all_compl <- rownames(events_all_complex[unique(subjectHits(hits_allJunct_dgv[intersect(which(percentOverlapA_all > 0.5),which(percentOverlapB_all > 0.5))])),])

remove_all <- c(remove_all_loss,remove_all_gain,remove_all_inv,remove_all_compl)
somatic_junctions_all_f1 <- somatic_junctions_all_f0[setdiff(rownames(somatic_junctions_all_f0),remove_all),]

## IDENTIFY ARTIFACTS BY LOCATION

highDensitySV <- function(sv_df,freq=10,size=100,offset=0){
highdens  <- list()
for(chr in paste("chr",c(1:22,"X"),sep="")) {
	message(chr)
	left <- sv_df$LeftPosition[which(sv_df$LeftChr == chr)]+offset
	names(left) <- rownames(sv_df)[which(sv_df$LeftChr == chr)]
	right <- sv_df$LeftPosition[which(sv_df$RightChr == chr)]+offset
	names(right) <- rownames(sv_df)[which(sv_df$RightChr == chr)]
	dat <- c(left,right)
	roundat <- round(dat/size)
	posids <- unique(names(which(sort(table(roundat)) > freq)))
	message(length(posids))
	for(i in posids){
		id<-paste(chr,paste(i,size/2-offset,sep=""),sep=":")
		highdens[[id]] <- unique(names(which(roundat == i)))
		}
	}
return(highdens)
}

filterRelated <- function(filteredJunctions,sv_df){

	samplerowname <- substr(rownames(sv_df[which(sv_df$Type == "complex"),]),0,24)
	relatedJunctions <- strsplit(sv_df[which(sv_df$Type == "complex"),"RelatedJunctions"],";")
	names(relatedJunctions) <- names(samplerowname) <- rownames(sv_df[which(sv_df$Type == "complex"),])

	filteredJunctions <- filteredJunctions
	results <-list()
	for(i in names(relatedJunctions)){
		results[[i]] <- length(intersect(paste(samplerowname[i],relatedJunctions[[i]],sep="."),filteredJunctions))
		}
	return(names(which(sort(unlist(results)) > 0)))
}

highdens<-highDensitySV(somatic_junctions_all_f1,freq=4,size=1000,offset=0)
highdens50<-highDensitySV(somatic_junctions_all_f1,freq=4,size=1000,offset=5)

highsamplefreq <- names(which(unlist(lapply(highdens,function(i) length(unique(substr(i,0,16))))) >=4))
hightumofreq <- names(which(unlist(lapply(highdens,function(i) length(table(substr(i,0,10)))))>1))
A<-unique(unlist(highdens[intersect(highsamplefreq,hightumofreq)]))

highsamplefreq50 <- names(which(unlist(lapply(highdens50,function(i) length(unique(substr(i,0,16))))) >=4))
hightumofreq50 <- names(which(unlist(lapply(highdens50,function(i) length(table(substr(i,0,10)))))>1))
B<-unique(unlist(highdens50[intersect(highsamplefreq50,hightumofreq50)]))
highDensfilter <- unique(c(A,B))

relatedJunctions <- filterRelated(highDensfilter,somatic_junctions_all_f1)

non_repeat_seq <- rownames(somatic_junctions_all_f1)[intersect(intersect(grep("Simple_repeat",somatic_junctions_all_f1$LeftRepeatClassification,invert=T),grep("Low_complexity",somatic_junctions_all_f1$LeftRepeatClassification,invert=T)),
	intersect(grep("Simple_repeat",somatic_junctions_all_f1$RightRepeatClassification,invert=T),grep("Low_complexity",somatic_junctions_all_f1$RightRepeatClassification,invert=T)))]


somatic_junctions_all_f2<-somatic_junctions_all_f1[setdiff(non_repeat_seq,c(highDensfilter,relatedJunctions)),]

# remove common variants found in DGV and split the data into the diferent TARGET tumors


all_somatic_junctions_all_f0 <- somatic_junctions_all_f2[grep("TARGET-10-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-0[4|9]",rownames(somatic_junctions_all_f2),perl=TRUE),]
aml_somatic_junctions_all_f0 <- somatic_junctions_all_f2[grep("TARGET-20-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-0[4|9]",rownames(somatic_junctions_all_f2),perl=TRUE),]
nbl_somatic_junctions_all_f0 <- somatic_junctions_all_f2[grep("TARGET-30-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-01",rownames(somatic_junctions_all_f2),perl=TRUE),]
os_somatic_junctions_all_f0 <- somatic_junctions_all_f2[grep("TARGET-40-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-01",rownames(somatic_junctions_all_f2),perl=TRUE),]
wt_somatic_junctions_all_f0 <- somatic_junctions_all_f2[grep("TARGET-50-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-01",rownames(somatic_junctions_all_f2),perl=TRUE),]
ccsk_somatic_junctions_all_f0 <- somatic_junctions_all_f2[grep("TARGET-51-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-01",rownames(somatic_junctions_all_f2),perl=TRUE),]

# Supplementary Table S3: nbl_somatic_junctions_all_f0


###################################
###################################
## Gene level enrichment analysis of SVs
# we are running enrichment analysis on each tumor and then focusing on NBL hits

feature_tab<- genes_tab
exons_tab <- exons_tab

upstr <- 100000         # cut off upstream mapping 
dnstr <- 25000          # cut off downstream mapping
copynumsize <- 2000000  # cut off for SV size in order to find overlaps
promoter=1000
offset=100
nonconsent <- c("TARGET-30-PASCHP","TARGET-30-PARCWT")
nbl_somatic_junctions_all_f00 <- nbl_somatic_junctions_all_f0[which(!nbl_somatic_junctions_all_f0$TARGET.USI %in% nonconsent), ]

results_NBL <- enrichSV(nbl_somatic_junctions_all_f00, feature_tab, exons_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100,copynumsize = 2000000)
results_ALL <- enrichSV(all_somatic_junctions_all_f0, feature_tab, exons_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100,copynumsize = 2000000)
results_AML <- enrichSV(aml_somatic_junctions_all_f0, feature_tab, exons_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100,copynumsize = 2000000)
results_OS <- enrichSV(os_somatic_junctions_all_f0, feature_tab, exons_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100,copynumsize = 2000000)
results_WT <- enrichSV(wt_somatic_junctions_all_f0, feature_tab, exons_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100,copynumsize = 2000000)
results_CCSK <- enrichSV(ccsk_somatic_junctions_all_f0, feature_tab, exons_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100,copynumsize = 2000000)


sv_disrupt <- merge2lists(results_NBL$disruptExonSamples,results_NBL$copynumSamples)
sv_proximal <- merge2lists(results_NBL$bothProximalSamples,results_NBL$disruptIntronSamples)
sv_summary_num <- sort(unlist(lapply(merge2lists(sv_disrupt,sv_proximal),length)))

# save filtered variants and additional annotation files
save(results_NBL,results_ALL,results_AML,results_OS,results_WT,results_CCSK,genes_tab,exons_tab,file="data/SV_analysis_Oct31_18.rda")


