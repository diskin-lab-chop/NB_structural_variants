## Code to generate figures for the manuscript :
# "Somatic structural variation targets neurodevelopmental genes and identifies SHANK2 as a tumor suppressor in neuroblastoma"
# Author: Gonzalo Lopez, PhD
# email: gonzalo.lopezgarcia@mssm.edu
# Date edited: Apr 16th, 2020


rm(list = ls(all.names = TRUE))
.rs.restartR()
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(devtools)
 library(VennDiagram)
 library(beeswarm)
 library(gplots)
library(data.table)
library(tidyr)
library(vioplot)
library(ComplexHeatmap)
library(cgdsr)

setwd("~/Box Sync/git/NB_structural_variants/")   # modify to local folder

source('R/my_stat_functions.R')
source('R/heatmap3.R')
source('R/01-nbl_somatic_SV_FUNCTIONS.R')


# clinical combined COG + TARGET
load("data/clinical_COG_20181129_plus_TARGET_20180331_plus_SCA_20180619_pheno.rda",verbose=T)
hr_unknown <- intersect(names(which(risk=="high")),names(which(mycna=="unknown")))


# load gene/exon info data UCSC/refseq
load("data/ucsc_hg19_refseq_genes_exon_df_Oct31_2018.rda",verbose=T)

# load SNP segmentation data data
load("data/cnv_segmentation_SNP_081318.rda",verbose=T)
# define NBL subtypes within the SNP dataset
snp_samples <- setdiff(unique(segment_snp$Sample),nonconsent)
HR_MNA_snp <- intersect(names(which(mycna == "amp")),snp_samples)
HR_NA_snp <- setdiff(intersect(names(which(risk == "high")),snp_samples),HR_MNA_snp)
INT_snp <- intersect(intrisk,snp_samples)
LOW_snp <- intersect(lowrisk,snp_samples)
LOWINT_snp <- c(INT_snp,LOW_snp)
HR_UNK_snp <- intersect(hr_unknown,snp_samples) # no need

## load SJBP, RDBP and CNBP filtered variants and segmentation breakpoint analyses 
load("data/SV_analysis_Oct31_18.rda",verbose=T)
load("data/BP_analysis_Nov15_19.rda",verbose=T)


# load dbgap metadata for CGI coverage figure
sraRunTab_wgs <- read.delim("data/SraRunTable_wgs_alltumors.txt",sep = "\t",as.is=TRUE)
sraRunTab_cgi <- sraRunTab_wgs[which(sraRunTab_wgs$Platform_s == "COMPLETE_GENOMICS"),]
tumor_libsize <- sraRunTab_cgi$MBases_l[grep("TARGET-30-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-01",sraRunTab_cgi$Sample_Name_s)]
names(tumor_libsize)<- substr(grep("TARGET-30-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-01",sraRunTab_cgi$Sample_Name_s,value=T),0,16)
blood_libsize <- sraRunTab_cgi$MBases_l[grep("TARGET-30-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-10",sraRunTab_cgi$Sample_Name_s)]
names(blood_libsize)<- substr(grep("TARGET-30-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-10",sraRunTab_cgi$Sample_Name_s,value=T),0,16)
ll<-list(tumor=tumor_libsize/3234.83,blood=blood_libsize/3234.83)

# Supplementary Figure S1A
vioplot(ll,las=1,border=c("black","black"),col=c("grey30","grey90"),names=c("Tumor","Blood"),ylab="Average depth")

# defoine NBL subtypes for WGS-CGI dataset
cgi_samples <- setdiff(substr(intersect(names(blood_libsize),names(tumor_libsize)),11,16),nonconsent)
HR_MNA <- intersect(names(which(mycna == "amp")),cgi_samples)
HR_NA <- intersect(intersect(names(which(risk=="high")),names(which(mycna !="amp" ))),cgi_samples)
INT <- intersect(names(which(risk=="intermediate")),cgi_samples)
LOW <- intersect(names(which(risk=="low")),cgi_samples)
LOWINT <- c(INT,LOW)

# load RNA-seq expression data for eQTL analysis and splicing 
load("data/kallisto_nbltpm_refseq_hg19_171.rda",verbose=TRUE)
sraRunTab_rna_nbl <- sraRunTab_rna[grep("TARGET-30-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-01",sraRunTab_rna$Sample_Name_s),]
srr <- sraRunTab_rna_nbl$Run_s[which(sraRunTab_rna_nbl$Center_Name_s == "NCI-KHAN")]
usi <- sraRunTab_rna_nbl$Sample_Name_s[which(sraRunTab_rna_nbl$Center_Name_s == "NCI-KHAN")]
r_expmat <- nblsym[,srr]
colnames(r_expmat) <- substr(usi,11,16)
r_expmat<-r_expmat[,setdiff(colnames(r_expmat),nonconsent)]
rna_samples <- colnames(r_expmat)

# load Human Exon array data for eQTL analysis and splicing 
idlist <- read.delim("data/usi_2targetID.txt",header=FALSE,as.is=TRUE)
convertid <- idlist$V1
names(convertid) <- idlist$V2
load("data/ln_mexp.rda",verbose=TRUE)
colnames(l_expmat) <- substr(convertid[colnames(l_expmat)],11,16)
l_expmat<-l_expmat[,setdiff(colnames(l_expmat),nonconsent)]
huex_samples <-  unname(colnames(l_expmat))

## Figure 1 clinical and copy number table
## combine datasets for clinical data table 
all_samples <- unique(c(rna_samples,huex_samples,snp_samples,cgi_samples))

wgs<-rna<-snp<-huex<-meth <- rep("white",length(all_samples))
names(wgs) <- names(rna) <- names(huex) <- names(snp) <- names(meth) <- all_samples

wgs[cgi_samples] <- "blue"
snp[snp_samples] <- "lightblue"
rna[rna_samples] <- "red"
huex[huex_samples] <- "purple"
stage_b <-stage
stage_b[which(stage_b == "4S")] <- 0

unknownmycn_snp <- genecn_snp["MYCN",intersect(hr_unknown,colnames(genecn_snp)) ]
unknownmycn_cgi <- genecn_cgi["MYCN",intersect(hr_unknown,colnames(genecn_cgi)) ]
mycna[names(unknownmycn_snp)[which(unknownmycn_snp > 1)]] <- "amp"
mycna[names(unknownmycn_snp)[which(unknownmycn_snp < 1)]] <- "non-amp"
mycna[intersect(hr_unk,colnames(genecn_cgi))] <- "non-amp"
hr_amp2 <- intersect(all_samples, names(which(mycna == "amp")))
hr_noa2 <- setdiff(intersect(all_samples, highrisk),hr_amp2)

group <- c(rep(0,length(intersect(all_samples,hr_amp2))),
			rep(1,length(intersect(all_samples,hr_noa2))),
				rep(2,length(intersect(all_samples,intrisk))),
					rep(3,length(intersect(all_samples,lowrisk))))
names(group)<-c(intersect(all_samples,hr_amp2),intersect(all_samples,hr_noa2),intersect(all_samples,intrisk),intersect(all_samples,lowrisk))
group_col <-group
group_col[which(group == 0)] <- "red" 
group_col[which(group == 1)] <- "orange" 
group_col[which(group == 2)] <- "darkgreen" 
group_col[which(group == 3)] <- "green" 

# Phonotype data.frame with clinical covariates
Pheno2 <- rbind(Pheno["gender_col",names(group_col)],age_col[names(group_col)],Pheno[c("vital_col","mycna_col","stage_col","risk_col"),names(group_col)],group_col)


# order samples for figure 1A according to data availability $ n=1025
a<- stage_b[names(sort(stage_b[names(group)],decreasing=T))]
		a<- a[names(sort(huex))]
			a<- a[names(sort(rna))]
				a<- a[names(sort(snp))]
					a<- a[names(sort(wgs))]
final_order <- names( sort(group[names(a)]))
final_order <- final_order[which(final_order %in% names(group))]


# plot mock heatmaps with phenotype bars for the figure 1A composition
fakemat <- matrix(rnorm(5*length(final_order)),ncol=length(final_order),nrow=5)
colnames(fakemat) <- final_order
Data_col <- rbind(meth,huex,rna,snp,wgs)
SCA_col2 <- matrix(ncol=length(final_order),nrow=nrow(SCA_col))
colnames(SCA_col2) <- final_order
rownames(SCA_col2) <- rownames(SCA_col)
SCA_col2[]<- "lightgrey" 
SCA_col2[,intersect(colnames(SCA_col),final_order) ] <- SCA_col[,intersect(colnames(SCA_col),final_order) ] 

png("Figures/Figure1/fig1a_sample/heat_data_10_17_18.png",height=500,width=1000)
heatmap.3(fakemat[,final_order], ColSideColors = t(Data_col[,final_order]),Colv=FALSE,Rowv=FALSE,lhei=c(0.5,1))
dev.off()

png("Figures/Figure1/fig1a_sample/heat_sca_10_17_18.png",height=600,width=1000)
heatmap.3(fakemat[,final_order], ColSideColors = t(SCA_col2[,final_order]),Colv=FALSE,Rowv=FALSE,lhei=c(0.5,1))
dev.off()

png("Figures/Figure1/fig1a_sample/heat_clinical_10_17_18.png",height=600,width=1000)
heatmap.3(fakemat[,final_order], ColSideColors = t(Pheno2[,final_order]),Colv=FALSE,Rowv=FALSE,lhei=c(0.5,1))
dev.off()

# write CGI segmentation and attributes for IGV visualization 
order_cgi<-1:length(cgi_samples)
names(order_cgi) <- final_order[which(final_order %in% segment_cgi$Sample)]
#write.table(segment_cgi[order(order_cgi[segment_cgi$Sample]),],row.names=F,quote=F,sep="\t")
cgi_attributes <- cbind(names(order_cgi),group[names(order_cgi)],mycna[names(order_cgi)],stage[names(order_cgi)],gender[names(order_cgi)],order_cgi)
colnames(cgi_attributes) <- c("TRACK_ID","group","mycna","stage","gender","order")
#write.table(cgi_attributes,row.names=F,quote=F,sep="\t")

# write SNP segmentation and attributes for IGV visualization 
order_snp<-1:length(final_order[which(final_order %in% segment_snp$Sample)])
names(order_snp) <- final_order[which(final_order %in% segment_snp$Sample)]
#write.table(segment_snp[order(order_snp[segment_snp$Sample]),],row.names=F,quote=F,sep="\t")
snp_attributes <- cbind(names(order_snp),group[names(order_snp)],mycna[names(order_snp)],stage[names(order_snp)],gender[names(order_snp)],order_snp)
colnames(snp_attributes) <- c("TRACK_ID","group","mycna","stage","gender","order")
#write.table(snp_attributes,row.names=F,quote=F,sep="\t")

# write supplementary table 2 (samples and clinical info)
data_clinical <- cbind(group[all_samples],stage[all_samples],mycna[all_samples],risk[all_samples],age[all_samples],status[all_samples],t(Data_col[,all_samples]))
#write.table(data_clinical,file="Suppl_table_02.txt",sep="\t",quote=F)

####################################################


# supplementary figure S3A
ll_sv <- list(
  NBL=results_NBL$sv_df,
  ALL=results_ALL$sv_df,
  AML=results_AML$sv_df,
  OS=results_OS$sv_df,
  WT=results_WT$sv_df
)

svTypeCols <- c("lightgrey","blue","purple","green","yellow","red")
names(svTypeCols)<- c("complex","deletion","distal-duplication","interchromosomal","inversion","tandem-duplication")

ll_type <- ll_type_m <-list()
for(i in names(ll_sv)){
  ll_sv[[i]][which(ll_sv[[i]]$Type == "probable-inversion"),"Type"] <- "inversion"
  ll_sv[[i]][which(ll_sv[[i]]$Type == "distal-duplication-by-mobile-element"),"Type"] <- "distal-duplication"
  ll_type[[i]]<- table(ll_sv[[i]]$Type)/sum(table(ll_sv[[i]]$Type))
  ll_type_m[[i]]<- table(ll_sv[[i]]$Type)/length(unique(as.character(ll_sv[[i]]$TARGET.USI)))
}

dev.new()
par(mfrow=c(2,1),mar=c(4,4,0,1))
datatab <- do.call(cbind,ll_type_m)
barplot(datatab,col=svTypeCols,las =1,beside = T,cex.axis=1.3,names=c("","","","",""))
legend("topleft",names(svTypeCols),fill=svTypeCols,bty='n',cex=1.3)
datatab <- do.call(cbind,ll_type)
barplot(datatab,col=svTypeCols,las =1,cex.axis=1.3,cex.names=1.3)


# supplementary figure S3B
sv_df <- results_NBL$sv_df
sv_df[,"TARGET.USI"] <- substr(sv_df$TARGET.USI ,11,16)
sv_df[which(sv_df$Type == "probable-inversion"),"Type"] <- "inversion"
sv_df[which(sv_df$Type == "distal-duplication-by-mobile-element"),"Type"] <- "distal-duplication"


ll_stype <- list(
  MNA = table(sv_df$Type[which(sv_df$TARGET.USI %in% HR_MNA)])/ length(which(sv_df$TARGET.USI %in% HR_MNA)),
  HR_NA = table(sv_df$Type[which(sv_df$TARGET.USI %in% HR_NA)])/length(which(sv_df$TARGET.USI %in% HR_NA)) ,
  LINT = table(sv_df$Type[which(sv_df$TARGET.USI %in% c(LOW,INT))])/length(which(sv_df$TARGET.USI %in% c(LOW,INT)))
)
ll_stype_m <- list(
  MNA = table(sv_df$Type[which(sv_df$TARGET.USI %in% HR_MNA)])/length(HR_MNA),
  HR_NA = table(sv_df$Type[which(sv_df$TARGET.USI %in% HR_NA)])/length(HR_NA) ,
  LINT = table(sv_df$Type[which(sv_df$TARGET.USI %in% c(LOW,INT))])/length(c(LOW,INT))
)
par(mfrow=c(2,1),mar=c(4,4,0.5,1))
datatab <- do.call(cbind,ll_stype_m)
barplot(datatab,col=svTypeCols,las =1,beside = T,cex.axis=1.3,names=c("","",""))
legend("topright",names(svTypeCols),fill=svTypeCols,bty='n')
datatab <- do.call(cbind,ll_stype)
barplot(datatab,col=svTypeCols,las =1,cex.axis=1.3,cex.names=1.3)


## join exonic and CNn into coding and intronic and proximal into noncoding
sv_exonic <- lapply(results_NBL$disruptExonSamples,function(x) substr(x,11,16))
sv_intronic <- lapply(results_NBL$disruptIntronSamples,function(x) substr(x,11,16))
sv_proximal <- lapply(results_NBL$bothProximalSamples,function(x) substr(x,11,16))
sv_copynum<- lapply(results_NBL$copynumSamples,function(x) substr(x,11,16))
sv_coding <- merge2lists(sv_exonic,sv_copynum)
sv_noncod <- merge2lists(sv_proximal,sv_intronic)
sv_summary <- merge2lists(sv_coding,sv_noncod)

j_exonic <- results_NBL$disruptExonJunctions
j_intronic <- results_NBL$disruptIntronJunctions
j_proximal <- results_NBL$bothProximalJunctions
j_copynum<- results_NBL$copynumJunctions
j_coding <- merge2lists(j_exonic,j_copynum)
j_noncod <- merge2lists(j_proximal,j_intronic)
j_summary <- merge2lists(j_coding,j_noncod)


# plot density of DiscordantMatePairAlignments in lieu of VAF as requested by reviewer 1
# (Supplementary Figure S1B)
dev.new()
plot(density(log2(ll_sv$AML$DiscordantMatePairAlignments)),col="white",las=1,xlim=c(0,15),main="",xaxt='n',xlab="",ylab=)
lines(density(log2(ll_sv$ALL$DiscordantMatePairAlignments)),col="salmon",lty=3,lwd=2)
lines(density(log2(ll_sv$AML$DiscordantMatePairAlignments)),col="pink",lty=3,lwd=2)
lines(density(log2(ll_sv$OS$DiscordantMatePairAlignments)),col="brown",lty=3,lwd=2)
lines(density(log2(ll_sv$WT$DiscordantMatePairAlignments)),col="green",lty=3,lwd=2)
lines(density(log2(ll_sv$NBL$DiscordantMatePairAlignments)),col="blue",lty=1,lwd=2)
legend("topright",c("ALL","AML","OS","WT","NBL"),lty=c(3,3,3,3,1),lwd=2,col=c("pink","salmon","brown","green","blue"),bty='n')
axis(1,at=seq(0,16,2),labels=2^seq(0,16,2),las=2)

## plot hitogram sizes per type (Supplementary Figure S1C-D)
sv_df_5t <- rbind(results_ALL$sv_df,results_AML$sv_df,results_NBL$sv_df,results_OS$sv_df,results_WT$sv_df)
sameChr<- which(sv_df_5t$Type %in% c("deletion","inversion","tandem-duplication","probable-inversion"))
sv_df_same_chr <- sv_df_5t[sameChr,]
sizes <- sv_df_same_chr$RightPosition[]-sv_df_same_chr$LeftPosition
par(mfrow=c(1,3))
h2 <- hist(log10(sizes[which(sv_df_same_chr$Type =="deletion")]+1),
	col=rgb(0,0,1,0.5),breaks=40,border=F,xlim=c(0,8),main="Deletion",las=1)
h3 <- hist(log10(sizes[which(sv_df_same_chr$Type == "tandem-duplication")]+1),
	col=rgb(0.62,0.12,0.94,0.5),breaks=50,border=F,xlim=c(0,8),main="Tandem-duplication",las=1)
h4 <- hist(log10(sizes[which(sv_df_same_chr$Type %in% c("inversion","probable-inversion"))]+1),
	col=rgb(1,0.84,0,0.5),breaks=50,border=F,xlim=c(0,8),main="Inversion",las=1)



## Figure 2A stack barplot ##
sample_order <- rep(0,length(unique(as.character(sv_df$TARGET.USI))))
names(sample_order) <- unique(as.character(sv_df$TARGET.USI))

types <- unique(sv_df$Type)
byType <- sapply(names(sample_order),function(i) table(sv_df[which(sv_df$TARGET.USI == i),"Type"])  )
byType<-list()
for(i in cgi_samples){
	byType[[i]] <- rep(0,length(types))
	names(byType[[i]]) <- types
	dat <- table(sv_df[which(sv_df$TARGET.USI == i),"Type"])
	byType[[i]][intersect(names(dat),types)] <- dat[intersect(names(dat),types)]
	}
totSV <- unlist(lapply(byType,sum))

sample_order[names(totSV)] <- totSV
sample_order <- c(sort(sample_order[HR_MNA]),sort(sample_order[HR_NA]),sort(sample_order[INT]),sort(sample_order[LOW]))


svTypeCols <- c("lightgrey","blue","purple","green","yellow","red")
names(svTypeCols)<- c("complex","deletion","distal-duplication","interchromosomal","inversion","tandem-duplication")

svBarPlorData <- matrix(ncol=length(sample_order),nrow=length(types))
colnames(svBarPlorData) <- names(sample_order)
rownames(svBarPlorData) <- types
svBarPlorData[]<-0
for(i in colnames(svBarPlorData)) svBarPlorData[names(byType[[i]]),i] <- byType[[i]]
colnames(svBarPlorData) <- substr(names(sample_order),11,16)
svBarPlorData<-svBarPlorData[intersect(names(svTypeCols),rownames(svBarPlorData)),]
par(family = "Arial",lwd = 0.3)
barplot(svBarPlorData,las=2,cex.names=0.7,col=svTypeCols[rownames(svBarPlorData)],lwd=1,cex.axis=1.3)
par(family = "Arial",lwd = 1)
abline(h=seq(0,3500,50),lty=3,lwd=0.7,col="grey")
legend("topright",rev(rownames(svBarPlorData)),fill=rev(svTypeCols[rownames(svBarPlorData)]),bty="n",cex=1.3,ncol=2)


#################################################
# overal overlap between segdata and SV junctions
# # # # # # # # # # # # # # # # # # # # # # # # # # 

# 
sv_noa <- sv_df[which(sv_df$TARGET.USI %in% HR_NA),]
sv_amp <- sv_df[which(sv_df$TARGET.USI %in% HR_MNA),]
sv_noa_types <- sv_amp_types <- list()
 types_temp <- rep(0,length(types))
 names(types_temp)<- types
for(chr in paste("chr",c(1:22,"X"),sep="")){
	sv_noa_types[[chr]] <- sv_amp_types[[chr]] <-types_temp
	noachr <- intersect(which(sv_noa$LeftChr == chr),which(sv_noa$RightChr == chr))
	ampchr <- intersect(which(sv_amp$LeftChr == chr),which(sv_amp$RightChr == chr))

	noatypescol <- sv_noa$Type[noachr]
	noatypescol[which(noatypescol == "probable-inversion")] <- "inversion"
	sv_noa_types[[chr]][names(table(noatypescol))] <-   table(noatypescol)

	amptypescol <- sv_noa$Type[ampchr]
	amptypescol[which(amptypescol == "probable-inversion")] <- "inversion"
	sv_amp_types[[chr]][names(table(amptypescol))] <-   table(amptypescol)
	}


tempsampl <- rep(0,length(cgi_samples))
names(tempsampl) <-cgi_samples
llsv <- list()
for(chr in paste("chr",c(1:22,"X"),sep="")) llsv[["all"]][[paste(chr,"S4na")]] <- tempsampl[HR_NA]
for(chr in paste("chr",c(1:22,"X"),sep="")) llsv[["all"]][[paste(chr,"MNa")]] <- tempsampl[HR_MNA]

par(mfrow=c(1,6),mar=c(5,1,3,1))
types <- c("complex","deletion","distal-duplication","inversion","interchromosomal","tandem-duplication")
results <- list()
for(type in types){
	tempsampl <- rep(0,length(cgi_samples))
	names(tempsampl) <-cgi_samples
	typelist<-c()
	llsv[[type]]<-list()
	for(chr in paste("chr",c(1:22,"X"),sep="")){
		if(type == "inversion"){ typelist <-c("inversion","probable-inversion")
		}else{typelist = type}
		sv_ids <- c(setdiff(intersect(which(sv_df$RightChr == chr),which(sv_df$Type %in% typelist)),intersect(which(sv_df$LeftChr == chr),which(sv_df$Type %in% typelist)) ),
			intersect(which(sv_df$LeftChr == chr),which(sv_df$Type %in% typelist)))
		chrdat <- na.omit(table(as.character(sv_df$TARGET.USI[sv_ids])))
		tempsampl[]<-0
		tempsampl[names(chrdat)] <- chrdat
		wt <- wilcox.test(tempsampl[HR_NA],tempsampl[HR_MNA])
		site <- sign(mean(rank(tempsampl)[HR_NA])-mean(rank(tempsampl)[HR_MNA]))
		results[[type]][[chr]] <- site*(-log10(wt$p.value))
		llsv[[type]][[paste(chr,"S4na")]]<-tempsampl[HR_NA]
		llsv[[type]][[paste(chr,"MNa")]]<-tempsampl[HR_MNA]
		llsv[["all"]][[paste(chr,"S4na")]] <- llsv[["all"]][[paste(chr,"S4na")]]+tempsampl[HR_NA]
		llsv[["all"]][[paste(chr,"MNa")]] <- llsv[["all"]][[paste(chr,"MNa")]]+tempsampl[HR_MNA]
		}
	results[[type]] <-rev(results[[type]])
	data <- unlist(results[[type]])
	color <- data
	color[] <- "white"
	color[which(data < 0)] <- "red"
	color[which(data > 0)] <- "orange"
	xmin<-min(results[[type]])
	if(min(results[[type]]) > -5) xmin <- -5
	xmax<-max(results[[type]])
	if(max(results[[type]]) < 5) xmax <- 5
	barplot(unlist(results[[type]]),las=1,names="",col=color,main=type,horiz=T,xlim=c(xmin,xmax))
	abline(v=c(-1.3,0,1.3),lty=c(3,1,3))
}
write.table(do.call(cbind,results),sep="\t",file=)
#write.table(do.call(cbind,results),sep="\t",file="~/Box Sync/My_CHOP/SV_paper_V2/Figures/Figure1/suppl_1/sv_difference_byChrbyType_pvalues.txt")

# Figures 2D-G
a<-b<-c<-d<-e<-f<-g<-h<-i<-j<-k<-l<-wtab<-wtcd<-wtef<-wtgh<-wtij<-wtkl<-list()
for(chr in paste("chr",c(1:22,"X"),sep="")){ 
	a[[chr]] <- mean(sort(llsv[["interchromosomal"]][[paste(chr,"MNa")]])[3:27])
	b[[chr]] <- mean(sort(llsv[["interchromosomal"]][[paste(chr,"S4na")]])[7:70])
	wtab[[chr]] <- wilcox.test(sort(llsv[["interchromosomal"]][[paste(chr,"MNa")]])[3:27],sort(llsv[["interchromosomal"]][[paste(chr,"S4na")]])[7:70])$p.value
	c[[chr]] <- mean(sort(llsv[["tandem-duplication"]][[paste(chr,"MNa")]])[3:27])
	d[[chr]] <- mean(sort(llsv[["tandem-duplication"]][[paste(chr,"S4na")]])[7:70])
	wtcd[[chr]] <- wilcox.test(sort(llsv[["tandem-duplication"]][[paste(chr,"MNa")]])[3:27],sort(llsv[["tandem-duplication"]][[paste(chr,"S4na")]])[7:70])$p.value
	e[[chr]] <- mean(sort(llsv[["deletion"]][[paste(chr,"MNa")]])[3:27])
	f[[chr]] <- mean(sort(llsv[["deletion"]][[paste(chr,"S4na")]])[7:70])
	wtef[[chr]]<- wilcox.test(sort(llsv[["deletion"]][[paste(chr,"MNa")]])[3:27],sort(llsv[["deletion"]][[paste(chr,"S4na")]])[7:70])$p.value
	g[[chr]] <- mean(sort(llsv[["inversion"]][[paste(chr,"MNa")]])[3:27])
	h[[chr]] <- mean(sort(llsv[["inversion"]][[paste(chr,"S4na")]])[7:70])
	wtgh[[chr]]<- wilcox.test(sort(llsv[["inversion"]][[paste(chr,"MNa")]])[3:27],sort(llsv[["inversion"]][[paste(chr,"S4na")]])[7:70])$p.value
	i[[chr]] <- mean(sort(llsv[["complex"]][[paste(chr,"MNa")]])[3:27])
	j[[chr]] <- mean(sort(llsv[["complex"]][[paste(chr,"S4na")]])[7:70])
	wtij[[chr]]<- wilcox.test(sort(llsv[["complex"]][[paste(chr,"MNa")]])[3:27],sort(llsv[["complex"]][[paste(chr,"S4na")]])[7:70])$p.value
	k[[chr]] <- mean(sort(llsv[["all"]][[paste(chr,"MNa")]])[3:27])
	l[[chr]] <- mean(sort(llsv[["all"]][[paste(chr,"S4na")]])[7:70])
	wtkl[[chr]]<- wilcox.test(sort(llsv[["all"]][[paste(chr,"MNa")]])[3:27],sort(llsv[["all"]][[paste(chr,"S4na")]])[7:70])$p.value
	}
par(mfrow=c(4,2))
barplot(rbind(unlist(a),unlist(b)),beside=T,las=2,col=c("red","orange"),ylim=c(0,1.5))
barplot(rbind(unlist(c),unlist(d)),beside=T,las=2,col=c("red","orange"),ylim=c(0,8),yaxt='n')
axis(2,at=seq(0,8,0.5))
barplot(rbind(unlist(e),unlist(f)),beside=T,las=2,col=c("red","orange"),ylim=c(0,8),yaxt='n')
axis(2,at=seq(0,8,0.5))
barplot(rbind(unlist(g),unlist(h)),beside=T,las=2,col=c("red","orange"),ylim=c(0,11))
axis(2,at=seq(0,11,0.5))
barplot(rbind(unlist(i),unlist(j)),beside=T,las=2,col=c("red","orange"),ylim=c(0,43),yaxt='n')
axis(2,at=seq(0,43,1))
barplot(rbind(unlist(k),unlist(l)),beside=T,las=2,col=c("red","orange"),ylim=c(0,78),yaxt='n')
axis(2,at=seq(0,78,1),las=1)



a <- sv_df$LeftChr[intersect(which(sv_df$Type == "complex"),which(sv_df$TARGET.USI %in% HR_MNA))]
b <- sv_df$RightChr[intersect(which(sv_df$Type == "complex"),which(sv_df$TARGET.USI %in% HR_MNA))]
length(intersect(which(a != "chr2"),which(b != "chr2")))/length(HR_MNA)

a <- sv_df$LeftChr[intersect(which(sv_df$Type == "complex"),which(sv_df$TARGET.USI %in% HR_NA))]
b <- sv_df$RightChr[intersect(which(sv_df$Type == "complex"),which(sv_df$TARGET.USI %in% HR_NA))]
length(intersect(which(a != "chr2"),which(b != "chr2")))/length(HR_NA)

###############
# Figures 2B
breaks <- results_BP$breaks
mnabps <- rep(0,length(HR_MNA));names(mnabps) <- HR_MNA; mnabps[names(table(breaks$sample)[HR_MNA])] <- table(breaks$sample)[HR_MNA]
nabps <- rep(0,length(HR_NA));names(nabps) <- HR_NA; nabps[names(table(breaks$sample)[HR_NA])] <- table(breaks$sample)[HR_NA]
intbps <- rep(0,length(INT));names(intbps) <- INT; intbps[names(table(breaks$sample)[INT])] <- table(breaks$sample)[INT]
lowbps <- rep(0,length(LOW));names(lowbps) <- LOW; lowbps[names(table(breaks$sample)[LOW])] <- table(breaks$sample)[LOW]
bpfreq<-c(sort(mnabps),sort(nabps),sort(intbps),sort(lowbps))
barplot(bpfreq,col=c(rep("red",length(HR_MNA)),rep("orange",length(HR_NA)),rep("darkgreen",length(INT)),rep("green",length(LOW)) ),las=1,border=NA) 

# Figures 2C
breaks <- results_BP_snp$breaks
mnabps <- rep(0,length(HR_MNA_snp));names(mnabps) <- HR_MNA_snp; mnabps[names(table(breaks$sample)[HR_MNA_snp])] <- table(breaks$sample)[HR_MNA_snp]
nabps <- rep(0,length(HR_NA_snp));names(nabps) <- HR_NA_snp; nabps[names(table(breaks$sample)[HR_NA_snp])] <- table(breaks$sample)[HR_NA_snp]
intbps <- rep(0,length(INT_snp));names(intbps) <- INT_snp; intbps[names(table(breaks$sample)[INT_snp])] <- table(breaks$sample)[INT_snp]
lowbps <- rep(0,length(LOW_snp));names(lowbps) <- LOW_snp; lowbps[names(table(breaks$sample)[LOW_snp])] <- table(breaks$sample)[LOW_snp]
bpfreq<-c(sort(mnabps),sort(nabps),sort(intbps),sort(lowbps))
barplot(bpfreq,col=c(rep("red",length(HR_MNA_snp)),rep("orange",length(HR_NA_snp)),rep("darkgreen",length(INT_snp)),rep("green",length(LOW_snp)) ),las=1,border=NA) 

# Figures 2H
a <- b <- wtab <- list()
tempamp <- rep(0,length(HR_MNA))
names(tempamp) <-HR_MNA
tempnoa <- rep(0,length(HR_NA))
names(tempnoa) <- HR_NA
for(chr in paste("chr",c(1:22,"X"),sep="")){ 
	a[[chr]] <- tempamp
	b[[chr]] <- tempnoa
 	resa <- table(results_BP$breaks$sample[intersect(which(results_BP$breaks$sample %in% HR_MNA),which(results_BP$breaks$chr == chr))])
 	resb <- table(results_BP$breaks$sample[intersect(which(results_BP$breaks$sample %in% HR_NA),which(results_BP$breaks$chr == chr))])
	a[[chr]][names(resa)] <- resa
	b[[chr]][names(resb)] <- resb
	wtab[[chr]]<-wilcox.test(a[[chr]],b[[chr]])$p.value
	a[[chr]] <- mean(sort(a[[chr]])[3:27])
	b[[chr]] <- mean(sort(b[[chr]])[7:70])
	}
par(mfrow=c(1,1),family = "Arial")
barplot(rbind(unlist(a),unlist(b)),beside=T,las=2,col=c("red","orange"),ylim=c(0,8),yaxt='n')
axis(2,at=seq(0,8,0.5))

# Figures 2I
a <- b <- wtab <- list()
tempamp <- rep(0,length(HR_MNA_snp))
names(tempamp) <-HR_MNA_snp
tempnoa <- rep(0,length(HR_NA_snp))
names(tempnoa) <- HR_NA_snp
for(chr in paste("chr",c(1:22,"X"),sep="")){ 
	a[[chr]] <- tempamp
	b[[chr]] <- tempnoa
 	resa <- table(results_BP_snp$breaks$sample[intersect(which(results_BP_snp$breaks$sample %in% HR_MNA_snp),which(results_BP_snp$breaks$chr == chr))])
 	resb <- table(results_BP_snp$breaks$sample[intersect(which(results_BP_snp$breaks$sample %in% HR_NA_snp),which(results_BP_snp$breaks$chr == chr))])
	a[[chr]][names(resa)] <- resa
	b[[chr]][names(resb)] <- resb
	wtab[[chr]]<-wilcox.test(a[[chr]],b[[chr]])$p.value
	a[[chr]] <- mean(sort(a[[chr]])[23:215])
	b[[chr]] <- mean(sort(b[[chr]])[42:380])
	}
barplot(rbind(unlist(a),unlist(b)),beside=T,las=2,col=c("red","orange"))
axis(2,at=seq(0,8,0.5))

## Show localization of most chr2 SV
## repeat chr2 genomic SVs
# Supplementary Figure S4C
ampchr2 <- intersect(which(sv_df$TARGET.USI %in% HR_MNA),which(sv_df$LeftChr == "chr2"))
hist(as.numeric(sv_df$LeftPosition[ampchr2]),breaks=100,col=rgb(1,0,0,0.5),border=rgb(1,0,0),las=1)
noachr2 <- intersect(which(sv_df$TARGET.USI %in% HR_NA),which(sv_df$LeftChr == "chr2"))
hist(as.numeric(sv_df$LeftPosition[noachr2]),breaks=100,col=rgb(1,0.65,0,0.5),border=rgb(1,0.65,0),add=T,las=1)


#######################################
### Obtain shortlisted genes for oncoprint and eQTL analyses

### Create a table of frequently altered genes by read depth
bp_summary <-merge2lists(results_BP$breakSamples, results_BP$proximalSamples)
bp_coding <- results_BP$breakSamples
bp_coding_4tum <-  merge2lists(lapply(lapply(results_BP_4tum$breakSamples, function(x) substr(x,11,16)),unique)	,
                  merge2lists(lapply(results_CN_4tum$gene_ampl, function(x) unique(substr(x,11,16))),
                              lapply(results_CN_4tum$gene_ddel, function(x) unique(substr(x,11,16)))))

bp_coding_4tum <- lapply(results_BP_4tum$breakSamples, function(x) substr(x,11,16))

bp_noncod<-results_BP$proximalSamples
genes_breaks <- sort(unlist(lapply(bp_coding,length)),decreasing=T) 
genes_bpprox  <- sort(unlist(lapply(bp_noncod,length)),decreasing=T) 
genes_bpsummary  <- sort(unlist(lapply(bp_summary,length)),decreasing=T) 

bp_summary_snp <-merge2lists(results_BP_snp$breakSamples, results_BP_snp$proximalSamples)
bp_coding_snp<-results_BP_snp$breakSamples
bp_noncod_snp<-results_BP_snp$proximalSamples
genes_breaks_snp <- sort(unlist(lapply(bp_coding_snp,length)),decreasing=T) 
genes_bpprox_snp  <- sort(unlist(lapply(bp_noncod_snp,length)),decreasing=T) 
genes_bpsummary_snp  <- sort(unlist(lapply(bp_summary_snp,length)),decreasing=T) 

ampl_ALL <- results_CN_4tum$gene_ampl[which( unlist(lapply(results_CN_4tum$gene_ampl,length) >1) )]
ampl_ALL_num <- unlist(lapply(ampl_ALL,length))
ampl_ALL_samples <- unlist(lapply(ampl_ALL,paste,collapse=" "))
ddel_ALL <- results_CN_4tum$gene_ddel[which( unlist(lapply(results_CN_4tum$gene_ddel,length) >1) )]
ddel_ALL_num <- unlist(lapply(ddel_ALL,length))
ddel_ALL_samples <- unlist(lapply(ddel_ALL,paste,collapse=" "))
tail(sort(unlist(lapply(ddel_ALL,function(x) length(unique(substr(x,11,16)))))),100)


ampl_WGS <- results_CN$gene_ampl[which( unlist(lapply(results_CN$gene_ampl,length) >1) )]
ampl_WGS_num <- unlist(lapply(ampl_WGS,length))
ampl_WGS_samples <- unlist(lapply(ampl_WGS,paste,collapse=" "))
ddel_WGS <- results_CN$gene_ddel[which( unlist(lapply(results_CN$gene_ddel,length) >1) )]
ddel_WGS_num <- unlist(lapply(ddel_WGS,length))
ddel_WGS_samples <- unlist(lapply(ddel_WGS,paste,collapse=" "))

ampl_SNP <- results_CN_snp$gene_ampl[which( unlist(lapply(results_CN_snp$gene_ampl,length) >1) )]
ampl_SNP_num <- unlist(lapply(ampl_SNP,length))
ampl_SNP_samples <- unlist(lapply(ampl_SNP,paste,collapse=" "))
ddel_SNP <- results_CN_snp$gene_ddel[which( unlist(lapply(results_CN_snp$gene_ddel,length) >1) )]
ddel_SNP_num <- unlist(lapply(ddel_SNP,length))
ddel_SNP_samples <- unlist(lapply(ddel_SNP,paste,collapse=" "))

all_cn <- merge2lists(merge2lists(ampl_WGS,ddel_WGS),merge2lists(ampl_SNP,ddel_SNP))

bp_summary_combined <- merge2lists(merge2lists(bp_summary_snp,bp_summary),all_cn)
bpgenes <- names(which(unlist(lapply(bp_summary_combined,length))>0))
genes_bpsummary_both  <- sort(unlist(lapply(bp_summary_combined,length))) 

coding_WGS_samples <- unlist(lapply(bp_coding,paste,collapse=" "))[bpgenes]
proximal_WGS_samples <- unlist(lapply(bp_noncod,paste,collapse=" "))[bpgenes]
coding_SNP_samples <- unlist(lapply(bp_coding_snp,paste,collapse=" "))[bpgenes]
proximal_SNP_samples <- unlist(lapply(bp_noncod_snp,paste,collapse=" "))[bpgenes]


bpgene_tab <- cbind(genes_tab[bpgenes,c(1,2,3,5)],genes_breaks[bpgenes],genes_bpprox[bpgenes],genes_bpsummary[bpgenes],coding_WGS_samples[bpgenes],proximal_WGS_samples[bpgenes],   
		ampl_WGS_num[bpgenes],ampl_WGS_samples[bpgenes],ddel_WGS_num[bpgenes],ddel_WGS_samples[bpgenes],
		genes_breaks_snp[bpgenes],genes_bpprox_snp[bpgenes],genes_bpsummary_snp[bpgenes],coding_SNP_samples[bpgenes],proximal_SNP_samples[bpgenes],
		ampl_SNP_num[bpgenes],ampl_SNP_samples[bpgenes],ddel_SNP_num[bpgenes],ddel_SNP_samples[bpgenes],
		genes_bpsummary_both[bpgenes])
colnames(bpgene_tab) <- c(colnames(genes_tab[bpgenes,c(1,2,3,5)]),
		"coding_CGI","proximal_CGI","summary_CGI","coding_CGI_samples","proximal_CGI_samples",
		"ampl_CGI","ampl_CGI_samples","ddel_CGI","ddel_CGI_samples",
		"coding_SNP","proximal_SNP","summary_SNP","coding_SNP_samples","proximal_SNP_samples",
		"ampl_SNP","ampl_SNP_samples","ddel_SNP","ddel_SNP_samples",
		"total_summary")

rownames(bpgene_tab) <- bpgenes
#write.table(bpgene_tab,file="Supplemental_Table_S6_1.txt",sep="\t",quote=F)

## create table based on SV only 
genes_copynum <- sort(unlist(lapply(sv_copynum,length)),decreasing=T)
genes_exonic <- sort(unlist(lapply(sv_exonic,length)),decreasing=T)
genes_proximal <-  sort(unlist(lapply(sv_proximal,length)),decreasing=T) 
genes_intronic <- sort(unlist(lapply(sv_intronic,length)),decreasing=T) 
genes_coding <- sort(unlist(lapply(sv_coding,length)),decreasing=T) 
genes_noncod <- sort(unlist(lapply(sv_noncod,length)),decreasing=T) 
genes_summary <- unlist(lapply(sv_summary,length))

a<-results_NBL$bothProximalTypesMat 
b<-results_NBL$copynumTypesMat
c<-results_NBL$disruptExonTypesMat
d<-results_NBL$disruptIntronTypesMat
genes <- unique(c(rownames(a),rownames(b),rownames(c),rownames(d)))
colns <- unique(c(colnames(a),colnames(b),colnames(c),colnames(d)))
allType <- data.frame(matrix(ncol=length(colns),nrow=length(genes)))
rownames(allType) <- genes
colnames(allType) <- colns
allType[]<-0
allType[rownames(a),colnames(a)] <- allType[rownames(a),colnames(a)] +a
allType[rownames(b),colnames(b)] <- allType[rownames(b),colnames(b)] +b
allType[rownames(c),colnames(c)] <- allType[rownames(c),colnames(c)] +c
allType[rownames(d),colnames(d)] <- allType[rownames(d),colnames(d)] +d

allgenes <- names(which(sort(genes_summary,decreasing=T) > 0))

samples_coding <- unlist(lapply(sv_coding,paste,collapse=" "))
samples_noncod <- unlist(lapply(sv_noncod,paste,collapse=" "))

topGenesTab <- cbind(genes_tab[allgenes,c("seqnames","start","end","strand")],genes_exonic[allgenes],genes_copynum[allgenes],genes_coding[allgenes],genes_proximal[allgenes],genes_intronic[allgenes],genes_noncod[allgenes],genes_summary[allgenes],allType[allgenes,],samples_coding[allgenes],samples_noncod[allgenes])
names(topGenesTab) <- c("seqnames","start","end","strand","exonic","copynumm","coding","proximal","intronic","non-coding","all",colnames(allType),"samples_coding","samples_noncod")
#write.table(topGenesTab,file="Supplemental_Table_S6_2.txt",sep="\t",quote=F)

# obtaine list of genes which have overlapping varinats and coverage based breakpoints
orth <- both <- sv <- bp <- all <- list()
for(i in names(sv_summary)){
	message(i)
	orth[[i]] <- unique(c(intersect(bp_noncod[[i]],sv_proximal[[i]]),intersect(bp_coding[[i]],c(sv_exonic[[i]],sv_intronic[[i]])),intersect(c(ampl_WGS[[i]],ddel_WGS[[i]]),sv_copynum[[i]])))
	both[[i]] <-intersect(c(bp_summary[[i]],ampl_WGS[[i]],ddel_WGS[[i]]),sv_summary[[i]])
	sv[[i]] <- unique(c(sv_proximal[[i]],sv_exonic[[i]],sv_intronic[[i]],sv_copynum[[i]]))
	bp[[i]] <- unique(c(bp_noncod[[i]],bp_coding[[i]],ampl_WGS[[i]],ddel_WGS[[i]]))
	all[[i]] <- unique(c(sv_proximal[[i]],sv_exonic[[i]],sv_intronic[[i]],sv_copynum[[i]],bp_noncod[[i]],bp_coding[[i]],ampl_WGS[[i]],ddel_WGS[[i]]))
}

#plotgenes <- setdiff(names(which(sort(unlist(lapply(orth,length)),decreasing=F) > 2)),passenger)
#cosmic <-  as.character(cosmic_genes$Gene.Symbol)

orthogonal <- names(which(sort(unlist(lapply(orth,length)),decreasing=F) > 0))
sv_bp_combined <- unlist(lapply(merge2lists(merge2lists(bp_summary,sv_summary),merge2lists(ampl_WGS,ddel_WGS)),length))

all_recurrent_noncod <- names(which(unlist(lapply(merge2lists(bp_noncod,sv_noncod),length)) >5))
all_recurrent_coding <- names(which(unlist(lapply(merge2lists(bp_coding,sv_coding),length)) >2))
shortlist <- intersect(names(which(sv_bp_combined > 3)),unique(c(intersect(all_recurrent_noncod,orthogonal),intersect(orthogonal,all_recurrent_coding))))
shortlist <- sort(grep("LOC|orf|-",shortlist,invert=T,value=T))

## cis-eQTL analysis
commongenes <- intersect(rownames(l_expmat),rownames(r_expmat))
a<-zrank(l_expmat[commongenes,intersect(highrisk,setdiff(colnames(l_expmat),colnames(r_expmat)))])
b<-zrank(r_expmat[commongenes,intersect(highrisk,colnames(r_expmat))])
lr_expmat <- t(apply(cbind(a,b),1,rank))

## supplementary Figure S18A 
res1 <- sv_eQTL(shortlist,merge2lists(sv_summary,bp_summary),r_expmat,highrisk) 
#write.table(do.call(rbind,res1),sep="\t",quote=F)
res2 <- sv_eQTL(shortlist,merge2lists(sv_summary,bp_summary),l_expmat,highrisk) 
#write.table(do.call(rbind,res2),sep="\t",quote=F)
res3 <- sv_eQTL(shortlist,merge2lists(sv_summary,bp_summary),lr_expmat,highrisk) 
#write.table(do.call(rbind,res3),sep="\t",quote=F)

resa <- data.frame(do.call(rbind,res1))
resb <- data.frame(do.call(rbind,res2))
resc <- data.frame(do.call(rbind,res3))
for(i in 1:3){
	resa[,i] <- as.numeric(as.character(resa[,i]))
	resb[,i] <- as.numeric(as.character(resb[,i]))
	resc[,i] <- as.numeric(as.character(resc[,i]))
}

logpa <- -log10(resa$wcox)*sign(log(resa$FC))
logpb <- -log10(resb$wcox)*sign(log(resb$FC))
logpc <- -log10(resc$wcox)*sign(log(resc$FC))


bardata <- rbind(logpa,logpb,logpc)
colnames(bardata) <- rownames(resa)
bardata<-bardata[,order(apply(bardata,2,mean,na.rm=T),decreasing=T)]
barplot(bardata,beside=T,las=2)
abline(h=c(-2,-1.30103,1.30103,2))

heatmap.3(bardata[3:1,],col=bluered(256),Colv=F,Rowv=F)
#write.table(data.frame(genes_tab[rownames(resa),],resa,logpa,resb,logpb,resc,logpc)[colnames(bardata),],file="",sep="\t",quote=F)


orth_rank <- unlist(lapply(orth,length))[shortlist]
snv_rank <- unlist(lapply(snv_list,length))[shortlist]
both_rank <- unlist(lapply(both,length))[shortlist]
sv_rank <- unlist(lapply(sv,length))[shortlist]
bp_rank <- unlist(lapply(bp,length))[shortlist]
all_rank <- unlist(lapply(all,length))[shortlist]
freqtab <- rbind(orth_rank, both_rank-orth_rank,sv_rank-both_rank,bp_rank-both_rank,all_rank)

# Supplementary table S6
supplTabS6 <- data.frame(orth_rank[shortlist],topGenesTab[shortlist,],bpgene_tab[shortlist,],sv_bp_combined[shortlist])
#write.table(supplTabS6,file="",sep="\t",quote=F)
shortUniqList <- rownames(supplTabS6)

# from shortlist we remove "passenger" genes
plotOrderGenes <- names(sort(orth_rank[shortUniqList[order(all_rank[shortUniqList],decreasing=T)]],decreasing=T))

### Incorporate SNV maf data 
snv_extended <- read.delim("data/extended_cgi_Aug_1_18.maf",as.is=T,sep="\t")
snv_list<-list()
for(i in plotOrderGenes){
  snv_list[[i]] <- unique(snv_extended$Sample[which(snv_extended$Hugo_Symbol == i )])
}

# plor stack barplot attached to oncoprint Figure 3G
par(mar=c(3,7,1,6))
bplot <- barplot(freqtab[1:4,rev(plotOrderGenes)],horiz=T,las=1,xlim=c(0,40),col=c("black","grey40","grey80","white"))
abline(v=seq(0,40,10),lty=3,lwd=.3)
legend("right",c("Orthogonal","Both","SJ-BP","RD-BP"),fill=c("black","grey40","grey80","white"),bty='n',cex=1.2)
#text(11+apply(freqtab[1:4,rev(plotOrderGenes)],2,sum),bplot,labels=paste(sprintf("%.1f",100*freqtab[4,rev(plotOrderGenes)]/135),"%",sep=""),cex=.9)

# ONCOPRINT Figure 3GFigure 3G
bp_coding2 <-bp_coding
bp_noncod2 <-bp_noncod
ampl_WGS2 <-ampl_WGS
ddel_WGS2 <-ddel_WGS
snv_list2 <- snv_list
names(bp_coding2) <-paste(names(bp_coding),"_2",sep="")
names(bp_noncod2) <-paste(names(bp_noncod),"_2",sep="")
names(ampl_WGS2)  <-paste(names(ampl_WGS),"_2",sep="")
names(ddel_WGS2) <- paste(names(ddel_WGS),"_2",sep="")
names(snv_list2) <- paste(names(snv_list2),"_3",sep="")


topgenes2x <- rep(NA,4*length(plotOrderGenes))
topgenes2x[seq(1,4*length(plotOrderGenes),4)] <- plotOrderGenes
topgenes2x[seq(2,4*length(plotOrderGenes),4)] <- paste(plotOrderGenes,"_2",sep="")
topgenes2x[seq(3,4*length(plotOrderGenes),4)] <- paste(plotOrderGenes,"_3",sep="")
topgenes2x[seq(4,4*length(plotOrderGenes),4)] <- paste(plotOrderGenes,"_4",sep="")

mat<-matrix(ncol=length(cgi_samples),nrow=length(topgenes2x))
colnames(mat) <- cgi_samples
rownames(mat) <- topgenes2x
for(i in topgenes2x){
	for(j in cgi_samples){
		mutlist <- list()
		if(j %in% sv_exonic[[i]] ){
			mutlist[["EXON"]] <- "EXON"
		}
		if(j %in% sv_copynum[[i]] ){
			mutlist[["CN"]] <- "CN"
		}
		if(j %in% sv_proximal[[i]] ){
			mutlist[["PROX"]] <- "PROX"
		}
		if(j %in% sv_intronic[[i]] ){
			mutlist[["INT"]] <- "INT"
		}
		if(j %in% ampl_WGS2[[i]] ){
			mutlist[["AMP"]] <- "AMP"
		}
		if(j %in% ddel_WGS2[[i]] ){
			mutlist[["DDEL"]] <- "DDEL"
		}
		if(j %in% snv_list2[[i]] ){
			mutlist[["MISS"]] <- "MISS"
		}
		if(j %in% bp_coding2[[i]] ){
			mutlist[["BPC"]] <- "BPC"
		}
		if(j %in% bp_noncod2[[i]] ){
			mutlist[["BPP"]] <- "BPP"
		}
		mat[i,j] <- paste(unlist(unname(mutlist)),sep="",collapse=";")
	}
}
mat[sort(c(seq(1,2*length(plotOrderGenes),4), seq(1,2*length(plotOrderGenes),4)+1)),] <- mat[sort(c(seq(1,2*length(plotOrderGenes),4), seq(1,2*length(plotOrderGenes),4)+1)),]
mat[sort(c(seq(3,2*length(plotOrderGenes),4), seq(3,2*length(plotOrderGenes),4)+1)),] <- mat[sort(c(seq(3,2*length(plotOrderGenes),4), seq(3,2*length(plotOrderGenes),4)+1)),]

mat2 <- mat
mat2[grep("_4",rownames(mat2)),] <- "white"


alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "grey95", col = NA))
    },
    white= function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "white", col = NA))
    },
    EXON = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "black", col = NA))
    },
    PROX = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "green4", col = NA))
    },
    CN= function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h*0.6, gp = gpar(fill = "pink", col = NA))
    },
    INT = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h*0.5, gp = gpar(fill = "gold", col = NA))
    },
    MISS= function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "black", col = NA))
    },
    BPC = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "black", col = NA))
    },
    BPP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "green4", col = NA))
    },
    AMP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h*0.6, gp = gpar(fill = "red", col = NA))
    },
#    DEL = function(x, y, w, h) {
#        grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.4, "mm"), gp = gpar(fill = "lightblue", col = NA))
#    },
    DDEL = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.2, "mm"), h*0.6, gp = gpar(fill = "cyan", col = NA))
    }
)
#col = c(EXON = "black", PROX = "green",INT="gold",CN="pink",GAIN="pink",AMP="red",DEL="lightblue",DDEL="blue",BPC="black",BPP="green",BG1="grey98", BG2="grey92")
col = c(EXON = "black", PROX = "green",INT="gold",CN="pink",AMP="red",DDEL="blue",BPC="black",BPP="green",MISS="black",white="white")

a<-oncoPrint(mat2[topgenes2x,c(HR_MNA,HR_NA,INT,LOW)], get_type = function(x) strsplit(x, ";")[[1]], row_order=NULL,
    alter_fun = alter_fun, col = col, 
    column_title = "SV analysis", 
    row_names_gp = gpar(fontsize = 10),
    pct_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Alternations", at = c("EXON","CN", "PROX","INT"), 
    labels = c("Exonic SV","Copy Number SV", "Proximal SV","Intronic SV" )))
a

allorder <-1:length(cgi_samples)
names(allorder) <- c(HR_MNA,HR_NA,INT,LOW)[a@column_order]


# plot MNA subtype oncoprint
b<-oncoPrint(mat2[topgenes2x,names(sort(allorder[HR_MNA]))], get_type = function(x) strsplit(x, ";")[[1]], row_order=NULL,
    alter_fun = alter_fun, col = col, 
    column_title = "SV analysis", column_order=NULL,
    row_names_gp = gpar(fontsize = 10),
    pct_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Alternations", at = c("EXON","CN", "PROX","INT"), 
    labels = c("Exonic SV","Copy Number SV", "Proximal SV","Intronic SV" )))
b


# plot HRNA subtype oncoprint
c<-oncoPrint(mat2[topgenes2x,names(sort(allorder[HR_NA]))], get_type = function(x) strsplit(x, ";")[[1]], row_order=NULL,
    alter_fun = alter_fun, col = col, 
    column_title = "SV analysis", column_order=NULL,
    row_names_gp = gpar(fontsize = 10),
    pct_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Alternations", at = c("EXON","CN", "PROX","INT"), 
    labels = c("Exonic SV","Copy Number SV", "Proximal SV","Intronic SV" )))
c

d<-oncoPrint(mat2[topgenes2x,names(sort(allorder[INT]))], get_type = function(x) strsplit(x, ";")[[1]], row_order=NULL,
    alter_fun = alter_fun, col = col, 
    column_title = "SV analysis", column_order=NULL,
    row_names_gp = gpar(fontsize = 10),
    pct_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Alternations", at = c("EXON","CN", "PROX","INT"), 
    labels = c("Exonic SV","Copy Number SV", "Proximal SV","Intronic SV" )))

# plot INT subtype oncoprint
d

e<-oncoPrint(mat2[topgenes2x,names(sort(allorder[LOW]))], get_type = function(x) strsplit(x, ";")[[1]], row_order=NULL,
    alter_fun = alter_fun, col = col, 
    column_title = "SV analysis", column_order=NULL,
    row_names_gp = gpar(fontsize = 10),
    pct_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Alternations", at = c("EXON","CN", "PROX","INT"), 
    labels = c("Exonic SV","Copy Number SV", "Proximal SV","Intronic SV" )))
# plot LOW subtype oncoprint
e


# generate Phenotype top bars for each subtype
fakemat <- matrix(ncol=length(HR_MNA),nrow=5)
fakemat[] <- 0
colnames(fakemat) <- names(sort(allorder[HR_MNA]))
rownames(fakemat) <- rownames(Pheno2)[2:6]
heatmap.3(fakemat, ColSideColors = t(Pheno2[2:6,colnames(fakemat)]),Colv=FALSE,Rowv=FALSE,lhei=c(0.5,1))

fakemat <- matrix(ncol=length(HR_NA),nrow=5)
fakemat[] <- 0
colnames(fakemat) <- names(sort(allorder[HR_NA]))
rownames(fakemat) <- rownames(Pheno2)[2:6]
heatmap.3(fakemat, ColSideColors = t(Pheno2[2:6,colnames(fakemat)]),Colv=FALSE,Rowv=FALSE,lhei=c(0.5,1))

fakemat <- matrix(ncol=length(INT),nrow=5)
fakemat[] <- 0
colnames(fakemat) <- names(sort(allorder[INT]))
rownames(fakemat) <- rownames(Pheno2)[2:6]
heatmap.3(fakemat, ColSideColors = t(Pheno2[2:6,colnames(fakemat)]),Colv=FALSE,Rowv=FALSE,lhei=c(0.5,1))

fakemat <- matrix(ncol=length(LOW),nrow=5)
fakemat[] <- 0
colnames(fakemat) <- names(sort(allorder[LOW]))
rownames(fakemat) <- rownames(Pheno2)[2:6]
heatmap.3(fakemat, ColSideColors = t(Pheno2[2:6,colnames(fakemat)]),Colv=FALSE,Rowv=FALSE,lhei=c(0.5,1))


## obtain altered genes from segmentation data
# SNP

sv_samples <- merge2lists(sv_coding,sv_noncod)
bp_coding <- results_BP$breakSamples
bp_proximal <- results_BP$proximalSamples
bp_samples <- merge2lists(bp_coding ,bp_proximal)

common_wgs_snp <- intersect(snp_samples,cgi_samples)
bp_proximal_snp <- results_BP_snp$proximalSamples
bp_proximal_snp_only <- lapply(bp_proximal_snp,function(x) setdiff(x,common_wgs_snp))
bp_proximal_snp_only <- bp_proximal_snp_only[which(lapply(bp_proximal_snp_only,length) > 0)]

bp_coding_snp <- results_BP_snp$breakSamples
bp_coding_snp_only <- lapply(bp_coding_snp,function(x) setdiff(x,common_wgs_snp))
bp_coding_snp_only <- bp_coding_snp_only[which(lapply(bp_coding_snp_only,length) > 0)]

bp_samples_snp_only <- merge2lists(bp_coding_snp_only,bp_proximal_snp_only)

all_samples <- merge2lists(sv_samples,bp_coding_snp)

a<-length(all_samples$SHANK2)
b<-length(all_samples$DLG2)
tt <-length(unique(c(cgi_samples,snp_samples)))

a<-length(intersect(all_samples$SHANK2,c(HR_NA_snp, HR_NA)))
b<-length(intersect(all_samples$DLG2,c(HR_NA_snp, HR_NA)))
tt <- length(c(HR_NA_snp, HR_NA))

# SHANK2 and DLG2 altered samples are mutually exclusive
fisher.test(matrix(c(0,a,b,tt-a-b),2,2),alternative="less")


a <- sort(unlist(lapply(bp_samples_snp,length)),decreasing=T)
tab <- cbind(genes_tab[names(a),c(6,1,2,3)],a,unlist(lapply(bp_coding_snp,length))[names(a)],unlist(lapply(bp_proximal_snp,length))[names(a)])
tab[which(is.na(tab),arr.ind=T)] <- 0
colnames(tab) <- c( "Symbol" ,"Chr","Start","End","Combined","Coding","Proximal")
#write.table(tab,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/BP_SNP_altered_genes_Oct5_18.txt",sep="\t",quote=F,row.names=F)

a <- sort(unlist(lapply(bp_samples,length)),decreasing=T)
tab <- cbind(genes_tab[names(a),c(6,1,2,3)],a,unlist(lapply(bp_coding,length))[names(a)],unlist(lapply(bp_proximal,length))[names(a)])
tab[which(is.na(tab),arr.ind=T)] <- 0
colnames(tab) <- c( "Symbol" ,"Chr","Start","End","Combined","Coding","Proximal")
#write.table(tab,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/BP_CGI_altered_genes_Oct5_18.txt",sep="\t",quote=F,row.names=F)

# regional gene sets associated with MYCN, TERT and other regions
chr2p24 <- c("DDX1","FAM49A","FAM84A","GACAT3","KCNS3","LINC00276","LINC01804","LPIN1","MYCN","MYCNOS","MYCNUT","NBAS","NT5C1B","NT5C1B-RDH14","RDH14","SMC6","VSNL1")
chr5p15 <- c("CLPTM1L","CTD-3080P12.3","LINC01511","LOC100506688","LPCAT1","SLC12A7","SLC6A18","SLC6A19","SLC6A3","TERT")
chr17p21 <- c("LINC00854","DHX8","LINC00910","ARL4D")
chr11q1314 <- c("PACS1", "PPFIA1", "SHANK2", "NUMA1", "GAB2", "DLG2")

# Figure 3A
par(mfrow=c(2,3),mar=c(3,4,2,1))
b<-sort(unlist(lapply(sv_coding,length)),decreasing=T)
plot(b,cex=.2,col="white",main="SJ-BP coding",xlab="",ylab="",las=1,cex.axis=1.2)
points(which(!names(b) %in% c(chr2p24,chr5p15)),b[which(!names(b) %in% c(chr2p24,chr5p15))],cex=.4)
points(which(names(b) %in% chr2p24),b[which(names(b) %in% chr2p24)],col="red",pch=3)
points(which(names(b) %in% chr5p15),b[which(names(b) %in% chr5p15)],col="blue",cex=1,pch=4)
#text(b,labels=names(b),pos=4)

# Figure 3B
a<-sort(unlist(lapply(bp_coding,length)),decreasing=T)
plot(a,cex=.2,col="white",main="RD-BP coding (WGS)",xlab="",ylab="",las=1,cex.axis=1.2)
points(which(!names(a) %in% c(chr2p24,chr5p15)),a[which(!names(a) %in% c(chr2p24,chr5p15))],cex=.4)
points(which(names(a) %in% chr2p24),a[which(names(a) %in% chr2p24)],col="red",pch=3)
points(which(names(a) %in% chr5p15),a[which(names(a) %in% chr5p15)],col="blue",cex=1,pch=4)
#text(a,labels=names(a),pos=4)

# Figure 3C
b<-sort(unlist(lapply(bp_coding_snp_only,length)),decreasing=T)
plot(b,cex=.2,col="white",main="CN-BP coding (SNP)",xlab=" ",ylab="",las=1,cex.axis=1.2)
points(which(!names(b) %in% c(chr2p24,chr5p15)),b[which(!names(b) %in% c(chr2p24,chr5p15))],cex=.4)
points(which(names(b) %in% chr2p24),b[which(names(b) %in% chr2p24)],col="red",pch=3)
points(which(names(b) %in% chr5p15),b[which(names(b) %in% chr5p15)],col="blue",cex=1,pch=4)
#text(b,labels=names(b),pos=4)

# Figure 3D
c<-sort(unlist(lapply(sv_noncod,length)),decreasing=T)
plot(c,cex=.2,col="white",main="SJ-BP non-coding",xlab=" ",ylab="",las=1,cex.axis=1.2)
points(which(!names(c) %in% c(chr2p24,chr5p15)),c[which(!names(c) %in% c(chr2p24,chr5p15))],cex=.4)
points(which(names(c) %in% chr2p24),c[which(names(c) %in% chr2p24)],col="red",pch=3)
points(which(names(c) %in% chr5p15),c[which(names(c) %in% chr5p15)],col="blue",cex=1,pch=4)
#text(c,labels=names(c),pos=4)

# Figure 3E
c<-sort(unlist(lapply(bp_noncod,length)),decreasing=T)
plot(c,cex=.2,col="white",main="RD-BP non-coding (WGS)",xlab="",ylab="",las=1,cex.axis=1.2)
points(which(!names(c) %in% c(chr2p24,chr5p15)),c[which(!names(c) %in% c(chr2p24,chr5p15))],cex=.4)
points(which(names(c) %in% chr2p24),c[which(names(c) %in% chr2p24)],col="red",pch=3)
points(which(names(c) %in% chr5p15),c[which(names(c) %in% chr5p15)],col="blue",cex=1,pch=4)
#text(c,labels=names(c),pos=4)

# Figure 3F
c<-sort(unlist(lapply(bp_proximal_snp_only,length)),decreasing=T)
plot(c,cex=.2,col="white",main="CN-BP non-coding (SNP)",xlab="",ylab="",las=1,cex.axis=1.2)
points(which(!names(c) %in% c(chr2p24,chr5p15)),c[which(!names(c) %in% c(chr2p24,chr5p15))],cex=.4)
points(which(names(c) %in% chr2p24),c[which(names(c) %in% chr2p24)],col="red",pch=3)
points(which(names(c) %in% chr5p15),c[which(names(c) %in% chr5p15)],col="blue",cex=1,pch=4)
#text(c,labels=names(c),pos=4)

plot(1,1,col="white")
legend("center",c("2p24 locus","5p15 locus","other loci"),ncol=4,pch=c(3,4,0,21),pt.cex=c(1,1,1,.4),col=c("red","blue","black"))

### PATHWAY ENRICHMENT ANALYSIS , Figure A-C & supplementary S19

a<-sort(unlist(lapply(sv_coding,length)),decreasing=T)
b<-sort(unlist(lapply(bp_coding,length)),decreasing=T)
c<-sort(unlist(lapply(bp_coding_snp_only,length)),decreasing=T)
  
d<-sort(unlist(lapply(sv_noncod,length)),decreasing=T)
e<-sort(unlist(lapply(bp_noncod,length)),decreasing=T)
f<-sort(unlist(lapply(bp_proximal_snp_only,length)),decreasing=T)
# supplementary table S5
write.table(names(which(a>2)),row.names=FALSE,quote=F,col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_sv_coding_nov6.txt")
write.table(names(which(b>2)),row.names=FALSE,quote=F,col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_bp_coding_nov6.txt")
write.table(names(which(c>2)),row.names=FALSE,quote=F,col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_bp_codingSnp_nov6.txt")
write.table(names(which(d>3)),row.names=FALSE,quote=F,col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_sv_noncod_nov6.txt")
write.table(names(which(e>3)),row.names=FALSE,quote=F,col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_bp_noncod_nov6.txt")
write.table(names(which(f>3)),row.names=FALSE,quote=F,col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_bp_noncodSnp_nov6.txt")

write.table(cbind(a,a/135),row.names=TRUE,quote=F,sep="\t",col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/table4/table_4_rank_sv_coding_nov6.txt")
write.table(cbind(b,b/135),row.names=TRUE,quote=F,sep="\t",col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/table4/table_4_rank_bp_coding_nov6.txt")
write.table(cbind(c,c/915),row.names=TRUE,quote=F,sep="\t",col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/table4/table_4_rank_bp_codingSnp_nov6.txt")
write.table(cbind(d,d/135),row.names=TRUE,quote=F,sep="\t",col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/table4/table_4_rank_sv_noncod_nov6.txt")
write.table(cbind(e,e/135),row.names=TRUE,quote=F,sep="\t",col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/table4/table_4_rank_bp_noncod_nov6.txt")
write.table(cbind(f,f/915),row.names=TRUE,quote=F,sep="\t",col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/table4/table_4_rank_bp_noncodSnp_nov6.txt")
# other tumors...


topgenecats <- c("GO: Biological Process","GO: Cellular Component","GO: Molecular Function","Pathway","Disease")
davidcats <- c("GOTERM_CC_DIRECT","GOTERM_BP_DIRECT","BIOCARTA","GOTERM_MF_DIRECT","KEGG_PATHWAY")

par(mar=c(3,30,1,3),mfrow=c(3,1))
## coding SV
data <- read.delim("data/FUNCTION/FUNCTION_topgenes_sv_coding_nov6.txt",sep="\t",as.is=T)
data <- data[which(data$Category %in% topgenecats),]
data <- data[rev(order(data$q.value.Bonferroni)[1:20]),]
categoriesnames<-gsub(": Cellular Component","_CC",data$Category);categoriesnames<-gsub(": Molecular Function","_MF",categoriesnames);categoriesnames<-gsub(": Biological Process","_BP",categoriesnames);
datanames1 <- apply(cbind(categoriesnames,data$Name),1,paste,collapse=": ") 
barplot(rbind(-log10(data$q.value.Bonferroni), -log10(data$p.value) +log10(data$q.value.Bonferroni) ),names=datanames1,horiz=T,las=1)
abline(v=seq(2,8,2),lty=3,lwd=0.5)
## coding BP
data <- read.delim("data/FUNCTION/FUNCTION_topgenes_bp_coding_nov6.txt",sep="\t",as.is=T)
data <- data[which(data$Category %in% topgenecats),]
data <- data[rev(order(data$q.value.Bonferroni)[1:20]),]
categoriesnames<-gsub(": Cellular Component","_CC",data$Category);categoriesnames<-gsub(": Molecular Function","_MF",categoriesnames);categoriesnames<-gsub(": Biological Process","_BP",categoriesnames);
datanames2 <- apply(cbind(categoriesnames,data$Name),1,paste,collapse=": ") 
barplot(rbind(-log10(data$q.value.Bonferroni), -log10(data$p.value) +log10(data$q.value.Bonferroni) ),names=datanames2,horiz=T,las=1)
abline(v=seq(2,8,2),lty=3,lwd=0.5)
## coding BP SNP
data <- read.delim("data/FUNCTION/FUNCTION_topgenes_bp_codingSnp_nov6.txt",sep="\t",as.is=T)
data <- data[which(data$Category %in% topgenecats),]
data <- data[rev(order(data$q.value.Bonferroni)[1:20]),]
categoriesnames<-gsub(": Cellular Component","_CC",data$Category);categoriesnames<-gsub(": Molecular Function","_MF",categoriesnames);categoriesnames<-gsub(": Biological Process","_BP",categoriesnames);
datanames3 <- apply(cbind(categoriesnames,data$Name),1,paste,collapse=": ") 
barplot(rbind(-log10(data$q.value.Bonferroni), -log10(data$p.value) +log10(data$q.value.Bonferroni) ),names=datanames3,horiz=T,las=1)
abline(v=seq(2,8,2),lty=3,lwd=0.5)

write.table(c(rev(datanames1),
rev(datanames2),
rev(datanames3)),row.names=F)


par(mar=c(3,40,1,3),mfrow=c(3,1))

## non-coding SV
data <- read.delim("data/FUNCTION/FUNCTION_topgenes_sv_noncod_nov6.txt",sep="\t",as.is=T)
data <- data[which(data$Category %in% topgenecats),]
data <- data[rev(order(data$q.value.Bonferroni)[1:20]),]
categoriesnames<-gsub(": Cellular Component","_CC",data$Category);categoriesnames<-gsub(": Molecular Function","_MF",categoriesnames);categoriesnames<-gsub(": Biological Process","_BP",categoriesnames);
datanames1 <- apply(cbind(categoriesnames,data$Name),1,paste,collapse=": ") 
barplot(rbind(-log10(data$q.value.Bonferroni), -log10(data$p.value) +log10(data$q.value.Bonferroni) ),names=datanames1,horiz=T,las=1)
abline(v=seq(2,10,2),lty=3,lwd=0.5)
## non coding BP
data <- read.delim("data/FUNCTION/FUNCTION_topgenes_bp_noncod_nov6.txt",sep="\t",as.is=T)
data <- data[which(data$Category %in% topgenecats),]
data <- data[rev(order(data$q.value.Bonferroni)[1:20]),]
categoriesnames<-gsub(": Cellular Component","_CC",data$Category);categoriesnames<-gsub(": Molecular Function","_MF",categoriesnames);categoriesnames<-gsub(": Biological Process","_BP",categoriesnames);
datanames2 <- apply(cbind(categoriesnames,data$Name),1,paste,collapse=": ") 
barplot(rbind(-log10(data$q.value.Bonferroni), -log10(data$p.value) +log10(data$q.value.Bonferroni) ),names=datanames2,horiz=T,las=1)
abline(v=seq(2,10,2),lty=3,lwd=0.5)
## non coding SNP
data <- read.delim("data/FUNCTION/FUNCTION_topgenes_bp_noncodSnp_nov6.txt",sep="\t",as.is=T)
data <- data[which(data$Category %in% topgenecats),]
data <- data[rev(order(data$q.value.Bonferroni)[1:20]),]
categoriesnames<-gsub(": Cellular Component","_CC",data$Category);categoriesnames<-gsub(": Molecular Function","_MF",categoriesnames);categoriesnames<-gsub(": Biological Process","_BP",categoriesnames);
datanames3 <- apply(cbind(categoriesnames,data$Name),1,paste,collapse=": ") 
barplot(rbind(-log10(data$q.value.Bonferroni), -log10(data$p.value) +log10(data$q.value.Bonferroni) ),names=datanames3,horiz=T,las=1)
abline(v=seq(2,10,2),lty=3,lwd=0.5)



## join exonic and CNn into coding and intronic and proximal into noncoding

all_coding <- merge2lists( lapply(results_ALL$disruptExonSamples,function(x) substr(x,11,16)),lapply(results_ALL$copynumSamples,function(x) substr(x,11,16)))
aml_coding <- merge2lists( lapply(results_AML$disruptExonSamples,function(x) substr(x,11,16)),lapply(results_AML$copynumSamples,function(x) substr(x,11,16)))
os_coding <- merge2lists( lapply(results_OS$disruptExonSamples,function(x) substr(x,11,16)),lapply(results_OS$copynumSamples,function(x) substr(x,11,16)))
wt_coding <- merge2lists( lapply(results_WT$disruptExonSamples,function(x) substr(x,11,16)),lapply(results_WT$copynumSamples,function(x) substr(x,11,16)))



all4tumors <- merge2lists(merge2lists(all_coding,aml_coding),merge2lists(os_coding,wt_coding))

write.table(names(which(unlist(lapply(all_coding,length)) > 3)),sep="\t",quote=F,row.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_sv_coding_ALL.txt")
write.table(names(which(unlist(lapply(aml_coding,length)) > 2)),sep="\t",quote=F,row.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_sv_coding_AML.txt")
write.table(names(which(unlist(lapply(os_coding,length)) > 3)),sep="\t",quote=F,row.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_sv_coding_OS.txt")
write.table(names(which(unlist(lapply(wt_coding,length)) > 2)),sep="\t",quote=F,row.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_sv_coding_WT.txt")

write.table(names(which(unlist(lapply(all4tumors,length)) > 5)),sep="\t",quote=F,row.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_sv_coding_all4tumors.txt")
xx<-sort(unlist(lapply(all4tumors,length)),decreasing=T)
nsamples4tum <- c(as.character(unique(results_ALL$sv_df$TARGET.USI)),as.character(unique(results_AML$sv_df$TARGET.USI)),as.character(unique(results_OS$sv_df$TARGET.USI)),as.character(unique(results_WT$sv_df$TARGET.USI)))
write.table(cbind(xx,xx/length(nsamples4tum)),row.names=TRUE,quote=F,sep="\t",col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/table4/table_4_rank_sv_coding_all4tumors.txt")

g <- sort(unlist(lapply(bp_coding_4tum,length)),decreasing=T)
g2 <- g[intersect(names(g),names(all4tumors))]
write.table(names(which(g2 > 5)),quote=F,row.names=F)

write.table(names(which(g2 > 5)),sep="\t",quote=F,row.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/topgenes_bp_coding_all4tumors.txt")
write.table(cbind(g2 ,g2 /length(nsamples4tum)),row.names=TRUE,quote=F,sep="\t",col.names=F,file="~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/table4/table_4_rank_bp_coding_all4tumors.txt")

out <-rankset.enrich(rank(g2[which(g2 > 3)]),msigDBsymbol$c2.all.v5.2)
head(out,20)

dev.new()
par(mar=c(3,30,1,3))
data <- read.delim("~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/FUNCTION_topgenes_sv_coding_all4tumors.txt",sep="\t",as.is=T)
data <- data[which(data$Category %in% topgenecats),]
data <- data[rev(order(data$q.value.Bonferroni)[1:60]),]
categoriesnames<-gsub(": Cellular Component","_CC",data$Category);categoriesnames<-gsub(": Molecular Function","_MF",categoriesnames);categoriesnames<-gsub(": Biological Process","_BP",categoriesnames);
datanames1 <- apply(cbind(categoriesnames,data$Name),1,paste,collapse=": ") 
barplot(rbind(-log10(data$q.value.Bonferroni), -log10(data$p.value) +log10(data$q.value.Bonferroni) ),names=datanames1,horiz=T,las=1,cex.names=.8,xaxt='n')
axis(1,at=seq(0,18,3))
abline(v=seq(0,18,3),lty=3,lwd=0.5)

dev.new()
par(mar=c(3,30,1,3))
data <- read.delim("~/Box Sync/My_CHOP/SV_paper_V2/data_sheets/Function/FUNCTION_topgenes_bp_coding_all4tumors.txt",sep="\t",as.is=T)
data <- data[which(data$Category %in% topgenecats),]
data <- data[rev(order(data$q.value.Bonferroni)[1:60]),]
categoriesnames<-gsub(": Cellular Component","_CC",data$Category);categoriesnames<-gsub(": Molecular Function","_MF",categoriesnames);categoriesnames<-gsub(": Biological Process","_BP",categoriesnames);
datanames1 <- apply(cbind(categoriesnames,data$Name),1,paste,collapse=": ") 
barplot(rbind(-log10(data$q.value.Bonferroni), -log10(data$p.value) +log10(data$q.value.Bonferroni) ),names=datanames1,horiz=T,las=1,cex.names=.8,xaxt='n')
axis(1,at=seq(0,18,3))
abline(v=seq(0,18,3),lty=3,lwd=0.5)

#

###### expression and survival analysis 
library(limma)
library(corto)
library(edgeR)
load("data/MSigDB_v6.1_human.rda",verbose=T)

synapse <- msigDBsymbol$msigdb.v6.1.symbols.gmt$GO_SYNAPSE
neurogenesys <-  msigDBsymbol$msigdb.v6.1.symbols.gmt$GO_NEUROGENESIS
neuronal<-  msigDBsymbol$msigdb.v6.1.symbols.gmt$GO_NEURON_PART
psd <- msigDBsymbol$msigdb.v6.1.symbols.gmt$GO_POSTSYNAPSE

# autism geneset from toppgene
autism <- read.delim("data/autism_disorder_genes.txt",as.is=T)[,2]
# COSMIC cancer gene census
cosmic_genes <- as.character(read.delim(paste("data/cancer_gene_census_v85.csv",sep=""),sep=",")$Gene.Symbol)


## humanExon array differential expression and GSEA
# Figure 5D-F and supplementary S20
expmat<-l_expmat

hr <- intersect(highrisk,colnames(l_expmat))
low <- intersect(intersect(names(which(risk == "low")),mycn_wt),colnames(l_expmat))

l_groups <- c(rep(1,length(hr)),rep(2,length(low)))
names(l_groups) <-  c(hr,low)

samples <- as.factor(c(l_groups)  )
design <- model.matrix(~ -1+samples)
colnames(design) <- c("hr","low")
contrast.matrix <- makeContrasts(hr-low,levels=design)
	# Fit a linear model to the data (limma)
l_fit <- lmFit(l_expmat[,names(l_groups)], design)
l_fit2 <- contrasts.fit(l_fit, contrast.matrix)
l_fit2 <- eBayes(l_fit2)
l_fdr <- apply(l_fit2$p.value,2,p.adjust,method="fdr")
colnames(l_fdr) <- c("hr_low")
l_fc <- cbind(apply(l_expmat[,hr],1,mean,na.rm=T) - apply(l_expmat[,low],1,mean,na.rm=T))
colnames(l_fc)<-c("hr_low")

l_neuro <-  intersect(intersect(c(neuronal,synapse),plotOrderGenes),rownames(l_expmat))
l_autism <- intersect(intersect(autism,plotOrderGenes),rownames(expmat))
l_cosmic_genes<- intersect(intersect(cosmic_genes,plotOrderGenes),rownames(expmat))
plot(l_fc[,1],-log10(l_fit2$p.value[,1]),pch=19,col=rgb(0, 0, 0, 0.1),xlim=c(-2,2),frame.plot=FALSE,las=1  )
points(l_fc[l_autism,1],-log10(l_fit2$p.value[,1])[l_autism],col="red",cex=.4,pch=19)
points(l_fc[l_neuro,1],-log10(l_fit2$p.value[,1])[l_neuro],col="green",cex=1,pch=21,lwd=1.2)
points(l_fc[l_cosmic_genes,1],-log10(l_fit2$p.value[,1])[l_cosmic_genes],col="purple",cex=.7,pch=21,lwd=1.2)
text(l_fc[l_autism,1],-log10(l_fit2$p.value[,1])[l_autism],labels=l_autism,col="blue")
text(l_fc[l_neuro,1],-log10(l_fit2$p.value[,1])[l_neuro],labels=l_neuro,col="green")
text(l_fc[l_cosmic_genes,1],-log10(l_fit2$p.value[,1])[l_cosmic_genes],labels=l_cosmic_genes,col="gold")
legend("topright",c("COSMIC census","Autism disorders","Neuron/Synapse","other"),pch=c(21,19,21,19),pt.lwd=c(1.2,1,1.2,1), 
	col=c("purple","red","green","#0000001A"),bty='n',pt.cex=c(0.7,0.4,1,1) )

ranklist <- -log10(l_fit2$p.value[,1])*sign(l_fc[,1])
g1 <- gsea(ranklist,autism,method="pareto")
plot_gsea(g1)

ranklist <- -log10(l_fit2$p.value[,1])*sign(l_fc[,1])
g1 <- gsea(ranklist,unique(c(neuronal,synapse)),method="pareto")
plot_gsea(g1)



## limma analysis for RNA-seq

hr <- intersect(highrisk,colnames(r_expmat))
low <- intersect(intersect(names(which(risk == "low")),mycn_wt),colnames(r_expmat))
expressed <- intersect(names(which(apply(r_expmat[,hr],1,mean,na.rm=T) > 0.1)),names(which(apply(r_expmat[,low],1,mean,na.rm=T) > 0.1)))

D <- calcNormFactors(r_expmat[expressed,])
y <- voom(r_expmat[expressed,],plot=TRUE,lib.size=colSums(r_expmat[expressed,])* D)
expmat <- y$E


r_groups <- c(rep(1,length(hr)),rep(2,length(low)))
names(r_groups) <-  c(hr,low)

samples <- as.factor(c(r_groups)  )
design <- model.matrix(~ -1+samples)
colnames(design) <- c("hr","low")
contrast.matrix <- makeContrasts(hr-low,levels=design)
	# Fit a linear model to the data (limma)
r_fit <- lmFit(expmat[,names(r_groups)], design)
r_fit2 <- contrasts.fit(r_fit, contrast.matrix)
r_fit2 <- eBayes(r_fit2)
r_fdr <- apply(r_fit2$p.value,2,p.adjust,method="fdr")
colnames(r_fdr) <- c("hr_low")
r_fc <- cbind(log(apply(r_expmat[expressed,hr],1,mean,na.rm=T)/apply(r_expmat[expressed,low],1,mean,na.rm=T)))
colnames(r_fc)<-c("hr_low")

l_neuro <-  intersect(intersect(c(neuronal,synapse),plotOrderGenes),expressed)
l_autism <- intersect(intersect(autism,plotOrderGenes),expressed)
l_cosmic_genes<- intersect(intersect(cosmic_genes,plotOrderGenes),expressed)
plot(r_fc[,1],-log10(r_fit2$p.value[,1]),pch=19,col=rgb(0, 0, 0, 0.1),frame.plot=FALSE,las=1 )
points(r_fc[l_autism,1],-log10(r_fit2$p.value[,1])[l_autism],col="red",cex=.4,pch=19)
points(r_fc[l_neuro,1],-log10(r_fit2$p.value[,1])[l_neuro],col="green",cex=1,pch=21,lwd=1.2)
points(r_fc[l_cosmic_genes,1],-log10(r_fit2$p.value[,1])[l_cosmic_genes],col="purple",cex=.7,pch=21,lwd=1.2)
text(r_fc[l_autism,1],-log10(r_fit2$p.value[,1])[l_autism],labels=l_autism,col="blue")
text(r_fc[l_neuro,1],-log10(r_fit2$p.value[,1])[l_neuro],labels=l_neuro,col="green")
text(r_fc[l_cosmic_genes,1],-log10(r_fit2$p.value[,1])[l_cosmic_genes],labels=l_cosmic_genes,col="gold")
legend("topright",c("COSMIC genes","Autism disorders","Neuron/Synapse","other"),pch=c(21,19,21,19),pt.lwd=c(1.2,1,1.2,1), 
	col=c("purple","red","green","#0000001A"),bty='n',pt.cex=c(0.7,0.4,1,1) )

ranklist <- -log10(r_fit2$p.value[,1])*sign(r_fc[,1])
g1 <- gsea(ranklist,autism,method="pareto")
plot_gsea(g1)

ranklist <- -log10(r_fit2$p.value[,1])*sign(r_fc[,1])
g1 <- gsea(ranklist,unique(c(neuronal,synapse)),method="pareto")
plot_gsea(g1)



## eQTL analysis (supplementary Figure S18)

expmat<-r_expmat
expmat<-l_expmat

lowintrisk <- names(risk)[which(risk %in% c("low","intermediate"))]
hr <- intersect(highrisk,colnames(expmat))
mna <- intersect(hr,mycn_amp)
na <- intersect(hr,mycn_wt)
lir <- intersect(lowintrisk,colnames(expmat))
low <- intersect(intersect(names(which(risk == "low")),mycn_wt),lir)
int  <- intersect(intersect(names(which(risk == "intermediate")),mycn_wt),lir)


genelist <- shortlist
genelist<- plotOrderGenes
genelist<- intersect(plotOrderGenes,rownames(expmat))

samplelist<- intersect(c(HR_MNA,HR_NA),colnames(expmat))
sv_summary_exp <- lapply(sv_summary,function(x) intersect(x,samplelist))
sv_summary_exp <- sv_summary_exp[intersect(names(sv_summary_exp),rownames(expmat))]
a<-sv_eQTL(genelist,sv_summary_exp,expmat,samplelist)
do.call(rbind,a)


res_amp<-res_noa<-list()
for(g in genelist){
	wt <- wilcox.test(expmat[g,mna],expmat[g,low])
	res_amp[[g]] <-c(wt$p.value,wt$statistic)
	wt <- wilcox.test(expmat[g,na],expmat[g,low])
	res_noa[[g]] <-c(wt$p.value,wt$statistic)
}
a<-do.call(rbind,res_amp)
a<-a[order(a[,2]),]
b<-do.call(rbind,res_noa)
b<-b[order(b[,2]),]

alog <- -log10(a[,1])*sign(log(a[,2] /a[which(a[,1] == max(a[,1])),2]))
blog <- -log10(b[,1])*sign(log(b[,2] /b[which(b[,1] == max(b[,1])),2]))

x <- -log10(p.adjust(a[,1]))
x <- x[which(x > 0.5) ]
y<- -log10(a[,1])[names(x)]
plot(x,y,type='o')
abline(lm(y~x))
abline(v=1)
fdr_cut_a <-2.90418

x <- -log10(p.adjust(b[,1]))
x <- x[which(x > 0.5) ]
y<- -log10(b[,1])[names(x)]
plot(x,y)
abline(lm(y~x))
abline(v=1)
#locator()
fdr_cut_b <-2.997203


plot(alog,blog[names(alog)],col="white")
#text(alog,blog[names(alog)],labels=names(alog),cex=.6,pos=1)
rect(-fdr_cut_a,-fdr_cut_b,fdr_cut_a,fdr_cut_b)
rect(-100,-fdr_cut_b,100,-100,col=rgb(0, 0, 1, 0.1, names = NULL, maxColorValue = 1),border=NA )
rect(-100,fdr_cut_b,100,100,col=rgb(1, 0, 0, 0.1, names = NULL, maxColorValue = 1),border=NA )
rect(-100,-100,-fdr_cut_a,100,col=rgb(0, 0, 1, 0.1, names = NULL, maxColorValue = 1),border=NA )
rect(100,-100,fdr_cut_a,100,col=rgb(1, 0, 0, 0.1, names = NULL, maxColorValue = 1),border=NA )
#points(alog[synapse],blog[synapse],cex=.8,pch=19,col=rgb(0, 1, 1, 0.7, names = NULL, maxColorValue = 1))
points(alog,blog[names(alog)],pch=19,col=rgb(0, 0, 0, 0.1, names = NULL, maxColorValue = 1))
#points(alog[autism],blog[autism],cex=.8,pch=19,col=rgb(1, 0, 0, 0.7, names = NULL, maxColorValue = 1))
#text(alog[autism],blog[autism],labels=autism,col="red4",cex=.6,pos=1)
points(alog[autism],blog[autism],col="black",cex=.4,pch=19)
points(alog[unique(c(neuronal,synapse))],blog[unique(c(neuronal,synapse))],col="green",cex=1,pch=21,lwd=1.2)
points(alog[cosmic_genes_short],blog[cosmic_genes_short],col="purple",cex=.7,pch=21,lwd=1.2)
#text(alog[unique(c(neuronal,autism,synapse,"TERT"))],blog[unique(c(neuronal,autism,synapse,"TERT"))],labels=unique(c(neuronal,autism,synapse,"TERT")),col="black",cex=.7,pos=1)
#text(alog[neuronal],blog[neuronal],labels=neuronal,col="black",cex=.7,pos=1)
#text(alog[autism],blog[autism],labels=autism,col="black",cex=.7,pos=1)
legend("bottomright",c("COSMIC census","Autism disorders","Neuron/Synapse","other"),pch=c(21,19,21,19),pt.lwd=c(1.2,1,1.2,1), 
	col=c("purple","black","green","#0000001A"),bty='n',pt.cex=c(0.7,0.4,1,1) )
text(alog,blog[names(alog)],labels=names(alog),col="black",cex=.7,pos=1)




#### Km plots based on expressio Supplementary Figure S22 DE
source("R/my_survival.r")
seqc_survival <- read.delim("data/seqc_survival_clinical.txt",row.names=1,as.is=TRUE)
seqc_expmat <- fread("data/GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt.gz",header = TRUE,data.table = FALSE)
s_expmat <- as.matrix(seqc_expmat[2:ncol(seqc_expmat)])
rownames(s_expmat) <- seqc_expmat$RefSeqID
colnames(s_expmat) <- substr(colnames(s_expmat) ,6,12)

s_surv <- seqc_survival$OS_d
s_event <- seqc_survival$OS_bin
s_stage <- seqc_survival$INSS
s_risk <- seqc_survival$high_risk
s_mycna <- rep(NA,nrow(seqc_survival))
s_mycna[which(seqc_survival$MYCN == "nonMNA")] <- 0
s_mycna[which(seqc_survival$MYCN == "MNA")] <- 1
s_age <- seqc_survival$age
names(s_surv) <- names(s_event) <- names(s_risk)<-names(s_stage) <-names(s_mycna)<-names(s_age)<-rownames(seqc_survival)
s_expmat <- 2^s_expmat

# S22D
km_r2(s_expmat["NM_012309",],s_surv,s_event,varname=NULL,main="")
# S22E
km_r2(s_expmat["NM_012309",names(which(s_risk != "HR" ))],s_surv,s_event,varname=NULL,main="")



#  #  # END

