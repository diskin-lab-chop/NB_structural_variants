


readDepthCopynum <- function(segdat, feature_tab, mapped_genes, ampl_cut = 2, ddel_cut = -2){

feature_tab_GR = with(feature_tab, GRanges(seqnames, IRanges(start=start, end=end),Rle(strand)))
mcols(feature_tab_GR)$name2 <-rownames(feature_tab)

genecn <- list()
for(id in  unique(segdat$Sample) ){
	message(id)
	segdat_i<-	segdat[which(segdat$Sample == id),]
	segdat_i_GR = with(segdat[which(segdat$Sample == id),], GRanges(Chromosome, IRanges(start=Start, end=End)))

	overlap_gr <- GenomicAlignments::findOverlaps(feature_tab_GR,segdat_i_GR,ignore.strand=TRUE,type="within")
	genenames <- rownames(feature_tab)[queryHits(overlap_gr)]
	segcn <- segdat_i[subjectHits(overlap_gr),"Seg_CN"]
	names(segcn)<-genenames
	dupgenes <- genenames[which(duplicated(genenames))]
	singlegenes <- setdiff(rownames(feature_tab),dupgenes)
	if(length(dupgenes) > 0){ 
		genecn[[id]]<-c(segcn[singlegenes], sapply(dupgenes,function(i) mean(segcn[which(names(segcn) == i)]) ))[rownames(feature_tab)]
		}else{
		genecn[[id]]<-segcn[rownames(feature_tab)]
		}
	}
genecn<-do.call(cbind,genecn)
ampl <- apply(genecn,1,function(x) which(x > ampl_cut))
ampl <-sapply(names(which(unlist(lapply(ampl,length)) > 0)), function(i) names(ampl[[i]]))
ddel <- apply(genecn,1,function(x) which(x < ddel_cut))
ddel <-sapply(names(which(unlist(lapply(ddel,length)) > 0)), function(i) names(ddel[[i]]))
geneCopyNum <- list(gene_CN=genecn, gene_ampl=ampl, gene_ddel =ddel)
return(geneCopyNum)
}



readDepthBreaks <- function(segdata,cutoff=0.152,segsize=10000,lowcov=NULL){

# Identify all breakpoints given a CN delta
#segdata <- your_seg_data 	# columns: sample, chromosome, start, stop, CN
#cutoff <- 0.152		# which represents a 10% copy number change for a log2 ratio
#cutoff <- 0.2		# which represents a 10% copy number change for diploid regions

segdata <- segdata[which(segdata[,4]-segdata[,3] > segsize),]
if(length(grep("chr",segdata[1,2])) != 1) segdata[,2] <- paste("chr",segdata[,2],sep="")

cn_breaks  <- list()
for(i in unique(segdata[,1])){
	sample_id <- i
	message(i)
	segdata_i <- segdata[which(segdata[,1] == i),]
	for(chr in paste("chr",c(1:22,"X"),sep="") ){
		segdata_i_chr <- segdata_i[which(segdata_i[,2] == chr),]
		if(nrow(segdata_i_chr) > 1){
			delta_all <-  segdata_i_chr[1:nrow(segdata_i_chr)-1,6] - segdata_i_chr[2:nrow(segdata_i_chr),6]
			brkpos1 <- segdata_i_chr[which( abs(delta_all) > cutoff) , 4]
			brkpos2 <- segdata_i_chr[which( abs(delta_all) > cutoff)+1, 3]
			delta <- delta_all[which( abs(delta_all) > cutoff)]
			if(length(brkpos1) > 1) {
				cn_breaks[[paste(i,chr)]]<- cbind(rep(sample_id,length(brkpos1)),rep(chr,length(brkpos1)),brkpos1,brkpos2,delta)
			}else if(length(brkpos1) == 1){
				cn_breaks[[paste(i,chr)]]<- cbind(sample_id,chr,brkpos1,brkpos2,delta)
			}
		}
	}
}

breakpoints <- data.frame(do.call(rbind,cn_breaks))
colnames(breakpoints) <- c("sample","chr","start","stop","delta")
segdataB<-breakpoints[which(!breakpoints[,"chr"] %in% c("chrY","chrM")),]
breakpoints[,"chr"] <- as.character(breakpoints[,"chr"] )
breakpoints[,"start"] <- as.numeric(as.character(breakpoints[,"start"] ))
breakpoints[,"stop"] <- as.numeric(as.character(breakpoints[,"stop"] ))
breakpoints[,"sample"] <- as.character(breakpoints[,"sample"] )
breakpoints[,"delta"] <- as.numeric(as.character(breakpoints[,"delta"] ))

if(!is.null(lowcov)){
	low_cov_GR = with(low_cov, GRanges(chr, IRanges(start=start, end=stop)))
	breakpoints_GR = with(breakpoints, GRanges(chr, IRanges(start=start, end=stop)))
	overlapgr <- GenomicAlignments::findOverlaps(breakpoints_GR,low_cov_GR,ignore.strand=TRUE)
	breakpoints<-breakpoints[setdiff(1:nrow(breakpoints),queryHits(overlapgr)),]
}

return(breakpoints)
}


########################################
##### recurrently altered genes using breakpoint data

#feature_tab<- rbind(mir_tab,genes_tab)
#exons_tab <- rbind(mir_tab[,colnames(exons_tab)],exons_tab)
#upstr <- 100000
#dnstr <- 25000
#copynumsize <- 2000000

enrichBP <- function(breaks,feature_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100){

feature_tab[which(feature_tab$strand == "+"),"start"] <- feature_tab[which(feature_tab$strand == "+"),"start"] -promoter
feature_tab[which(feature_tab$strand == "+"),"end"] <- feature_tab[which(feature_tab$strand == "+"),"end"] + offset
feature_tab[which(feature_tab$strand == "-"),"start"] <- feature_tab[which(feature_tab$strand == "-"),"start"] -offset
feature_tab[which(feature_tab$strand == "-"),"end"] <- feature_tab[which(feature_tab$strand == "-"),"end"] + promoter

features_tab_GR = with(feature_tab[,c("seqnames","start","end","strand")], GRanges(seqnames, IRanges(start=start, end=end),Rle(strand)))
mcols(features_tab_GR)$name2 <-rownames(feature_tab)


message("# Create ranges object for upstream and dnstream gene features")
upstream <- feature_tab
startup <- feature_tab$start[which(feature_tab$strand == "+")] - upstr
startup[which(startup<0)]<- 0
upstream[which(feature_tab$strand == "+"),"start"] <- startup
upstream[which(feature_tab$strand == "+"),"end"] <- feature_tab$start[which(feature_tab$strand == "+")] 
upstream[which(feature_tab$strand == "-"),"end"] <- feature_tab$end[which(feature_tab$strand == "-")] + upstr
upstream[which(feature_tab$strand == "-"),"start"] <- feature_tab$end[which(feature_tab$strand == "-")]
dnstream <- feature_tab
dnstream[which(feature_tab$strand == "+"),"start"] <- feature_tab$start[which(feature_tab$strand == "+")] 
dnstream[which(feature_tab$strand == "+"),"end"] <- feature_tab$start[which(feature_tab$strand == "+")] + dnstr
dnstream[which(feature_tab$strand == "-"),"end"] <- feature_tab$end[which(feature_tab$strand == "-")] 
startdn<-feature_tab$end[which(feature_tab$strand == "-")] - dnstr
startdn[which(startdn<0)]<- 0
dnstream[which(feature_tab$strand == "-"),"start"] <- startdn

upstream_GR = with(rbind(upstream)[,c("seqnames","start","end","strand")], GRanges(seqnames, IRanges(start=start, end=end),Rle(strand)))
mcols(upstream_GR)$name2 <-rownames(rbind(upstream))
dnstream_GR = with(rbind(dnstream)[,c("seqnames","start","end","strand")], GRanges(seqnames, IRanges(start=start, end=end),Rle(strand)))
mcols(dnstream_GR)$name2 <-rownames(rbind(dnstream))

#segments_cgi_cn  <- segdat[intersect(which(segdat[,4]-segdat[,3] < 2000000),which(abs(segdat[,6]) > 1)),]
#segments_cgi_cn[,"Chromosome"] <- paste("chr",segments_cgi_cn[,"Chromosome"],sep="")
#segments_cgi_cn  <- segments_cgi_cn[which(segments_cgi_cn$Chromosome != "chrY"),]
#segments_cgi_cn_GR = with(segments_cgi_cn[,c("Chromosome","Start","End")], GRanges(Chromosome, IRanges(start=Start, end=End)))

breaks_GR = with(breaks[,c("chr","start","stop")], GRanges(chr, IRanges(start=start, end=stop)))

overlapGenes <- GenomicAlignments::findOverlaps(features_tab_GR,breaks_GR,ignore.strand=TRUE)

overlapUpstream <- GenomicAlignments::findOverlaps(upstream_GR,breaks_GR,ignore.strand=TRUE)
overlapDnstream <- GenomicAlignments::findOverlaps(dnstream_GR,breaks_GR,ignore.strand=TRUE)


genedat <- cbind(rownames(breaks)[subjectHits(overlapGenes)],genes_tab[queryHits(overlapGenes),"gene_id"])
breakSegments<-list()
for(i in unique(genedat[,2])){
	breakSegments[[i]] <-unique(genedat[which(genedat[,2] == i),1])
}
breakSamples <- lapply(breakSegments,function(x) unique(sapply(x, function(i) strsplit(i,":")[[1]][1])) )
#sort(unlist(lapply(breakSamples,length)))

updat <- cbind(rownames(breaks)[subjectHits(overlapUpstream)],upstream[queryHits(overlapUpstream),"gene_id"])
upstreamSegments<-list()
for(i in unique(updat[,2])){
	upstreamSegments[[i]] <-unique(updat[which(updat[,2] == i),1])
}
upstreamSamples <- lapply(upstreamSegments,function(x) unique(sapply(x, function(i) strsplit(i,":")[[1]][1])) )
#sort(unlist(lapply(upstreamSamples,length)))

dndat <- cbind(rownames(breaks)[subjectHits(overlapDnstream)],upstream[queryHits(overlapDnstream),"gene_id"])
downstreamSegments<-list()
for(i in unique(dndat[,2])){
	downstreamSegments[[i]] <-unique(dndat[which(dndat[,2] == i),1])
}
downstreamSamples <- lapply(downstreamSegments,function(x) unique(sapply(x, function(i) strsplit(i,":")[[1]][1])) )
#sort(unlist(lapply(downstreamSamples,length)))

proximalSegments <-  merge2lists(upstreamSegments,downstreamSegments)
proximalSamples <-  merge2lists(upstreamSamples,downstreamSamples)
#sort(unlist(lapply(proximalSamples,length)))

results<-list()
results[["breaks"]] <-  breaks
results[["feature_tab"]] <-  feature_tab
results[["proximalSegments"]] <-  proximalSegments
results[["proximalSamples"]] <-  proximalSamples
results[["breakSegments"]] <-  breakSegments
results[["breakSamples"]] <-  breakSamples

return(results)

}


########################################
##### recurrently altered genes


#feature_tab<- rbind(mir_tab,genes_tab)
#exons_tab <- rbind(mir_tab[,colnames(exons_tab)],exons_tab)
#upstr <- 100000
#dnstr <- 25000
#copynumsize <- 2000000


enrichSV <- function(sv_df,feature_tab,exons_tab, upstr=upstr, dnstr=dnstr, promoter=1000, offset=100,copynumsize = 2000000){

# sv_df: 	A data.frame with junctions from CGI fomr multiple samples
# genes_tab:	Gene definitions with columns: seqnames,start,end,width,strand  and rownames = gene symbol ifs
#       seqnames    start      end width strand gene_id
#A1BG       chr19 58858172 58864865  6694      -    A1BG
# upstr:	bases upstream each gene to look for proximal SV
# dnstr:	bases downstream each gene to look for proximal SV

cn_df <- sv_df[which(sv_df$Type %in% c("tandem-duplication","deletion")),]
cn_df <- cn_df[which(cn_df$RightPosition - cn_df$LeftPosition < copynumsize),]
# consider segmental those svs that are duppl and deletions and complex small withing a chromosome
segids <- c(intersect(intersect(which(sv_df$Type == "complex"),which(unlist(lapply(apply(sv_df[,c("LeftChr","RightChr")],1,unique),length)) == 1)),
			which(sv_df$RightPosition - sv_df$LeftPosition < 1000000)),
				which(sv_df$Type %in% c("tandem-duplication","deletion","inversion","probable-inversion")))
seg_df  <- sv_df[segids,]
typelist <- unique(sv_df$Type)


message("# Generating all ranges objects")
# We will create genomic ranges objects with left and right junction sides 
leftJunction <- data.frame(sv_df$LeftChr, sv_df$LeftPosition, sv_df$LeftPosition + sv_df$LeftLength,sv_df$LeftStrand)
rightJunction <- data.frame(sv_df$RightChr, sv_df$RightPosition, sv_df$RightPosition + sv_df$RightLength,sv_df$RightStrand)
copynumEvent <- data.frame(cn_df$LeftChr, cn_df$LeftPosition, cn_df$RightPosition, cn_df$LeftStrand)
segmentEvent <- data.frame(seg_df$LeftChr, seg_df$LeftPosition, seg_df$RightPosition, seg_df$LeftStrand)
colnames(leftJunction)<-colnames(rightJunction) <-colnames(copynumEvent) <-colnames(segmentEvent) <- c("chr","start","end","strand")



leftJunctionGR = with(leftJunction[,c("chr","start","end","strand")], GRanges(chr, IRanges(start=start, end=end),Rle(strand)))
strand(leftJunctionGR) <- leftJunction$strand
mcols(leftJunctionGR)$name2 <-rownames(sv_df)
mcols(leftJunctionGR)$type <- sv_df$Type

rightJunctionGR = with(rightJunction[,c("chr","start","end","strand")], GRanges(chr, IRanges(start=start, end=end),Rle(strand)))
strand(rightJunctionGR) <- leftJunction$strand
mcols(rightJunctionGR)$name2 <-rownames(sv_df)
mcols(rightJunctionGR)$type <- sv_df$Type

segmentEventGR = with(segmentEvent[,c("chr","start","end","strand")], GRanges(chr, IRanges(start=start, end=end),Rle(strand)))
strand(segmentEventGR) <- segmentEvent$strand
mcols(segmentEventGR)$name2 <-rownames(seg_df)
mcols(segmentEventGR)$type <- seg_df$Type

copynumEventGR = with(copynumEvent[,c("chr","start","end","strand")], GRanges(chr, IRanges(start=start, end=end),Rle(strand)))
strand(copynumEventGR) <- copynumEvent$strand
mcols(copynumEventGR)$name2 <-rownames(cn_df)
mcols(copynumEventGR)$type <- cn_df$Type

# create genomic ranges with gene annotation coordinates
feature_tab[which(feature_tab$strand == "+"),"start"] <- feature_tab[which(feature_tab$strand == "+"),"start"] -promoter
feature_tab[which(feature_tab$strand == "+"),"end"] <- feature_tab[which(feature_tab$strand == "+"),"end"] + offset
feature_tab[which(feature_tab$strand == "-"),"start"] <- feature_tab[which(feature_tab$strand == "-"),"start"] -offset
feature_tab[which(feature_tab$strand == "-"),"end"] <- feature_tab[which(feature_tab$strand == "-"),"end"] + promoter

features_tab_GR = with(feature_tab[,c("seqnames","start","end","strand")], GRanges(seqnames, IRanges(start=start, end=end),Rle(strand)))
mcols(features_tab_GR)$name2 <-rownames(feature_tab)

exons_tab_GR = with(exons_tab[,c("seqnames","start","end")], GRanges(seqnames, IRanges(start=start, end=end)))
mcols(exons_tab_GR)$name2 <- exons_tab$gene_id

message("# Create ranges object for upstream and dnstream gene features")
upstream <- feature_tab
startup <- feature_tab$start[which(feature_tab$strand == "+")] - upstr
startup[which(startup<0)]<- 0
upstream[which(feature_tab$strand == "+"),"start"] <- startup
upstream[which(feature_tab$strand == "+"),"end"] <- feature_tab$start[which(feature_tab$strand == "+")] 
upstream[which(feature_tab$strand == "-"),"end"] <- feature_tab$end[which(feature_tab$strand == "-")] + upstr
upstream[which(feature_tab$strand == "-"),"start"] <- feature_tab$end[which(feature_tab$strand == "-")]
dnstream <- feature_tab
dnstream[which(feature_tab$strand == "+"),"start"] <- feature_tab$start[which(feature_tab$strand == "+")] 
dnstream[which(feature_tab$strand == "+"),"end"] <- feature_tab$start[which(feature_tab$strand == "+")] + dnstr
dnstream[which(feature_tab$strand == "-"),"end"] <- feature_tab$end[which(feature_tab$strand == "-")] 
startdn<-feature_tab$end[which(feature_tab$strand == "-")] - dnstr
startdn[which(startdn<0)]<- 0
dnstream[which(feature_tab$strand == "-"),"start"] <- startdn

upstream_GR = with(rbind(upstream)[,c("seqnames","start","end","strand")], GRanges(seqnames, IRanges(start=start, end=end),Rle(strand)))
mcols(upstream_GR)$name2 <-rownames(rbind(upstream))
dnstream_GR = with(rbind(dnstream)[,c("seqnames","start","end","strand")], GRanges(seqnames, IRanges(start=start, end=end),Rle(strand)))
mcols(dnstream_GR)$name2 <-rownames(rbind(dnstream))

# identify segments that don't overlap with any exon 
segmentOverlaps = GenomicAlignments::findOverlaps(exons_tab_GR,segmentEventGR,ignore.strand=TRUE,type="any",minoverlap=1)
segmentIntron <- setdiff(1:length(segmentEventGR@seqnames), unique(as.data.frame(segmentOverlaps)[,2]))
segmentExon <- intersect(1:length(segmentEventGR@seqnames), unique(as.data.frame(segmentOverlaps)[,2]))
segmentIntronIds <- rownames(seg_df[segmentIntron,])
segmentExonIds <- rownames(seg_df[segmentExon,])

###  IDENTIFICATION OF OVERLAPS BETWEEN FEATURES AND SVs

message("# find disrupting overlaps with left and right junction")
leftOverlaps = GenomicAlignments::findOverlaps(features_tab_GR,leftJunctionGR,ignore.strand=TRUE,type="any",minoverlap=1)
leftResults <-as.data.frame(leftOverlaps)
leftHits <- sort(table(queryHits(leftOverlaps)))[which( sort(table(queryHits(leftOverlaps))) > 1)]
leftSamples <- sapply(names(leftHits),function(i) substr(rownames(sv_df)[leftResults$subjectHits[which(leftResults$queryHits == i)]],0,24),simplify=FALSE)
leftJunctions <- sapply(names(leftHits),function(i) rownames(sv_df)[leftResults$subjectHits[which(leftResults$queryHits == i)]],simplify=FALSE)
leftTypes <- sapply(names(leftHits),function(i) sv_df[leftResults$subjectHits[which(leftResults$queryHits == i)],"Type"],simplify=FALSE)
names(leftJunctions) <- names(leftTypes) <- rownames(feature_tab)[as.numeric(names(leftJunctions))]
for(i in names(leftTypes)) names(leftTypes[[i]]) <- leftJunctions[[i]]
leftTypesExon <- sapply(leftTypes,function(i) i[setdiff(names(i),segmentIntronIds)],simplify=FALSE)
leftTypesIntron <- sapply(leftTypes,function(i) i[intersect(names(i),segmentIntronIds)],simplify=FALSE)

#find disrupting overlaps with right junction
rightOverlaps = GenomicAlignments::findOverlaps(features_tab_GR,rightJunctionGR,ignore.strand=TRUE,type="any",minoverlap=1)
rightResults <-as.data.frame(rightOverlaps)
rightHits <- sort(table(queryHits(rightOverlaps)))[which( sort(table(queryHits(rightOverlaps))) > 1)]
rightSamples <- sapply(names(rightHits),function(i) substr(rownames(sv_df)[rightResults$subjectHits[which(rightResults$queryHits == i)]],0,24),simplify=FALSE)
rightJunctions <- sapply(names(rightHits),function(i) rownames(sv_df)[rightResults$subjectHits[which(rightResults$queryHits == i)]],simplify=FALSE)
rightTypes <- sapply(names(rightHits),function(i) sv_df[rightResults$subjectHits[which(rightResults$queryHits == i)],"Type"],simplify=FALSE )
names(rightJunctions) <-names(rightTypes) <- rownames(feature_tab)[as.numeric(names(rightJunctions))]
for(i in names(rightTypes)) names(rightTypes[[i]]) <- rightJunctions[[i]]
# separate intronic SVs from the rest
rightTypesExon <- sapply(rightTypes,function(i) i[setdiff(names(i),segmentIntronIds)],simplify=FALSE)
rightTypesIntron <- sapply(rightTypes,function(i) i[intersect(names(i),segmentIntronIds)],simplify=FALSE)

disruptTypesExon <- list()
for(i in unique(c(names(leftTypesExon),names(rightTypesExon)))){
	disruptTypesExon[[i]] <- c(leftTypesExon[[i]],rightTypesExon[[i]])[which(!duplicated(c(names(leftTypesExon[[i]]),names(rightTypesExon[[i]]))))]
	}
disruptJunctionsExon <- merge2lists(sapply(leftTypesExon,names,simplify=FALSE),sapply(rightTypesExon,names,simplify=FALSE))
disruptSamplesExon <- lapply(disruptJunctionsExon, function(i) unique(substr(i,0,24)))

disruptTypesIntron <- list()
for(i in unique(c(names(leftTypesIntron),names(rightTypesIntron)))){
	disruptTypesIntron[[i]] <- c(leftTypesIntron[[i]],rightTypesIntron[[i]])[which(!duplicated(c(names(leftTypesIntron[[i]]),names(rightTypesIntron[[i]]))))]
	}
disruptJunctionsIntron <- merge2lists(sapply(leftTypesIntron,names,simplify=FALSE),sapply(rightTypesIntron,names,simplify=FALSE))
disruptSamplesIntron <- lapply(disruptJunctionsIntron, function(i) unique(substr(i,0,24)))



#sort(unlist(lapply(disruptSamplesExon,length)))
#sort(unlist(lapply(disruptSamplesIntron,length)))

message("# find focal copy number containing whole genes")
copynumOverlaps = GenomicAlignments::findOverlaps(features_tab_GR,copynumEventGR,ignore.strand=TRUE,type="within",minoverlap=1)
copynumResults <- as.data.frame(copynumOverlaps)
copynumHits <- sort(table(queryHits(copynumOverlaps)))[which( sort(table(queryHits(copynumOverlaps))) > 1)]
copynumSamples <- sapply(names(copynumHits),function(i) substr(rownames(cn_df)[copynumResults$subjectHits[which(copynumResults$queryHits == i)]],0,24),simplify=FALSE)
copynumJunctions <- sapply(names(copynumHits),function(i) rownames(cn_df)[copynumResults$subjectHits[which(copynumResults$queryHits == i)]],simplify=FALSE)
copynumTypes <- sapply(names(copynumHits),function(i) cn_df[copynumResults$subjectHits[which(copynumResults$queryHits == i)],"Type"] ,simplify=FALSE)
names(copynumJunctions) <- names(copynumTypes)  <- names(copynumSamples) <- rownames(feature_tab)[as.numeric(names(copynumJunctions))]
for(i in names(copynumTypes)) names(copynumTypes[[i]]) <- copynumJunctions[[i]]
copynumSamples<-lapply(copynumSamples,unique,simplify=FALSE)

message("#find upstream overlaps with left junction")
leftUpstreamOverlaps = GenomicAlignments::findOverlaps(upstream_GR,leftJunctionGR,ignore.strand=TRUE)
leftResults <-as.data.frame(leftUpstreamOverlaps)
leftHits <- sort(table(queryHits(leftUpstreamOverlaps)))[which( sort(table(queryHits(leftUpstreamOverlaps))) > 1)]
leftUpstreamSamples <- sapply(names(leftHits),function(i) unique(substr(rownames(sv_df)[leftResults$subjectHits[which(leftResults$queryHits == i)]],0,24)),simplify=FALSE)
leftUpstreamTypes <- sapply(names(leftHits),function(i) sv_df[leftResults$subjectHits[which(leftResults$queryHits == i)],"Type"],simplify=FALSE )
leftUpstreamJunctions <- sapply(names(leftHits),function(i) rownames(sv_df)[leftResults$subjectHits[which(leftResults$queryHits == i)]],simplify=FALSE)
names(leftUpstreamJunctions) <- names(leftUpstreamTypes)<- names(leftUpstreamSamples) <- rownames(feature_tab)[as.numeric(names(leftUpstreamJunctions))]
for(i in names(leftUpstreamTypes)) names(leftUpstreamTypes[[i]]) <- leftUpstreamJunctions[[i]]

message("#find dnstream overlaps with right junction")
rightUpstreamOverlaps = GenomicAlignments::findOverlaps(upstream_GR,rightJunctionGR,ignore.strand=TRUE)
rightResults <-as.data.frame(rightUpstreamOverlaps)
rightHits <- sort(table(queryHits(rightUpstreamOverlaps)))[which( sort(table(queryHits(rightUpstreamOverlaps))) > 2)]
rightUpstreamSamples <- sapply(names(rightHits),function(i) unique(substr(rownames(sv_df)[rightResults$subjectHits[which(rightResults$queryHits == i)]],0,24)),simplify=FALSE)
rightUpstreamTypes <- sapply(names(rightHits),function(i) sv_df[rightResults$subjectHits[which(rightResults$queryHits == i)],"Type"],simplify=FALSE )
rightUpstreamJunctions <- sapply(names(rightHits),function(i) rownames(sv_df)[rightResults$subjectHits[which(rightResults$queryHits == i)]],simplify=FALSE)
names(rightUpstreamJunctions) <- names(rightUpstreamTypes)<- names(rightUpstreamSamples) <- rownames(feature_tab)[as.numeric(names(rightUpstreamJunctions))]
for(i in names(rightUpstreamTypes)) names(rightUpstreamTypes[[i]]) <- rightUpstreamJunctions[[i]]

message("#find Dnstream overlaps with left junction")
leftDnstreamOverlaps = GenomicAlignments::findOverlaps(dnstream_GR,leftJunctionGR,ignore.strand=TRUE)
leftResults <-as.data.frame(leftDnstreamOverlaps)
leftHits <- sort(table(queryHits(leftDnstreamOverlaps)))[which( sort(table(queryHits(leftDnstreamOverlaps))) > 1)]
leftDnstreamSamples <- sapply(names(leftHits),function(i) unique(substr(rownames(sv_df)[leftResults$subjectHits[which(leftResults$queryHits == i)]],0,24)),simplify=FALSE)
leftDnstreamTypes <- sapply(names(leftHits),function(i) sv_df[leftResults$subjectHits[which(leftResults$queryHits == i)],"Type"] ,simplify=FALSE)
leftDnstreamJunctions <- sapply(names(leftHits),function(i) rownames(sv_df)[leftResults$subjectHits[which(leftResults$queryHits == i)]],simplify=FALSE)
names(leftDnstreamJunctions) <- names(leftDnstreamTypes)<- names(leftDnstreamSamples) <- rownames(feature_tab)[as.numeric(names(leftDnstreamJunctions))]
for(i in names(leftDnstreamTypes)) names(leftDnstreamTypes[[i]]) <- leftDnstreamJunctions[[i]]

message("#find dnstream overlaps with right junction")
rightDnstreamOverlaps = GenomicAlignments::findOverlaps(dnstream_GR,rightJunctionGR,ignore.strand=TRUE)
rightResults <- as.data.frame(rightDnstreamOverlaps)
rightHits <- sort(table(queryHits(rightDnstreamOverlaps)))[which( sort(table(queryHits(rightDnstreamOverlaps))) > 2)]
rightDnstreamSamples <- sapply(names(rightHits),function(i) unique(substr(rownames(sv_df)[rightResults$subjectHits[which(rightResults$queryHits == i)]],0,24)),simplify=FALSE)
rightDnstreamTypes <- sapply(names(rightHits),function(i) sv_df[rightResults$subjectHits[which(rightResults$queryHits == i)],"Type"] ,simplify=FALSE)
rightDnstreamJunctions <- sapply(names(rightHits),function(i) rownames(sv_df)[rightResults$subjectHits[which(rightResults$queryHits == i)]],simplify=FALSE)
names(rightDnstreamJunctions) <- names(rightDnstreamTypes) <- names(rightDnstreamSamples) <- rownames(feature_tab)[as.numeric(names(rightDnstreamJunctions))]
for(i in names(rightDnstreamTypes)) names(rightDnstreamTypes[[i]]) <- rightDnstreamJunctions[[i]]

message("# join left and right proximal junction overlaps")
bothProximalJunctions <- merge2lists(merge2lists(leftUpstreamJunctions,rightUpstreamJunctions),merge2lists(leftDnstreamJunctions,rightDnstreamJunctions))
bothProximalJunctions <- sapply(intersect(names(bothProximalJunctions),c(names(disruptJunctionsExon),names(copynumJunctions) )),function(i) 
	setdiff(bothProximalJunctions[[i]],c(disruptJunctionsExon[[i]],copynumJunctions[[i]]) ),simplify=FALSE)
bothProximalJunctions <- bothProximalJunctions[names(which(lapply(bothProximalJunctions,length) > 0))]
bothProximalJunctions <- lapply(bothProximalJunctions,unique)
bothProximalSamples <- sapply(names(bothProximalJunctions),function(i) unique(substr(bothProximalJunctions[[i]],0,24)) ,simplify=FALSE)
bothProximalTypes<-sapply(names(bothProximalJunctions),function(i) table(sv_df[bothProximalJunctions[[i]],"Type"]),simplify=FALSE)
message("# For every SV group create a matrix genes vs. SV type")


TypesMatExon <- matrix(ncol=length(typelist),nrow=length(disruptTypesExon))
rownames(TypesMatExon) <- names(disruptTypesExon)
colnames(TypesMatExon) <- typelist
TypesMatExon[]<-0
for(i in names(disruptTypesExon)) TypesMatExon[i,names(table(disruptTypesExon[[i]]))] <- table(disruptTypesExon[[i]])

# segmental SVs within introns
TypesMatIntron <- matrix(ncol=length(unique(unlist(disruptTypesIntron))),nrow=length(disruptTypesIntron))
rownames(TypesMatIntron) <- names(disruptTypesIntron)
colnames(TypesMatIntron) <- unique(unlist(disruptTypesIntron))
TypesMatIntron[]<-0
for(i in names(disruptTypesIntron)) TypesMatIntron[i,names(table(disruptTypesIntron[[i]]))] <- table(disruptTypesIntron[[i]])

# focal segmental SVs comprising whole genes
TypesMatCopynum <- matrix(ncol=length(unique(unlist(copynumTypes))),nrow=length(copynumTypes))
rownames(TypesMatCopynum) <- names(copynumTypes)
colnames(TypesMatCopynum) <- unique(unlist(copynumTypes))
TypesMatCopynum[]<-0
for(i in names(copynumTypes)) TypesMatCopynum[i,names(table(copynumTypes[[i]]))] <- table(copynumTypes[[i]])

# proximal
TypesMatProximal <- matrix(ncol=length(unique(unlist(bothProximalTypes))),nrow=length(bothProximalTypes))
rownames(TypesMatProximal) <- names(bothProximalTypes)
colnames(TypesMatProximal) <- unique(unlist(bothProximalTypes))
TypesMatProximal[]<-0
for(i in names(bothProximalTypes)) TypesMatProximal[i,names(table(bothProximalTypes[[i]]))] <- table(bothProximalTypes[[i]])

message("# Returning 'results'")

results<-list()
results[["sv_df"]] <- sv_df
results[["features"]] <- feature_tab
results[["exons"]] <- exons_tab

results[["bothProximalTypesMat"]] <- TypesMatProximal
results[["bothProximalSamples"]] <- bothProximalSamples
results[["bothProximalJunctions"]] <- bothProximalJunctions

results[["copynumTypesMat"]] <- TypesMatCopynum
results[["copynumSamples"]] <- copynumSamples
results[["copynumJunctions"]] <- copynumJunctions

results[["disruptExonTypesMat"]] <- TypesMatExon
results[["disruptExonJunctions"]] <- disruptJunctionsExon
results[["disruptExonSamples"]] <- disruptSamplesExon

results[["disruptIntronTypesMat"]] <- TypesMatIntron
results[["disruptIntronJunctions"]] <- disruptJunctionsIntron
results[["disruptIntronSamples"]] <- disruptSamplesIntron

results[["param"]] <- list(upstr=upstr, dnstr=dnstr, promoter=promoter, offset=offset,copynumsize = copynumsize)
return(results)
}



# ontain copy number for a given segment in the genome
segment_CN <- function(segdat,chr,begin,finish){
# This function reports the average segment copy number provided the boundaries
chrseg <-  segdat[which(segdat[,2] == chr),]
#chrsegment<- intersect(which(chrseg[,3] <= finish),which(chrseg[,4] >= begin))
#chrseg <- chrseg[chrsegment,]
colnames(chrseg) <- c("Sample","Chromosome","Start","End","Num_markers","Seg_CN")
chrcn<-chrseg
chrcn[which(chrcn$Start < begin),"Start"] <- begin
chrcn[which(chrcn$Start > finish),"Start"] <- finish
chrcn[which(chrcn$End < begin),"End"] <- begin
chrcn[which(chrcn$End > finish),"End"] <- finish
partCN <- list()
for(i in unique(chrseg$Sample)){
	chrcnsample <- chrcn[which(chrcn$Sample == i),]
	if(nrow(chrcnsample) == 0){
		part <- chrseg[which(chrcn$Sample == i),]
		indup <- which(abs((part$End-begin)[which((part$End-begin)<0)]) == min(abs((part$End-begin)[which((part$End-begin)<0)])))
		inddn <- which((part$Start-finish) == min(abs((part$Start-finish)[which((part$Start-finish)>0)])))
		partCN[i]<- log2((mean(2*2^part[indup,"Seg_CN"],2*2^part[inddn,"Seg_CN"]))/2)
		}else if(nrow(chrcnsample) == 1){
		partCN[i]<-chrcnsample$Seg_CN
		}else{
		seglength <- ((chrcnsample$End-chrcnsample$Start)/(max(chrcnsample$End)-min(chrcnsample$Start)))*(2*2^chrcnsample$Seg_CN)
		partCN[i]<-log2(sum(seglength)/2)
		}
	}
return(unlist(partCN))
}



breakpoint_finder2 <- function(segdata,cutoff=0.152,segsize=10000){

# Identify all breakpoints given a CN delta
#segdata <- your_seg_data 	# columns: sample, chromosome, start, stop, CN
#cutoff <- 0.152		# which represents a 10% copy number change for a log2 ratio
#cutoff <- 0.2		# which represents a 10% copy number change for diploid regions

segdata <- segdata[which(segdata[,4]-segdata[,3] > segsize),]

cn_breaks  <- list()
for(i in unique(segdata[,1])){
	sample_id <- i
	message(i)
	segdata_i <- segdata[which(segdata[,1] == i),]
	for(chr in c(1:22,"X")){
		segdata_i_chr <- segdata_i[which(segdata_i[,2] == chr),]
		if(nrow(segdata_i_chr) > 1){
			delta_all <-  segdata_i_chr[1:nrow(segdata_i_chr)-1,6] - segdata_i_chr[2:nrow(segdata_i_chr),6]
			brkpos1 <- segdata_i_chr[which( abs(delta_all) > cutoff) , 4]
			brkpos2 <- segdata_i_chr[which( abs(delta_all) > cutoff)+1, 3]
			delta <- delta_all[which( abs(delta_all) > cutoff)]
			if(length(brkpos1) > 1) {
				cn_breaks[[paste(i,chr)]]<- cbind(rep(sample_id,length(brkpos1)),paste("chr",rep(chr,length(brkpos1)),sep=""),brkpos1,brkpos2,delta)
			}else if(length(brkpos1) == 1){
				cn_breaks[[paste(i,chr)]]<- cbind(sample_id,paste("chr",chr,sep=""),brkpos1,brkpos2,delta)
			}
		}
	}
}

breakpoints <- data.frame(do.call(rbind,cn_breaks))
colnames(breakpoints) <- c("sample","chr","start","stop","delta")
segdataB<-breakpoints[which(!breakpoints[,"chr"] %in% c("chrY","chrM")),]
breakpoints[,"chr"] <- as.character(breakpoints[,"chr"] )
breakpoints[,"start"] <- as.numeric(as.character(breakpoints[,"start"] ))
breakpoints[,"stop"] <- as.numeric(as.character(breakpoints[,"stop"] ))
breakpoints[,"sample"] <- as.character(breakpoints[,"sample"] )
breakpoints[,"delta"] <- as.numeric(as.character(breakpoints[,"delta"] ))

return(breakpoints)
}



##### Functions for the analysis of somatic SVs from the Complete Genomics

################################################################################
# this function runs eQTL analysis for a list of genes provided a list of mutant samples for that gene and the expression set

sv_eQTL <- function(genelist,junctions,expmat,samplelist){

eqtl <- list()
common <- intersect(samplelist,colnames(expmat))
for(gene in genelist){
	if(gene %in% rownames(expmat)){
		mut_samples <- intersect(unique(substr(junctions[[gene]],11,16)),common)
		if(length(mut_samples) > 1){
		ll <- list(mut=expmat[gene,mut_samples], wt = expmat[gene,setdiff(common,mut_samples)])
		eqtl[[gene]][["wcox"]] <- wilcox.test(ll$mut,ll$wt)$p.value
		eqtl[[gene]][["FC"]] <- mean(ll$mut)/mean(ll$wt)
		eqtl[[gene]][["N"]] <- length(mut_samples)
		}else if(length(mut_samples) == 1){
		ll <- list(mut=expmat[gene,mut_samples], wt = expmat[gene,setdiff(common,mut_samples)])
		eqtl[[gene]][["wcox"]] <- NA
		eqtl[[gene]][["FC"]] <- ll$mut/mean(ll$wt)
		eqtl[[gene]][["N"]] <- 1
		}else{
		eqtl[[gene]][["wcox"]] <- NA
		eqtl[[gene]][["FC"]] <- NA
		eqtl[[gene]][["N"]] <- 0
		}
	}else{
	eqtl[[gene]][["wcox"]] <- NA
	eqtl[[gene]][["FC"]] <- NA
	}
}
return(eqtl)
}

################################################################################
# create a function to merge 2 lists with common names

merge2lists <- function(x,y){
	# given two lists with (or without common names) returns a combined list of unique elements
bothTypes<-list()
for(i in unique(c(names(x),names(y)))){
	if(length(y[[i]]) == 0 & length(x[[i]]) > 0){			bothTypes[[i]] <- x[[i]]
	}else if(length(y[[i]]) > 0 & length(x[[i]]) == 0){	bothTypes[[i]] <- y[[i]]
	}else if(length(y[[i]]) > 0 & length(x[[i]]) > 0){	bothTypes[[i]] <- c(x[[i]],y[[i]])
	}
}
bothTypes<-lapply(bothTypes,unique)
return(bothTypes)
}

merge2listsDup <- function(x,y){
	# given two lists with (or without common names) returns a combined list of unique elements
bothTypes<-list()
for(i in unique(c(names(x),names(y)))){
	if(length(y[[i]]) == 0 & length(x[[i]]) > 0){			bothTypes[[i]] <- x[[i]]
	}else if(length(y[[i]]) > 0 & length(x[[i]]) == 0){	bothTypes[[i]] <- y[[i]]
	}else if(length(y[[i]]) > 0 & length(x[[i]]) > 0){	bothTypes[[i]] <- c(x[[i]],y[[i]])
	}
}
return(bothTypes)
}



##########
###### FUNCTIONS for genome view of genes and structural variants
gene_model <- function(chr, start, stop, refseq,
	gene=NULL,
	interval=100000,
	addgrid=TRUE,
	addtext=TRUE,
	cex.label=1,
	cex.axis=1,
	axis.line=1.5,
	name="name",
	isotype=c("NM_","NR_")
	){

isored <- unique(unlist(sapply(isotype,function(i) grep(i,refseq$name[which(refseq$name2 == gene)],value=T))))

isotype_ids <- unique(unlist(sapply(isotype,function(i) grep(i,refseq$name))))
rows <- intersect(	intersect(which(refseq$chrom == chr),which(refseq$txStart > start)),
	intersect(which(refseq$txEnd < stop),isotype_ids))

if(length(isored) == 0 ){ gene = NULL}

dat <- refseq[rows,]
iso_dat<-list()
for(i in 1:nrow(dat)) iso_dat[[dat[i,"name"]]] <- cbind(strsplit(dat[i,"exonStarts"],",")[[1]],strsplit(dat[i,"exonEnds"],",")[[1]])

ylimit <- ypos <- -length(iso_dat)/2 -1
plot(x=NULL,y=NULL,xlim=range(c(start,stop)),ylim=range(c(ylimit,-0.5)),
	xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
if(addgrid == TRUE){ grid() }
ypos <- ypos

for(iso in names(iso_dat)){
	ypos <- ypos + 0.5
	exon_rs <- refseq[which(refseq$name == iso),]

	longexon <- cbind(as.numeric(strsplit(exon_rs$exonStarts,",")[[1]]),rep(ypos,exon_rs$exonCount),
				as.numeric(strsplit(exon_rs$exonEnds,",")[[1]]),rep(ypos,exon_rs$exonCount))
	lines(matrix(c(exon_rs$txStart,exon_rs$txEnd,ypos,ypos),2,2),lwd=2)
	if(iso %in% isored){ bordercolor <- "red"; bgcolor<-"salmon"	
	}else{bordercolor <- "black"; bgcolor<-"grey"}
	for(i in 1:nrow(longexon)){ 
		polygon(rbind(
			c(longexon[i,1],ypos+0.2),
			c(longexon[i,1],ypos-0.2),
			c(longexon[i,3],ypos-0.2),
			c(longexon[i,3],ypos+0.2)
			),lwd=1,col=bgcolor,border=bordercolor)
		}
	strandpos <- longexon[1,1]-(stop-start)/100
	if(exon_rs$strand == "-" ){ points(strandpos,ypos,pch="-",cex=1,col="red")
	}else if(exon_rs$strand == "+" ) {points(strandpos,ypos,pch="+",cex=1,col="blue")}
	if(addtext == TRUE){
		namelab<-exon_rs[,name]
		text(exon_rs$txEnd,ypos,label=namelab,pos=4,cex=cex.label)
		}
	}
xlabs <- seq(round(start/10000)*10000,stop,interval)
axis(1, at = xlabs,cex.axis=2,cex=2,labels=FALSE, lwd.ticks=1.5)
mtext(xlabs,at=xlabs,side=1,line=axis.line,cex=cex.axis)

}


plotIdeTrak <- function(chr,start,stop){
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
plotTracks(ideoTrack, from = start, to = stop)
}



sv_model <- function(tab, chr, start, stop, 
	interval=100000, 
	addgrid=TRUE,
	addtext=TRUE,
	cex.label=1,
	cex.axis=1){

ylimit <- ypos <- -nrow(tab)/3 -1
plot(x=NULL,y=NULL,xlim=range(c(start,stop)),ylim=range(c(ylimit,- 1/3)),
	xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
if(addgrid == TRUE){ grid() }
ypos <- ypos

for(i in 1:nrow(tab)){ 
	ypos <- ypos + 1/3
	if(tab[i,"Type"] == "deletion"){ bordercolor <- "blue"; bgcolor<-"lightblue"	
	}else if(tab[i,"Type"] %in% c("tandem-duplication","distal-duplication") ){bordercolor <- "red"; bgcolor<-"salmon"
	}else if(tab[i,"Type"] %in% c("inversion","probable-inversion")){bordercolor <- "black"; bgcolor<-"yellow"
	}else{bordercolor <- "black"; bgcolor<-"grey"}
	polygon(rbind(
		c(tab[i,2],ypos+0.15),
		c(tab[i,2],ypos-0.15),
		c(tab[i,3],ypos-0.15),
		c(tab[i,3],ypos+0.15)
		),lwd=1,col=bgcolor,border=bordercolor)
	if(addtext == TRUE){
		sname <- substr(rownames(tab)[i],11,16)
		textpos <-tab[i,3]
		if(tab[i,3] > stop){
			textpos <- stop - (stop-start)/8
			sname<-paste(sname,tab[i,3],sep="-")
		}else if(tab[i,2] < start){
			textpos <- tab[i,3]
			sname<-paste(sname,tab[i,2],sep="-")
			}
		if(tab[i,"Type"] %in% c("interchromosomal","distal-duplication","complex")){
			text(textpos,ypos-0.05,label=paste(sname,tab[i,4],tab[i,5],sep="-"),pos=4,cex=cex.label)
		}else{
			text(textpos,ypos-0.05,label=sname,pos=4,cex=cex.label)
			}
		}
	}
xlabs <- seq(round(start/10000)*10000,stop,interval)
axis(1, at = xlabs,cex.axis=3,cex=3,labels=FALSE, lwd.ticks=1.5)
mtext(xlabs,at=xlabs,side=1,line=1.5,cex=cex.axis)
}


################################################################################
# this function runs eQTL analysis for a list of genes provided a list of mutant samples for that gene and the expression set

sv_eQTL <- function(genelist,junctions,expmat,samplelist){

eqtl <- list()
common <- intersect(samplelist,colnames(expmat))
for(gene in genelist){
	if(gene %in% rownames(expmat)){
		mut_samples <- intersect(unique(substr(junctions[[gene]],0,16)),common)
		if(length(mut_samples) > 1){
		ll <- list(mut=expmat[gene,mut_samples], wt = expmat[gene,setdiff(common,mut_samples)])
		eqtl[[gene]][["wcox"]] <- sprintf("%.2e",wilcox.test(ll$mut,ll$wt)$p.value)
		eqtl[[gene]][["FC"]] <- round(mean(ll$mut)/mean(ll$wt),3)
		eqtl[[gene]][["N"]] <- length(mut_samples)
		}else if(length(mut_samples) == 1){
		ll <- list(mut=expmat[gene,mut_samples], wt = expmat[gene,setdiff(common,mut_samples)])
		eqtl[[gene]][["wcox"]] <- NA
		eqtl[[gene]][["FC"]] <- round(ll$mut/mean(ll$wt),3)
		eqtl[[gene]][["N"]] <- 1
		}else{
		eqtl[[gene]][["wcox"]] <- NA
		eqtl[[gene]][["FC"]] <- NA
		eqtl[[gene]][["N"]] <- 0
		}
	}else{
	eqtl[[gene]][["wcox"]] <- NA
	eqtl[[gene]][["FC"]] <- NA
	}
}
return(eqtl)
}

