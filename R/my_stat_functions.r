
rowVars <- function (x)rowSums((x - rowMeans(x, na.rm=T))^2, na.rm=T)/ncol(x)

myttest <-function (x, y, mu = 0, alternative = "two.sided", welch=T)
{
    lx <- ncol(x)
    ly <- ncol(y)
    x.var <- rowVars(x)
    y.var <- rowVars(y)
    if(welch) {
      t <- as.vector((rowMeans(x, na.rm=T) - rowMeans(y, na.rm=T))/sqrt((x.var/lx+y.var/ly)))
      df <- as.vector((x.var/lx+y.var/ly)^2/( (x.var/lx)^2/(lx-1) + (y.var/ly)^2/(ly-1)))
      df[which(df >  (lx+ly-2))] <- lx+ly-2
      df[which(df <  min(lx,ly))] <- min(lx, ly)
    }
    else{
      t <- as.vector((rowMeans(x, na.rm=T) - rowMeans(y, na.rm=T))/
        (sqrt(((lx - 1) * x.var + (ly - 1) * y.var)/(lx + ly - 2))*sqrt(1/lx + 1/ly)))
      df <- lx + ly -2
    }
    names(t) <- rownames(x)
    p <- as.vector(switch(pmatch(alternative,
           c("two.sided", "greater", "less")),
           pt(abs(t), df, lower.tail = F) * 2,
           pt(t, df, lower.tail = F),
           pt(t, df, lower.tail = T)))
    names(p) <- rownames(x)
    list(statistic = t,
         p.value = p)
}

# parametric test for the difference of two proportions
z.prop <- function(x){
  x1<-x[1]
  x2<-x[2]
  n1<-x[3]
  n2<-x[4]
  numerator = (x1/n1) - (x2/n2)
  p.common = (x1+x2) / (n1+n2)
  denominator = sqrt(p.common * (1-p.common) * (1/n1 + 1/n2))
  z.prop.ris = numerator / denominator
  return(z.prop.ris)
}


resimat <- function(x,y){
# residual calculation for two matrices
# from Alex Lachmann code
rowc <- intersect(rownames(x),rownames(y))
colc <- intersect(colnames(x),colnames(y))
message(paste("residual calculation on",length(rowc),"columns and",length(rows), "rows",sep=" "))
lres <- sapply(rowc, function(i) lm(x[i,colc]~y[i,colc])$res)
return(lres)
}

# requires (limma)
limmafit <- function(expmat,classes){
	expmat<-expmat[,c(classes$amp,classes$noamp)]
	groups <- c("amp","noamp")
	samples <- as.factor(c(rep(1,length(classes$amp)),rep(2,length(classes$noamp))))
	design <- model.matrix(~ -1+samples)
	colnames(design) <- groups
	contrast.matrix <- makeContrasts(amp-noamp,levels=design)
	fit <- lmFit(expmat, design)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	return(fit2)
	}

stouff<-function(n){1-pnorm(abs(sum(n))/sqrt(length(n)))}
stouff.nes <-function(n){sum(n)/sqrt(length(n))}
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
## apply fisher transformation to a vector of correlation coefficients
ftrans <- function(rho) log((1+rho)/(1-rho))/2
fscore <- function(zv,mrv) 2*zv*mrv/(abs(zv)+abs(mrv))

p2correct <- function(x,F=0.01,method="fdr"){
# x: a vector of p.values	
# method: any method accepted by p.ajust
return( x[names(sort(abs(p.adjust(x,method=method)-F))[1])])
}

# Brown's method to correct correlation p-values??
browns_method<-function (pval,rhos){
# pval = vector of p-values to correct
# rhos = vector of correlation coefficients
int.p.chisq <- -2 * sum(log(pval)) 
mean.chisq <- 2 * length(pval)
# measure covariance
pos.cor <-  rhos[which(rhos > 0)] 
neg.cor <- rhos[which(rhos < 0)] 
cov.chisq <- c(pos.cor*(3.25+0.75*pos.cor), neg.cor*(3.27+0.71*neg.cor)) 
#cov.chisq <- comb.nbs.cor*(3.25+0.75*comb.nbs.cor)
var.chisq <- 4*length(pval) + 2*sum(cov.chisq)
# get f & c
f = (2*(mean.chisq)^2)/var.chisq
c = var.chisq/(2*(mean.chisq))
# get p-value after correction
corr.p.chisq <- int.p.chisq/c
final.pvalue <- pchisq(corr.p.chisq ,df=f,lower=FALSE) # corrected p-value
return(final.pvalue)
}

## fisher's exact test for two vector of character elements

fet<-function(a,b,tot=NULL,alternative='t'){
## a =vector a
## b = vector b
## tot = background list
## alternative 't' = two-sided; 'g'=greater; 'l'=less

if(is.null(tot)) {
	tot <- unique(c(a,b))
	warning("No background list provided; a+b-(ab) used as background")
	}
common <- length(intersect(intersect(a,b),tot))
na <- length(intersect(setdiff(a,intersect(a,b)),tot))
nb <- length(intersect(setdiff(b,intersect(a,b)),tot))
rest<- length(tot) - na - nb + common
res <- fisher.test(matrix(c(common,na,nb,rest),2,2),alternative=alternative)
return(list(p.value=res$p.value,oddsr=res$estimate,call=matrix(c(common,na,nb,rest),2,2)))
}

multi_fet <- function(set, collection, tot=NULL, alternative='t'){
if(is.null(tot)) tot <- unique(unlist(collection))
res<-list()
for(i in names(collection)){
	f <- fet(set,collection[[i]],tot,alternative=alternative)
	res[[i]] <- c(f$p.value,f$oddsr, as.numeric(f$call))
	}
ret <- data.frame(do.call(rbind,res))
fdr <- p.adjust(ret[,1],method="fdr")
ret <- cbind(fdr,ret)
colnames(ret) <- c("FDR","fet","OddsR","AB","BB","AA","tot")
return(ret[order(ret[,2]),])
}


hyp <-function(a,b,tot=NULL,alternative='t'){
## a =vector a
## b = vector b
## tot = background list
## alternative 't' = two-sided; 'g'=greater; 'l'=less

if(is.null(tot)) {
  tot <- unique(c(a,b))
  warning("No background list provided; a+b-(ab) used as background")
  }
q <- length(intersect(intersect(a,b),tot))
k <- length(intersect(a,tot))
m <- length(intersect(b,tot))
n<- length(setdiff(tot,a))
res <- phyper(q,m,n,k,lower.tail=FALSE)
return(res)
}


multi_hyp <- function(set, collection, tot=NULL,alternative='t'){
if(is.null(tot)) tot <- unique(unlist(collection))
res<-list()
for(i in names(collection)){
  h <- hyp(set,collection[[i]],tot,alternative=alternative)
  res[[i]] <- h
  }
ret <- unlist(res)
hypfdr <- p.adjust(ret,method="fdr")
ret <- cbind(hypfdr,ret)
colnames(ret) <- c("hypFDR","hyp")
return(ret[order(ret[,2]),])
}


val2col<-function(z,nbreaks=256,col.range=NULL,color=NULL){
# z = a numeric vector or matrix 
# nbreaks= number of colors in range
# col.range, a vector with three values to set the color scale: c(minimum, median, maximum)
	#	an option is: col.range <- c(quantile(z,0.01),0,quantile(z,0.99))
	#
# color, a vector with three values to set the extreme colors, default: c("darkblue","white","salmon") 
#
if(is.null(col.range)){
	maximum= max(z)+ max(z)/257
	minimum= min(z)- min(z)/257
	#maximum= max(z)
	#minimum= min(z)
	med <- median(z)
}else{
	z[which(z < col.range[1])] <- col.range[1]
	minimum <- col.range[1]-col.range[1]/100
	med<-col.range[2]
	z[which(z > col.range[3])] <- col.range[3]
	maximum <- col.range[3]+col.range[3]/100
}
if(is.null(color)) color<-c("darkblue","white","salmon")

dn_breaks <- seq(minimum, med, length = nbreaks/2)
up_brakes <- seq(med, maximum, length = nbreaks/2)
breaks <- c(dn_breaks,up_brakes[2:length(up_brakes)])
ncol <- length(breaks)
col <- colorpanel(ncol,color[1],color[2],color[3])

if(is.null(dim(z))){
	CUT <- cut(z, breaks=breaks)
	colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
	names(colorlevels)<-names(z)
	}else if(length(dim(z)) == 2){
	CUT <- cut(z, breaks=breaks)
	colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
	colorlevels[which(is.na(colorlevels))] <- color[2]
	colorlevels<-matrix(colorlevels,nrow(z),ncol(z))
	colnames(colorlevels) <-colnames(z)
	rownames() <-rownames(z)
	}
return(colorlevels)
}


## calculate the paired correlations of two matrices with same dimensions
colCors = function(x, y , d="row") {
if(d == "row"){
	x<-t(x)
	y<-t(y)
	}
sqr = function(x) x*x
if(!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y)))
stop("Please supply two matrices of equal size.")
x = sweep(x, 2, colMeans(x))
y = sweep(y, 2, colMeans(y))
cor = colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
return(cor)
}

## filter matrrix by CV
filterCV <- function (expset) 
{
    repet <- tapply(rownames(expset), factor(rownames(expset)), 
        length)
    if (max(repet) > 1) {
        d <- expset[rownames(expset) %in% (names(repet)[repet == 1]), ]
        d1 <- expset[rownames(expset) %in% (names(repet)[repet > 1]), ]
        d2 <- tapply(1:nrow(d1), factor(rownames(d1)), function(pos, d1, cv) {
            list(name = rownames(d1)[pos][which.max(cv[pos])], value = d1[pos, ][which.max(cv[pos]), ])
        }, d1 = d1, cv = frcv(d1))
        d3 <- t(sapply(d2, function(x) x$value))
        if (nrow(d3) != length(d2)) d3 <- t(d3)
        rownames(d3) <- sapply(d2, function(x) x$name)
        expset <- rbind(d, d3)
    }
    expset
}

frcv <-function (x) sqrt(frvarna(x))/frmeanna(x)

frvarna <- function (x) {
    ave <- as.vector(frmeanna(x))
    pos <- which(is.na(x))
    largo <- frlengthna(x)
    x[pos] <- rep(ave, ncol(x))[pos]
    (x - ave)^2 %*% rep(1, ncol(x))/(largo - 1)
}

frmeanna <- function (x) {
    largo <- frlengthna(x)
    x[is.na(x)] <- 0
    res <- x %*% rep(1, ncol(x))/largo
    names(res) <- rownames(x)
    res
}

frlengthna <- function (x) 
{
    r <- x/x
    r[x == 0] <- 1
    r[!is.finite(r)] <- 0
    r %*% rep(1, ncol(r))
}

## apply double transformation to expression matrix

zrank<-function(mexp){
  # transform matrix into ranked by columns and 
  # transform matrix into Z-score by columns
rankexp<-t(apply(mexp,1,rank))
zrankexp<-scale(rankexp)
return(zrankexp)
}

# transformation of a rows/cols of a matrix into a normal distrib
ntrans <- function(mexp,side=1){
if(side == 1){
	mexp <- t(apply(mexp,1,rank,ties.method="random"))
	mexp.n<-mexp
	normal <- sort(rnorm(ncol(mexp)))
	for(i in rownames(mexp)) mexp.n[i,] <- normal[mexp[i,]]
}else if(side == 2){
	mexp <- apply(mexp,2,rank,ties.method="random")
	mexp.n<-mexp
	normal <- sort(rnorm(nrow(mexp)))
	for(i in colnames(mexp)) mexp.n[,i] <- normal[mexp[,i]]
	}
return(mexp.n)
}

ztrans<-function(mexp,side=1){
  # transform matrix into Z-score by rows
mexp.z<-mexp
mexp.mean<-apply(mexp,side,mean)
mexp.sd<-apply(mexp,side,sd)
if(side == 1) for(i in 1:nrow(mexp))mexp.z[i,] <- (mexp[i,]-mexp.mean[i])/mexp.sd[i]
if(side == 2) for(i in 1:ncol(mexp))mexp.z[,i] <- (mexp[,i]-mexp.mean[i])/mexp.sd[i]
	return(mexp.z)
}

zrank_lim<-function (mexp,limit=100){
  # transform matrix into limited 1:n ranked by columns
  # transform matrix into Z-score by columns
br<-c(1:limit)
for(i in 1:limit) br[i]<-i*ncol(mexp)/limit
mexp.rank<-apply(mexp,1,rank)
mexp.rank<-apply(mexp.rank,1,function(x) as.numeric(cut(x,breaks=c(0,br), labels=c(1:limit))))
mexp.mean<-apply(mexp.rank,2,mean)
mexp.sd<-apply(mexp.rank,2,sd)
for(i in 1:ncol(mexp.rank))mexp.rank[,i] <- (mexp.rank[,i]-mexp.mean[i])/mexp.sd[i]
return(mexp.rank)
}

na2median<-function(mexp,sd_noise=0.0001){
# transform character matrix into numeric
#	mexp: a numeric matrix
#	sd_noise: standard deviation of white noise added to the median (0 if no noise to be added)

b <- t(apply(a,1,as.numeric))		
# obtain which array index have NA values
na_index <- which(apply(t(b),1,is.na), arr.ind = T)
# transform NA values into the median of its tow adding white noise
	# obtain vector with median values for each row with ne
median_list <- apply(b[na_index[,1],],1, median, na.rm=TRUE)
	# substitute NA values by the medians + white noise
b[na_index] <- median_list + rnorm(length(median_list),0,sd)
return(b)
}

# return a list with for vectors containing the names of elements splited by quantiles

quantile_names <- function(x, out = "v"){
# output options: 
#	v=vector; l=list
xcat <- factor(cut(x, quantile(x), include.lowest = TRUE), 
                  labels = c(1:4))
if(out == "l"){
output<-list(fst = names(x)[which(xcat == 1)],
snd = names(x)[which(xcat == 2)],
trd = names(x)[which(xcat == 3)],
fth = names(x)[which(xcat == 4)]) 
}else if(out =="v"){
output<-xcat
names(output)<-names(x)
} else{stop (paste("Invalid option",out))}
return(output)
}


IQM <- function(x,lowQ=0.1,upQ=0.9){
# Obtains inter quantile average for a defined 'x' vector and both lower and upper quantiles
# 'x' numeric vector 
# 'lowQ' lower quantile
# 'upQ' upper quantile
  rx <- rank(x,ties.method ='random')
  qt1<-quantile(rx,lowQ)
  qt2<-quantile(rx,upQ)
  inter_quantile_mean <- mean(x[intersect(which(rx > qt1),which(rx < qt2))])
  return(inter_quantile_mean)
}


IQSD <- function(x,lowQ=0.1,upQ=0.9){
# Obtains inter quantile average for a defined 'x' vector and both lower and upper quantiles
# 'x' numeric vector 
# 'lowQ' lower quantile
# 'upQ' upper quantile
  rx <- rank(x,ties.method ='random')
  qt1<-quantile(rx,lowQ)
  qt2<-quantile(rx,upQ)
  inter_quantile_mean <- sd(x[intersect(which(rx > qt1),which(rx < qt2))])
  return(inter_quantile_mean)
}



cor.dif <- function(na=NULL,ra=NULL,nb=NULL,rb=NULL, pval=0.05, minsize=6){
# Using the Fisher r-to-z transformation, this function calculates z values that can be applied to assess
# the significance of the difference between two correlation coefficients, ra and rb, found in two independent samples.
# If ra is greater than rb, the resulting value of z will have a positive sign; if ra is smaller than rb, the sign of z will be negative.
# na = length of sample set a
# ra = vector of correlations a
# nb = length of sample set b
# rb = vector of correlations b
        #if(na=NULL || ra=NULL || nb=NULL || rb=NULL){stop ("missing data! should provide 1:4 items")}
        if(max(round(ra, 7) ) > 1 || min(round(ra, 7)) < -1){stop("ra must fall between +1.0 and -1.0, inclusive.")}
        if(max(round(rb, 7)) > 1 || min(round(rb, 7)) < -1){stop("rb must fall between +1.0 and -1.0, inclusive.")}
        if(na < minsize || nb < minsize ){stop("n must be equal to or greater than 4.")}
        if(floor(na) < na) {stop("n_a must be an integer value.")}
        if(floor(nb) < nb) {stop("n_b must be an integer value.")}

# first we transform correlation matrix into a matrix of z values
        raplus = 1*ra+1 + 1e-5
        raminus = 1-ra + 1e-5
        rbplus = 1*rb+1 + 1e-5
        rbminus = 1-rb + 1e-5

        za = (log(raplus)-log(raminus))/2
        zb = (log(rbplus)-log(rbminus))/2

        se = sqrt((1/(na-3))+(1/(nb-3)))
        z = (za-zb)/se

#  pvalue estimation R version
        p <- 1-pnorm(abs(z))

        ## pvalue estimation javascript version
#z2 = abs(z)
#p2 =(((((.000005383*z2+.0000488906)*z2+.0000380036)*z2+.0032776263)*z2+.0211410061)*z2+.049867347)*z2+1
#p2 = p2^-16
#p1 = p2/2
#p<-p1[1,]
#z<-z[1,]

        return(list(
        score=z,
        p.value=p))

}

