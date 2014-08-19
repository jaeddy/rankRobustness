
### Install required packages ###
#################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")
# install.packages("MASS")
# biocLite("hgfocus.db")
# biocLite("sva")
# biocLite("limma")
# packrat::install_github("rafalib","ririzarr")
# packrat::install_github("dagdata","genomicsclass")

### Load libraries & data ###
#############################

library(rafalib) # useful for some of the plotting
library(dagdata)
library(genefilter) # for t-tests
data(GSE5859) # expression data
library(hgfocus.db) # chromosome data
library(sva) # for ComBat and SVA
library(limma)


### Pull out information about chromosome & sex ###
###################################################

chr <- mget(featureNames(e), hgfocusCHRLOC)
chr <- sapply(chr, function(x) {
    tmp <- names(x[1])
    ifelse(is.null(tmp), NA, paste0("chr", tmp))
    })
y <- colMeans(exprs(e)[which(chr == "chrY"), ])
hist(y, nc = 30)
sex <- ifelse(y < 5, "F", "M")


### Specify batches based on month ###
######################################

batch <- format(pData(e)$date, "%y%m")
ind <- which(batch %in% c("0506", "0510"))

set.seed(1)
N <- 12; N1 <-3; 
M <- 12; M1 <- 9
ind <- c(sample(which(batch == "0506" & sex == "F"), N1),
         sample(which(batch == "0510" & sex == "F"), N - N1),
         sample(which(batch == "0506" & sex == "M"), M1),
         sample(which(batch == "0510" & sex == "M"), M - M1))
table(batch[ind], sex[ind])

### Visualize the data with confounding ###
###########################################

set.seed(1)
tt <- genefilter::rowttests(exprs(e)[, ind], factor(batch[ind]))
ind1 <- which(chr == "chrY") ##real differences
ind2 <- setdiff(c(order(tt$dm)[1:25], order(-tt$dm)[1:25]), ind1)
ind0 <- setdiff(sample(seq(along = tt$dm), 50), c(ind2, ind1))
geneindex <- c(ind2, ind0, ind1)
mat <- exprs(e)[geneindex, ind]
mat <- mat - rowMeans(mat)
icolors <- rev(brewer.pal(11, "RdYlBu"))
mypar(1,1)
image(t(mat), xaxt = "n", yaxt = "n", col = icolors)

### Set up variables ###
########################

dat <- exprs(e)[, ind]
X <- sex[ind] # the covariate
Z <- batch[ind]
tt <- genefilter::rowttests(dat, factor(X))

HLIM <- c(0, 1500)
mypar(1, 2)
hist(tt$p[!chr %in% c("chrX", "chrY")], nc = 20, xlab = "p-value", ylim = HLIM,
     main = "")
hist(tt$p[chr %in% c("chrY")], nc = 20, xlab = "p-value", ylim = c(0, 12), 
     main = "")


### Run ComBat ###
##################

mod <- model.matrix(~ X)
cleandat <- ComBat(dat, Z, mod)
tt <- genefilter::rowttests(cleandat, factor(X))
mypar(1, 2)
hist(tt$p[!chr %in% c("chrX", "chrY")], nc = 20, xlab = "p-value", ylim = HLIM,
     main = "")
hist(tt$p[chr %in% c("chrY")], nc = 20, xlab = "p-value", ylim = c(0, 12), 
     main = "")


## Run SVA ###
##############

svafit <- sva(dat, mod)
svaX <- model.matrix(~X + svafit$sv)
lmfit <- lmFit(dat, svaX)
tt <- lmfit$coef[,2] * sqrt(lmfit$df.residual) / (2 * lmfit$sigma)
mypar(1, 2)
pval <- 2 * (1 - pt(abs(tt), lmfit$df.residual[1]))
hist(pval[!chr %in% c("chrX", "chrY")], xlab = "p-values", ylim = HLIM, 
     main = "")
hist(pval[chr %in% c("chrY")], nc = 20, xlab = "p-value", ylim = c(0, 12),
     main = "")
