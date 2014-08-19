
# open libraries
library(affy)
library(GEOquery)
library(stringr)

# create master list of GSE accession IDs
fileList <- read.table("data/file_list.txt", stringsAsFactors = FALSE)
labelFiles <- grepl("label", fileList$V1)
gseList <- fileList[!labelFiles, 1]
gseList <- gsub(".txt", "", gseList)

# create folder to store processed data sets
if (!file.exists("data/processed")) {
    dir.create("data/processed")
}

for (gse_i in gseList) {
    # download RAW files from GEO ftp
    gseID <- gseList[gse_i]
    
    # create a temporary directory to store .CEL files
    if (!file.exists("data/tmp")) {
        dir.create("data/tmp")
    }
    
    # will save GSE IDs when no RAW files exist; clear this file if it exists
    # from previous run
    if (!file.exists("data/no_raw.txt")) {
        file.remove("data/no_raw.txt")
    }
    
    # build the address string for the current GSE
    baseAddress <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/"
    
    createSubDirStr <- function(idStr) {
        if (str_length(idStr) < 6) {
            str_sub(idStr, -2, -1) <- "nnn/"
        } else {
            str_sub(idStr, -3, -1) <- "nnn/"
        }
        idStr
    }
    subDir <- createSubDirStr(gseID)
    rawTag <- paste0(gseID, "/suppl/", gseID, "_RAW.tar")
    rawAddress <- paste0(baseAddress, subDir, rawTag)
    
    destAddress <- paste0("data/tmp/", gseID, "_RAW.tar")
    download.file(rawAddress, 
                  destfile = destAddress,
                  method = "curl")
    
    if (!file.exists(destAddress)) {
        write(gseID, "data/no_raw.txt", append = TRUE)
        warning(paste("No RAW files found for", gseID))
    }
    
    # untar the downloaded file
    celAddress <- paste0("data/tmp/", gseID)
    untar(destAddress, exdir = celAddress)
    
    # read CEL files into AffyBatch object
    affyData <- ReadAffy(celfile.path = celAddress)
    samples <- gsub(".CEL.gz", "", sampleNames(affyData))
    
    # normalize expression data using RMA
    normData <- rma(affydata)
    normX <- exprs(normData)
    
    # convert expression data to data.frame object
    gseDat <- as.data.frame(normX)
    colnames(gseDat) <- samples
    
    # extract feature data from GEO soft file
    geoq <- getGEO(gseList[10], destdir = "data/")
    e <- geoq[[1]]
    
    fData <- featureData(e)
    probe <- featureNames(e)
    symbol <- fData$"Gene Symbol"
    entrez <- fData$"ENTREZ_GENE_ID"
    
    # combine feature information into a second data.frame
    featDat <- data.frame(probe = probe, symbol = symbol, entrez = entrez)
    
    # merge data frames
    mergeDat <- merge(featDat, gseDat, by.x = "probe", by.y = "row.names")
    
    # save the merged data frame to a new text file
    write.table(mergeDat, paste0("data/processed/", gseID, ".txt"))
    
    # delete the temporary folder
    unlink("data/tmp", recursive = "TRUE")
}

