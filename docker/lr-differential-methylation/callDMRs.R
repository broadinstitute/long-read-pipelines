suppressPackageStartupMessages({
    require(edgeR,quietly=T)
    require(bsseq,quietly=T)
    require(optparse,quietly=T)
    require(BiocParallel,quietly=T)
})

# Define the command line options
option_list <- list(
    make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Directory path that houses the bismark-formatted input files", metavar="character"),
    make_option(c("-o", "--outputdir"), type="character", default=NULL, help="Output directory", metavar="character"),
    make_option(c("-p", "--pedfile"), type="character", default=NULL, help="PED file for dataset", metavar="character"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of threads to use", metavar="integer")
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

message("Reading PED file")
# read in the PED file and process it to get affected sample(s)
ped <- read.table(opt$pedfile, header=F, stringsAsFactors=F)
affected <- ped[ped$V6==2,2] # get individual sample names of affected samples

message("Changing directory and reading methylation input files")
# change to the input file directory and find the bismark-formatted files
setwd(opt$dir)
file_list <- list.files(pattern = ".bismark.cov")

message("file list to read:")
print(file_list)

# read in the bismark-formatted files and set their attributes
bsobj.raw <- read.bismark(files=file_list)
message("completed reading bismark files")
bsobj <- orderBSseq(bsobj.raw) # need to make sure it's ordered for the smoothing function

samples=sub("\\.bismark.cov","",file_list)
sampleNames(bsobj) <- samples
message("sample names derived from file names:")
print(samples)
bsseq::pData(bsobj) <- data.frame(status = ifelse(samples %in% affected, "affected", "unaffected"))

message("running smoothing function")
bsobj.smooth <- BSmooth(bsobj, BPPARAM= MulticoreParam(workers = opt$threads, progressbar = TRUE), verbose=T)
message("completed smoothing")

sampleNames(bsobj.smooth)<- samples
BS.cov <- getCoverage(bsobj)
dmrs_df <- data.frame()

# switch working directory to the output directory
setwd(opt$outputdir)

# save the smoothed data for followup analysis
message("saving smoothed data to file")
saveRDS(bsobj.smooth,file="smoothed_bsobj.rds")

# loop over each affected sample to find DMRs
for (affected_sample in affected[affected %in% samples]) {
    message(paste("processing affected sample",affected_sample))
    # Identify samples in the same family as the affected sample
    same_family <- ped[ped$V1 == ped[ped$V2 == affected_sample,]$V1,]$V2

    # check that there are at least 2 samples not in the same family as the affected sample
    if (length(samples[!(samples %in% same_family)]) < 2) {
        message("Skipping: not enough samples unrelated to the affected sample")
        next
    }

    # lightly filter to remove sites with very low cov or too much missing data
    bsobj.smooth.filt <- bsobj.smooth[(which(BS.cov[,samples==affected_sample]>=5 & rowSums(BS.cov[,!(samples %in% same_family)]>=5)>=2)),]

    # run t-test with group2 = samples not in same family as affected sample
    bsobj.ttest<-BSmooth.tstat(bsobj.smooth.filt, group1=affected_sample, group2=samples[!(samples %in% same_family)], estimate.var="group2")

    # t-stat cutoff is empirical, selected after a lot of spot-checking
    dmrs0<-dmrFinder(bsobj.ttest, cutoff=c(-2.5,2.5))

    # filter: min #CpGs=5 and abs(meanDiff) at least 0.25
    dmrs <- subset(dmrs0, n >= 5 & abs(meanDiff) >= 0.25)

    # only run if there are dmrs for this sample
    if (nrow(dmrs) >= 1) {
        # add sample name to dmrs table
        dmrs$sample <- affected_sample

        # add dmrs to dmrs_df
        dmrs_df <- rbind(dmrs_df, dmrs)

        # plot many and save to file
        # set affected sample to red and control group (not in same family) to blue, others (same family) to grey
        pData(bsobj.smooth.filt)$col <- ifelse(samples==affected_sample, "red", "orange")
        pData(bsobj.smooth.filt[,!(samples %in% same_family)])$col <- "blue"
        pdf(file = paste(affected_sample,"_",dmrs$chr[1],"_","dmrs.pdf",sep=""), width = 10, height = 5)
        plotManyRegions(bsobj.smooth.filt, dmrs, extend = 5000, addRegions = dmrs)
        dev.off()
    }
}

#save DMR results to file
write.table(dmrs_df,file="dmrs.tsv",sep="\t",row.names=F,quote=F)
