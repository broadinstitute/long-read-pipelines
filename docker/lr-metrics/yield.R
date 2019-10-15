suppressPackageStartupMessages(library(tidyverse))
library(tidyverse)

options(scipen=8)

# Compute NX (N50, N90, etc.).  The vector of lengths, l, must be sorted in decreasing order.
nX <- function(l, X) {
    v = l[last(which(cumsum(as.numeric(l)) <= (X/100)*sum(as.numeric(l))))]
    return(v)
}

# Simple logging function
logger <- function(text) {
    cat(timestamp(prefix = "[", suffix = "]", quiet=TRUE), text, "\n")
}

# Read command-line info
if (interactive()) {
    file = "read_names_and_lengths.txt.gz"
    prefix = gsub(".txt.gz", "", file)
} else {
    args <- commandArgs(trailingOnly=TRUE)
    file = args[1]
    prefix = args[2]
}

read_sampling_size = 2000000

logger(paste("File:", file))
logger(paste("Prefix:", prefix))

logger("Load data...")
f = read_delim(file, delim=" ", col_names=c("read_name", "length"), col_types="ci")

logger("Compute basic stats...")
# basic info
num_reads = nrow(f)
num_bases = sum(f$length, na.rm=TRUE)
num_zmws = 0
movie_name = NA
data_type = "not_autodetected"

logger("Subsample data...")
N = ifelse(read_sampling_size < nrow(f), read_sampling_size, nrow(f))
d = f[1:N,]

if (grepl("/", d$read_name[1])) {
    # PacBio

    logger("Count number of ZMWs represented...")

    m = unique(as.integer(unlist(lapply(strsplit(f$read_name, "/", fixed=TRUE), function(x) { return(x[2]); }))))
    num_zmws = length(m)

    logger("Transform and facet PacBio data...")

    if (grepl("ccs", d$read_name[1])) {
        data_type = "pacbio_ccs"
    } else {
        data_type = "pacbio_uncorrected"
    }

    d = d %>% separate(read_name, c("movie_name", "zmw", "range"), "/") %>%
              select(-range) %>%
              mutate(zmw = as.integer(zmw))

    q = d %>% group_by(zmw) %>%
              summarize(s = as.integer(round(sum(length) / 1000) * 1000)) %>%
              group_by(s) %>%
              summarize(l = sum(s))

    l = d %>% group_by(zmw) %>% summarize(l = sum(length)) %>% pull(l) %>% sort(decreasing=TRUE)

    movie_name = d$movie_name[1]

    logger("Examine pass vs. read length distributions...")
    m = d %>% group_by(zmw) %>% summarize(np = n(), l = as.integer(median(length))) %>% group_by(np)

    ml = list()
    ml.last = ifelse(m %>% pull(np) %>% max() < 15, m %>% pull(np) %>% max(), 15)

    colors = rep("#93C700", ml.last)

    for (i in 1:ml.last) {
        ml[[i]] = m %>% filter(np == i) %>% pull(l)

        mad = median(abs(ml[[i]] - median(m$l)))
        mz = 0.6745*(median(ml[[i]]) - median(m$l))/mad

        colors[i] = ifelse(mz < -0.25, "#FFC801", "#93C700")
    }

    for (dev in c("pdf", "png")) {
        logger(paste("Plot passes vs. read length (", dev, ")...", sep=""))

        if (dev == "pdf") {
            pdf(paste(prefix, ".passes_vs_readlength.pdf", sep=""));
        } else {
            png(paste(prefix, ".passes_vs_readlength.png", sep=""), res=300, width=8, height=8, units="in");
        }

        title = paste("Pass vs. read length distributions for ", movie_name, "\n", "(", data_type, "; ", formatC(N, big.mark=","), " reads sampled)", sep="")

        par(mar=c(5, 6, 4, 2))
        boxplot(ml, outline=FALSE, frame=FALSE, varwidth=TRUE, col=colors, las=1, xlab="Number of passes", ylab="", yaxt="n", main=title)
        axis(2, at=axTicks(2), labels=formatC(axTicks(2), big.mark=','), las=1)
        mtext("Read length (bp)", 2, line=4)

        legend("bottomright", c(expression("Modified Z-score">="-0.25"), expression("Modified Z-score"<"-0.25")), fill=c("#93C700", "#FFC801"), bty="n")

        dev.off()
    }
} else {
    # Nanopore or unrecognized

    logger("Transform and facet unknown data...")

    q = d %>% summarize(s = as.integer(round(sum(length) / 1000) * 1000)) %>%
              group_by(s) %>%
              summarize(l = sum(s))

    l = d %>% summarize(l = sum(length)) %>% pull(l) %>% sort(decreasing=TRUE)
}

shortest_read_bp = min(l)
longest_read_bp = max(l)
mean_read_bp = round(mean(l))
stdev_read_bp = round(sd(l))
median_read_bp = round(median(l))

logger("Compute N10-N90...")
# read N10 ... N90
nv = c()

for (X in seq(1, 99, by=1)) {
    nv = c(nv, nX(l, X))
}

l.90 = last(which(q$s <= nv[90]))
l.70 = last(which(q$s <= nv[70]))
l.50 = last(which(q$s <= nv[50]))
l.30 = last(which(q$s <= nv[30]))
l.10 = last(which(q$s <= nv[10]))
l.0 = as.integer(nrow(q))

y.50_30 = 1.02*max(na.omit(q$l[l.50:l.30]))
y.90_50 = 1.02*max(na.omit(q$l[l.90:l.50]))
y.10_0 = 1.02*max(na.omit(q$l[l.10:nrow(q)]))

title = paste("Yield for ", movie_name, "\n", "(", data_type, "; ", formatC(N, big.mark=","), " reads sampled)", sep="")

for (dev in c("pdf", "png")) {
    logger(paste("Plot yield (", dev, ")...", sep=""))

    if (dev == "pdf") {
        pdf(paste(prefix, ".yield.pdf", sep=""));
    } else {
        png(paste(prefix, ".yield.png", sep=""), res=300, width=8, height=8, units="in");
    }

    par(mar=c(5, 9, 4, 2))

    plot(q$s, q$l, type="l", lwd=0.1, bty="n", xlab="Read Length (bp)", ylab="", xlim=c(0, 1.2*max(na.omit(q$s))), ylim=c(0, 1.1*max(na.omit(q$l))), xaxt="n", yaxt="n", main=title)
    axis(1, at=axTicks(1), labels=formatC(axTicks(1), big.mark=','))
    axis(2, at=axTicks(2), labels=formatC(axTicks(2), big.mark=','), las=1)
    mtext("Density\n(Bases Read per Read Length)", 2, line=ceiling(log10(max(axTicks(2)))) - 2)

    segments(nv[50], 100, nv[50], y.50_30, lty=3, lwd=1)
    segments(nv[50], y.50_30, nv[30], y.50_30, lty=3, lwd=1)
    arrows(q$s[l.30] - 1000, y.50_30, q$s[l.30], y.50_30, lwd=1, length=0.05)
    text(nv[50], y.50_30, labels=paste("50% of bases in reads:\n>", l.50, "kbp\n"), pos=4, cex=0.9, offset=-0.10, adj=c(1, 1))

    segments(nv[90], 100, nv[90], y.90_50, lty=3, lwd=1)
    segments(nv[90], y.90_50, nv[70], y.90_50, lty=3, lwd=1)
    arrows(q$s[l.70] - 1000, y.90_50, q$s[l.70], y.90_50, lwd=1, length=0.05)
    text(nv[90], y.90_50, labels=paste("90% of bases in reads:\n>", l.90, "kbp\n"), pos=4, cex=0.9, offset=-0.10, adj=c(1, 1))

    segments(nv[10], 100, nv[10], y.10_0, lty=3, lwd=1)
    segments(nv[10], y.10_0, q$s[nrow(q)], y.10_0, lty=3, lwd=1)
    arrows(q$s[nrow(q)] - 1000, y.10_0, q$s[nrow(q)], y.10_0, lwd=1, length=0.05)
    text(nv[10], y.10_0, labels=paste("10% of bases in reads:\n>", l.10, "kbp\n"), pos=4, cex=0.9, offset=-0.10, adj=c(1, 1))

    polygon(c(q$s[1:l.0], q$s[l.0] + 1000), c(q$l[1:l.0], 0), col="#334D5C", border=NA)
    polygon(q$s[1:(l.10 + 1)], c(q$l[1:l.10], 0), col="#45B29D", border=NA)
    polygon(q$s[1:(l.50 + 1)], c(q$l[1:l.50], 0), col="#EFC94C", border=NA)
    polygon(q$s[1:(l.90 + 1)], c(q$l[1:l.90], 0), col="#DF5A49", border=NA)

    points(q$s, q$l, type="l", lwd=1.5)

    dev.off();
}

# Write tables
logger("Write stats table...")
yield_stats = tribble(
    ~file, ~prefix, ~data_type, ~movie_name, ~num_reads, ~num_bases, ~num_zmws, ~read_sampling_size, ~shortest_read_bp, ~longest_read_bp, ~mean_read_bp, ~stdev_read_bp, ~median_read_bp, ~N10,    ~N50,   ~N90,
     file,  prefix,  data_type,  movie_name,  num_reads,  num_bases,  num_zmws,  read_sampling_size,  shortest_read_bp,  longest_read_bp,  mean_read_bp,  stdev_read_bp,  median_read_bp,  nv[10], nv[50],  nv[90]
)
write_tsv(yield_stats, paste(prefix, ".yield.stats.txt", sep=""))

logger("Write N(X) table...")
as_tibble(cbind(nx = 1:99, nv)) %>% write_tsv(paste(prefix, ".read_lengths.nX.txt", sep=""))

logger("Done.")
