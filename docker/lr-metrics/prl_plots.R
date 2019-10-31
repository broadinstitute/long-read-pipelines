source("utils.R")

default_args = c(
    "--sample_name", "HG02982",
    "--flowcell_name", "m64020_190307_140324",
    "--flowcell_name", "m64020_190816_194042",
    "--flowcell_name", "m64020_190818_014714",
    "testing/per_flowcell/m64020_190307_140324/unaligned/unmapped",
    "testing/per_flowcell/m64020_190816_194042/unaligned/unmapped",
    "testing/per_flowcell/m64020_190818_014714/unaligned/unmapped"
)

l = parse_args(default_args)

N = length(l[["files"]])
sample_name = l[["sample_name"]]

#m64020_190419_185501/unaligned/unmapped.read_metrics.prl_counts.txt
#m64020_190419_185501/unaligned/unmapped.read_metrics.prl_hist.txt
#m64020_190419_185501/unaligned/unmapped.read_metrics.prl_nx.txt
#m64020_190419_185501/unaligned/unmapped.read_metrics.prl_yield_hist.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users   16 Oct 19 23:55 unmapped.read_metrics.np_hist.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users  285 Oct 19 23:55 unmapped.read_metrics.prl_counts.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users 1.1K Oct 19 23:55 unmapped.read_metrics.prl_hist.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users  960 Oct 19 23:55 unmapped.read_metrics.prl_nx.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users 1.6K Oct 19 23:55 unmapped.read_metrics.prl_yield_hist.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users  231 Oct 19 23:55 unmapped.read_metrics.range_gap_hist.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users  207 Oct 19 23:55 unmapped.read_metrics.rl_counts.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users  881 Oct 19 23:55 unmapped.read_metrics.rl_hist.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users  797 Oct 19 23:55 unmapped.read_metrics.rl_nx.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users 1.3K Oct 19 23:55 unmapped.read_metrics.rl_yield_hist.txt
#-rw-r--r--  1 kiran CHARLES\Domain Users  161 Oct 19 23:55 unmapped.read_metrics.zmw_hist.txt

raw.unmapped.flag_stats.txt
raw.unmapped.read_metrics.np_hist.txt
raw.unmapped.read_metrics.prl_counts.txt
raw.unmapped.read_metrics.prl_hist.txt
raw.unmapped.read_metrics.prl_nx.txt
raw.unmapped.read_metrics.prl_yield_hist.txt
raw.unmapped.read_metrics.range_gap_hist.txt
raw.unmapped.read_metrics.rl_counts.txt
raw.unmapped.read_metrics.rl_hist.txt
raw.unmapped.read_metrics.rl_nx.txt
raw.unmapped.read_metrics.rl_yield_hist.txt
raw.unmapped.read_metrics.zmw_hist.txt

yield_hists = NULL

for (i in 1:N) {
    prefix = l[["files"]][i]
    flowcell_name = l[["flowcell_name"]][i]
    
    prl_counts = paste(prefix, "read_metrics.prl_counts.txt", sep=".")
    prl_hist = paste(prefix, "read_metrics.prl_hist.txt", sep=".")
    prl_nx = paste(prefix, "read_metrics.prl_nx.txt", sep=".")
    prl_yield_hist = paste(prefix, "read_metrics.prl_yield_hist.txt", sep=".")
    
    logger("Load data...")
    f_prl_counts = read_tsv(prl_counts)
    f_prl_hist = read_tsv(prl_hist)
    f_prl_nx = read_tsv(prl_nx)
    f_prl_yield_hist = read_tsv(prl_yield_hist)
    
    logger("Plotting yield...")
    nv = f_prl_nx$value
    
    l.90 = last(which(f_prl_yield_hist$polymerase_read_length <= nv[90]))
    l.70 = last(which(f_prl_yield_hist$polymerase_read_length <= nv[70]))
    l.50 = last(which(f_prl_yield_hist$polymerase_read_length <= nv[50]))
    l.30 = last(which(f_prl_yield_hist$polymerase_read_length <= nv[30]))
    l.10 = last(which(f_prl_yield_hist$polymerase_read_length <= nv[10]))
    l.0 = as.integer(nrow(f_prl_yield_hist))
    
    y.50_30 = 1.02*max(na.omit(f_prl_yield_hist$yield_bp[l.50:l.30]))
    y.90_50 = 1.02*max(na.omit(f_prl_yield_hist$yield_bp[l.90:l.50]))
    y.10_0 = 1.02*max(na.omit(f_prl_yield_hist$yield_bp[l.10:nrow(f_prl_yield_hist)]))
    
    v.90 = round(f_prl_yield_hist[l.90,] %>% pull(polymerase_read_length) / 1000)
    v.50 = round(f_prl_yield_hist[l.50,] %>% pull(polymerase_read_length) / 1000)
    v.10 = round(f_prl_yield_hist[l.10,] %>% pull(polymerase_read_length) / 1000)
    
    for (dev in devs) {
        if (dev == "pdf") {
            pdf(paste(flowcell_name, "polymerase_read_yield.pdf", sep="."), height=8, width=8)
        } else if (dev == "png") {
            png(paste(flowcell_name, "polymerase_read_yield.png", sep="."), height=8, width=8, units="in", res=300)
        }

        par(mar=c(5, 9, 4, 2))
    
        plot(f_prl_yield_hist$polymerase_read_length, f_prl_yield_hist$yield_bp, type="l", lwd=0.1, bty="n", xlab="Read Length (bp)", ylab="", xlim=c(0, 1.2*max(na.omit(f_prl_yield_hist$polymerase_read_length))), ylim=c(0, 1.1*max(na.omit(f_prl_yield_hist$yield_bp))), xaxt="n", yaxt="n", main=flowcell_name)
        axis(1, at=axTicks(1), labels=formatC(axTicks(1), big.mark=','))
        axis(2, at=axTicks(2), labels=formatC(axTicks(2), big.mark=','), las=1)
        mtext("Density\n(Bases Read per Read Length)", 2, line=ceiling(log10(max(axTicks(2)))) - 2)
    
        segments(nv[50], 100, nv[50], y.50_30, lty=3, lwd=1)
        segments(nv[50], y.50_30, nv[30], y.50_30, lty=3, lwd=1)
        arrows(f_prl_yield_hist$polymerase_read_length[l.30] - 1000, y.50_30, f_prl_yield_hist$polymerase_read_length[l.30], y.50_30, lwd=1, length=0.05)
        text(nv[50], y.50_30, labels=paste("50% of bases in reads:\n>", v.50, "kbp\n"), pos=4, cex=0.9, offset=-0.10, adj=c(1, 1))
    
        segments(nv[90], 100, nv[90], y.90_50, lty=3, lwd=1)
        segments(nv[90], y.90_50, nv[70], y.90_50, lty=3, lwd=1)
        arrows(f_prl_yield_hist$polymerase_read_length[l.70] - 1000, y.90_50, f_prl_yield_hist$polymerase_read_length[l.70], y.90_50, lwd=1, length=0.05)
        text(nv[90], y.90_50, labels=paste("90% of bases in reads:\n>", v.90, "kbp\n"), pos=4, cex=0.9, offset=-0.10, adj=c(1, 1))
    
        segments(nv[10], 100, nv[10], y.10_0, lty=3, lwd=1)
        segments(nv[10], y.10_0, f_prl_yield_hist$yield_bp[nrow(f_prl_yield_hist)], y.10_0, lty=3, lwd=1)
        arrows(f_prl_yield_hist$polymerase_read_length[nrow(f_prl_yield_hist)] - 1000, y.10_0, f_prl_yield_hist$polymerase_read_length[nrow(f_prl_yield_hist)], y.10_0, lwd=1, length=0.05)
        text(nv[10], y.10_0, labels=paste("10% of bases in reads:\n>", v.10, "kbp\n"), pos=4, cex=0.9, offset=-0.10, adj=c(1, 1))
    
        polygon(c(0, f_prl_yield_hist$polymerase_read_length[1:l.0], f_prl_yield_hist$polymerase_read_length[l.0] + 1000), c(0, f_prl_yield_hist$yield_bp[1:l.0], 0), col="#334D5C", border=NA)
        polygon(c(0, f_prl_yield_hist$polymerase_read_length[1:(l.10 + 1)]), c(0, f_prl_yield_hist$yield_bp[1:l.10], 0), col="#45B29D", border=NA)
        polygon(c(0, f_prl_yield_hist$polymerase_read_length[1:(l.50 + 1)]), c(0, f_prl_yield_hist$yield_bp[1:l.50], 0), col="#EFC94C", border=NA)
        polygon(c(0, f_prl_yield_hist$polymerase_read_length[1:(l.90 + 1)]), c(0, f_prl_yield_hist$yield_bp[1:l.90], 0), col="#DF5A49", border=NA)
    
        points(f_prl_yield_hist$polymerase_read_length, f_prl_yield_hist$yield_bp, type="l", lwd=0.5)
    
        if (dev != "screen") {
            dev.off()
        }
    }
    
    #f_prl_yield_hist = read_tsv(prl_yield_hist)
    f_prl_yield_hist$flowcell = flowcell_name
    if (is.null(yield_hists)) {
        yield_hists = f_prl_yield_hist
    } else {
        yield_hists = rbind(yield_hists, f_prl_yield_hist)
    }
}

for (dev in devs) {
    if (dev == "pdf") {
        pdf(paste(sample_name, "polymerase_read_yield.pdf", sep="."), height=8, width=8)
    } else if (dev == "png") {
        png(paste(sample_name, "polymerase_read_yield.png", sep="."), height=8, width=8, units="in", res=300)
    }
    
    par(mar=c(5, 9, 4, 2))
    plot(0, 0, type="n", xaxt="n", yaxt="n", xlim=c(0, max(yield_hists$polymerase_read_length)), ylim=c(0, max(yield_hists$yield_bp)), bty="n", xlab="Read Length (bp)", ylab="", main=sample_name)

    for (i in 1:N) {
        flowcell_name = l[["flowcell_name"]][i]
        df = yield_hists %>% filter(flowcell == flowcell_name)
        
        polygon(c(0, df$polymerase_read_length), c(0, df$yield_bp), col=gsub("$", "AA", colors[i]), border=NA)
        points(df$polymerase_read_length, df$yield_bp, type="l", lwd=0.5, col="#000000AA")
    
        axis(1, at=axTicks(1), labels=formatC(axTicks(1), big.mark=','))
        axis(2, at=axTicks(2), labels=formatC(axTicks(2), big.mark=','), las=1)
        mtext("Density\n(Bases Read per Read Length)", 2, line=ceiling(log10(max(axTicks(2)))) - 2)
    }
    
    legend("topright", l[["flowcell_name"]], fill=gsub("$", "AA", colors[1:N]), bty="n", cex=0.8)
    
    if (dev != "screen") {
        dev.off()
    }
}
