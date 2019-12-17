################################################################################
## usage example:
# source("plot.resources.R");
# resources = get_time_series_df("resources.log");
# plot_time_series_df(resources, "test.pdf", "resources usage by task")
################################################################################

get_time_series_df <- function(resources.log.file) {
    library("stringr")

    logs = readLines(resources.log.file)

    b = grep("Runtime Information", logs) + 1
    L = length(logs)
    i = grep("Read/Write IO:", logs[(L-5):L], fixed = T)
    if (2 == length(i)) {# last block is complete
        logs = logs[b:L]
    } else {# ditch last block which is incomplete
        logs = logs[b:(L-(6-i))]
    }
    rm(i, b, L)

    ##########
    loc = Sys.getlocale("LC_TIME")
    Sys.setlocale("LC_TIME", "C")
    # covert to dataframe of time series
    t = gsub("(\\[|\\])", "", sub("UTC ", "", sub("^  ", "", logs[1])))
    start.time = strptime(t, format = "%a %b %d %H:%M:%S %Y", tz = "UTC")

    b.arrays = grep("^  \\[", logs)
    res = vector(mode = "list", length = length(b.arrays))
    j = 1
    for (i in b.arrays) {
        block = logs[i: (i+4)]

        t = gsub("(\\[|\\])", "", sub("UTC ", "", sub("^  ", "", block[1])))
        time = strptime(t, format = "%a %b %d %H:%M:%S %Y", tz = "UTC")
        cpu = as.numeric(sub("%", "", str_match(block[2], "[0-9\\.]+%")[1,1]))
        mem_gb = as.numeric(sub(" GiB", "", str_match(block[3], "[0-9\\.]+ GiB")[1,1]))
        mem_percent = as.numeric(sub("%", "", str_match(block[3], "[0-9\\.]+%")[1,1]))
        disk_gb = as.numeric(sub(" GiB", "", str_match(block[4], "[0-9\\.]+ GiB")[1,1]))
        disk_percent = as.numeric(sub("%", "", str_match(block[4], "[0-9\\.]+%")[1,1]))
        # io = str_match(block[5], "[0-9\\.]+%")[1,1]
        res[[j]] = list(time = time,
                        cpu = cpu,
                        mem_gb = mem_gb,
                        mem_percent = mem_percent,
                        disk_gb = disk_gb,
                        disk_percent = disk_percent)
        j = j+1
    }
    Sys.setlocale("LC_TIME", loc)

    resources = do.call(rbind, lapply(res, data.frame))
    resources$"cpu" = as.numeric(resources$"cpu")
    resources$"mem_gb" = as.numeric(resources$"mem_gb")
    resources$"mem_percent" = as.numeric(resources$"mem_percent")
    resources$"disk_gb" = as.numeric(resources$"disk_gb")
    resources$"disk_percent" = as.numeric(resources$"disk_percent")

    rm(res, i,j)
    resources
}

plot_time_series_df <- function(resources.df, out.pdf, title) {
    library("ggplot2")
    library("reshape2")
    re = melt(resources.df, id = "time")
    p = ggplot(data = re, aes(x=time, y=value)) + geom_line() +
        facet_grid(variable ~ ., scales = "free") +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(p, filename = out.pdf)
}