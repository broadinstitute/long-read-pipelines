
# function to pretty-display cost of runnning a cromwell job retrieved from Big Query
summarize_cost <- function(big.query.json.file,
                           markdown.table.md,
                           plot.pdf) {
    library("jsonlite")
    library("knitr")
    library("dplyr")
    library("ggplot2")
    library("gridExtra")

    json = fromJSON(big.query.json.file)

    json$"cost" = as.numeric(json$"cost")
    json$"sku_description" = gsub(" running in Americas", "", json$"sku_description")

    non.trivial = json[json$"cost" > 0.01, ]
    non.trivial = non.trivial[order(non.trivial$"cost", decreasing = T), ]

    df = non.trivial[, c(4,3,2)]
    presentation = df
    presentation$"labels_value" = sub("cromwell-[0-9a-z]+-[0-9a-z]+-[0-9a-z]+-[0-9a-z]+-[0-9a-z]+,", "",
                                      df$"labels_value")
    names(presentation)= c("unitcost", "sku", "task")
    presentation$"sku" = as.factor(presentation$"sku")
    presentation$"task" = as.factor(presentation$"task")

    total.cost = sum(presentation$"unitcost")

    # save markdown table
    sink(markdown.table.md)
    md.table = kable(presentation, caption = sprintf("Total cost: %.2f", total.cost),
                     format = "markdown", row.names = F)
    print( md.table )
    sink()

    # save bar chart
    # https://stackoverflow.com/questions/47752037/pie-chart-with-ggplot2-with-specific-order-and-percentage-annotations
    per.task <- presentation %>%
        group_by(task) %>%
        summarise(cost = sum(unitcost)) %>%
        mutate(share=cost/sum(cost)*100.0) %>%
        arrange(desc(cost))
    # p = ggplot(per.task, aes("", share, fill = task)) +
    #     geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    #     coord_polar("y") +
    #     geom_text(aes(label = paste0(round(share), "%")),
    #               position = position_stack(vjust = 0.5)) +
    #     labs(x = NULL, y = NULL, fill = NULL,
    #          title = sprintf("Cost share per task class out of total cost: %.2f", total.cost)) +
    #     guides(fill = guide_legend(reverse = TRUE)) +
    #     theme_classic() +
    #     theme(axis.line = element_blank(),
    #           axis.text = element_blank(),
    #           axis.ticks = element_blank(),
    #           plot.title = element_text(hjust = 0.5))
    cost.labeller = vector()
    for(t in per.task$"task") {
        s = sprintf("%s: $%.2f", t, per.task$"cost"[ per.task$"task" == t ])
        cost.labeller <- c(cost.labeller, s)
    }
    names(cost.labeller) = per.task$"task"
    per.sku.per.task = NULL
    for (t in levels(presentation$"task")) {
        tmp <- presentation %>%
            filter(task == t) %>%
            group_by(sku) %>%
            summarise(cost = sum(unitcost)) %>%
            mutate(share=cost/sum(cost)*100.0) %>%
            arrange(desc(share))
        df = as.data.frame(tmp)
        df = subset(df, select=-c(cost))
        df$"sku" = as.factor(df$"sku")
        df$"task" = t
        if (is.null(per.sku.per.task)) {
            per.sku.per.task = df
        } else {
            per.sku.per.task = rbind.data.frame(per.sku.per.task, df)
        }
    }
    per.sku.per.task$"sku" = as.factor(per.sku.per.task$"sku")
    per.sku.per.task$"task" = as.factor(per.sku.per.task$"task")
    q = ggplot(per.sku.per.task, aes("", share, fill = sku)) +
        geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
        coord_polar("y") +
        geom_text(aes(label = paste0(round(share), "%")),
                  position = position_stack(vjust = 0.5)) +
        labs(x = NULL, y = NULL, fill = NULL,
             title = sprintf("Cost per task per SKU. Total cost: $%.2f", total.cost)) +
        guides(fill = guide_legend(reverse = TRUE, ncol = 2)) +
        facet_wrap(facets=. ~ task, ncol=2, labeller = as_labeller(cost.labeller))  +
        # theme_classic() +
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.position="bottom")
        # annotation_custom(tableGrob(mytable), xmin=35, xmax=50, ymin=-2.5, ymax=-1)
    ggsave(plot.pdf, plot = q)
    # ggsave(plot.pdf, arrangeGrob(p, q, nrow = 2))
}
